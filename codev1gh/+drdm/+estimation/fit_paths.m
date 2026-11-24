function fit = fit_paths(context,seeds,initialGains)
%FIT_PATHS 在固定双图观测上联合精化连续路径参数。
%
% 输入：context：模型与观测；seeds：路径初值；initialGains：联合模式下可选的复增益初值。
% 输出：fit：路径、加权代价、收敛状态、可辨识性及可选白化雅可比。
%
% 仅用 RDMp 时采用变量投影法求解复增益；联合模式优化全部参数并保留径对乘积。

assert(exist('lsqnonlin','file')==2,'drdm:OptimizationToolbox','Optimization Toolbox is required.');
assert(isstruct(seeds) && all(isfield(seeds,{'tau_s','nu_hz','gain_complex'})), ...
    'drdm:Paths','Seeds need tau_s, nu_hz and gain_complex.');
count = numel(seeds);
options = context.options;
coordinates = zeros(2,count);
if count>0
    coordinates = [[seeds.tau_s]/context.delay_scale_s;[seeds.nu_hz]/context.doppler_scale_hz];
end
assert(all(isfinite(coordinates(:))),'drdm:PathSeed','Nonfinite path seed.');
assert(all(all(coordinates>=context.lower & coordinates<=context.upper)), ...
    'drdm:PathBounds','Path seed is outside the fixed search domain.');
% 同一坐标下缓存传播响应，供目标函数和增益求解复用。
cachedCoordinates = [];
cachedDictionary = [];
cachedWaves = [];
[dictionary,initialKernels] = get_dictionary(coordinates);
pData = context.observations.rdmp.complex_values(:)/context.p_scale;
[gains,dictionaryDiagnostics] = drdm.estimation.solve_complex_gains(dictionary,pData,context.p_cholesky);
% 可选联合模型初值用于新增弱径：保留已有拟合及多参照线性增益初值。
% 原双参数调用继续采用 RDMp 最小二乘初始化。
if nargin>=3
    validateattributes(initialGains,{'numeric'},{'vector','finite','numel',count});
    assert(context.joint,'drdm:InitialGains','Explicit gain seeds require a joint fit.');
    gains = initialGains(:);
end
if count==0
    initialModel = drdm.estimation.model_vector(context,dictionary,initialKernels,gains);
    emptyResidual = drdm.estimation.whiten(initialModel-context.observed,context.cholesky);
    fit = struct('paths',seeds,'cost',sum(emptyResidual.^2),'residual',emptyResidual, ...
        'converged',true,'identifiable',true,'exitflag',1,'iterations',0, ...
        'dictionary',dictionaryDiagnostics,'jacobian_rank',0,'jacobian_condition',1, ...
        'at_boundary',false(0,1),'diagnostics_kind','local sensitivity only; uncalibrated');
    return;
end
lower = repmat(context.lower,count,1);
upper = repmat(context.upper,count,1);
q0 = coordinates(:);
if context.joint
    q0 = [q0;real(gains);imag(gains)];
    lower = [lower;-inf(2*count,1)];
    upper = [upper;inf(2*count,1)];
end
% 仅 RDMp 拟合时，复增益在每次坐标迭代中解析求解；双图模式联合优化。
solverOptions = optimoptions('lsqnonlin','Display','off', ...
    'MaxIterations',options.max_iterations,'MaxFunctionEvaluations',options.max_function_evaluations, ...
    'FunctionTolerance',options.function_tolerance,'StepTolerance',options.step_tolerance, ...
    'OptimalityTolerance',options.optimality_tolerance, ...
    'FiniteDifferenceStepSize',options.finite_difference_step,'TypicalX',ones(size(q0)));
solverOptions.FiniteDifferenceType = options.finite_difference_type;
[q,~,~,exitflag,output] = lsqnonlin(@objective,q0,lower,upper,solverOptions);
coordinates = reshape(q(1:2*count),2,count);
[dictionary,waves] = get_dictionary(coordinates);
[gains,dictionaryDiagnostics] = drdm.estimation.solve_complex_gains(dictionary,pData,context.p_cholesky);
if context.joint
    gains = q(2*count+(1:count))+1i*q(3*count+(1:count));
end
vector = drdm.estimation.model_vector(context,dictionary,waves,gains);
residual = drdm.estimation.whiten(vector-context.observed,context.cholesky);
% 用完整参数的数值雅可比检查局部可辨识性，并支持后续误差传播。
fullParameters = [coordinates(:);real(gains);imag(gains)];
fullLower = [repmat(context.lower,count,1);-inf(2*count,1)];
fullUpper = [repmat(context.upper,count,1);inf(2*count,1)];
jacobian = zeros(numel(residual),4*count);
for k = 1:4*count
    step = options.finite_difference_step*max(1,abs(fullParameters(k)));
    left = fullParameters;
    right = fullParameters;
    left(k) = max(fullLower(k),left(k)-step);
    right(k) = min(fullUpper(k),right(k)+step);
    jacobian(:,k) = (full_residual(right)-full_residual(left))/(right(k)-left(k));
end
singular = svd(jacobian,0);
rankJacobian = sum(singular>max(size(jacobian))*eps(max([singular;1])));
condition = inf;
if rankJacobian==4*count
    condition = singular(1)/singular(end);
end
paths = repmat(struct('tau_s',0,'nu_hz',0,'gain_complex',0),1,count);
for k = 1:count
    paths(k) = struct('tau_s',coordinates(1,k)*context.delay_scale_s, ...
        'nu_hz',coordinates(2,k)*context.doppler_scale_hz,'gain_complex',gains(k));
end
boundary = any(coordinates-context.lower<1e-3 | context.upper-coordinates<1e-3,1);
fit = struct('paths',paths,'cost',sum(residual.^2),'residual',residual, ...
    'converged',exitflag>0,'identifiable',dictionaryDiagnostics.rank==2*count && ...
    dictionaryDiagnostics.condition<=options.max_dictionary_condition && ...
    rankJacobian==4*count && condition<=options.max_jacobian_condition, ...
    'exitflag',exitflag,'iterations',output.iterations,'dictionary',dictionaryDiagnostics, ...
    'jacobian_rank',rankJacobian,'jacobian_condition',condition, ...
    'at_boundary',boundary,'diagnostics_kind','local sensitivity only; uncalibrated');
if isfield(context,'retain_jacobian') && context.retain_jacobian
    fit.whitened_jacobian = jacobian;
end

    function [a,x] = get_dictionary(c)
        if ~isequal(c,cachedCoordinates) || isempty(cachedDictionary)
            [cachedDictionary,cachedWaves] = drdm.estimation.unit_dictionary(context,c);
            cachedCoordinates = c;
        end
        a = cachedDictionary;
        x = cachedWaves;
    end

    function value = objective(parameters)
        c = reshape(parameters(1:2*count),2,count);
        [a,x] = get_dictionary(c);
        if context.joint
            gain = parameters(2*count+(1:count))+1i*parameters(3*count+(1:count));
        else
            gain = drdm.estimation.solve_complex_gains(a,pData,context.p_cholesky);
        end
        predicted = drdm.estimation.model_vector(context,a,x,gain);
        value = drdm.estimation.whiten(predicted-context.observed,context.cholesky);
    end

    function value = full_residual(parameters)
        c = reshape(parameters(1:2*count),2,count);
        [a,x] = get_dictionary(c);
        gain = parameters(2*count+(1:count))+1i*parameters(3*count+(1:count));
        predicted = drdm.estimation.model_vector(context,a,x,gain);
        value = drdm.estimation.whiten(predicted-context.observed,context.cholesky);
    end
end
