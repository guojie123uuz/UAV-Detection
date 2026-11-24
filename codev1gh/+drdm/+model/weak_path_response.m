function response = weak_path_response(context,references,location)
%WEAK_PATH_RESPONSE 构造候选弱径相对于强径模型的条件均值增量。
%
% 输入：context：联合上下文；references：强径；location：候选 tau_s 与 nu_hz。
% 输出：response：线性矩阵 linear、自项 self、各参照贡献及处理标识。
%
% 弱径复增益为 a 时，增量为 H[Re(a);Im(a)]+q|a|²，包含全部强弱和弱强有序径对。

assert(context.joint && isfield(context,'conditional_model'), ...
    'drdm:ConditionalModel','Use a joint conditional observation context.');
validateattributes(location.tau_s,{'numeric'},{'scalar','real','finite'});
validateattributes(location.nu_hz,{'numeric'},{'scalar','real','finite'});
candidate = struct('tau_s',location.tau_s,'nu_hz',location.nu_hz,'gain_complex',1);
paths = [references candidate];
coordinates = [[paths.tau_s]/context.delay_scale_s;[paths.nu_hz]/context.doppler_scale_hz];
assert(all(all(coordinates>=context.lower & coordinates<=context.upper)), ...
    'drdm:PathBounds','Reference and candidate coordinates must lie in the model domain.');
[dictionary,kernels] = drdm.model.conditional_kernels(context,coordinates);
count = numel(references);
weak = count+1;
np = size(dictionary,1);
ny = prod(kernels.map_size);
zeroP = complex(zeros(np,1));
zeroY = complex(zeros(ny,1));
direct = [drdm.estimation.pack_observations(context,dictionary(:,weak),zeroY), ...
    drdm.estimation.pack_observations(context,1i*dictionary(:,weak),zeroY)];
perReference = zeros(numel(context.observed),2,count);
for k = 1:count
    forward = kernels.pairs(:,weak,k)*conj(references(k).gain_complex);
    reverse = kernels.pairs(:,k,weak)*references(k).gain_complex;
    % a*u+conj(a)*v = real(a)*(u+v)+imag(a)*1i*(u-v)。
    perReference(:,1,k) = drdm.estimation.pack_observations(context,zeroP,forward+reverse);
    perReference(:,2,k) = drdm.estimation.pack_observations(context,zeroP,1i*(forward-reverse));
end
response = struct('linear',direct+sum(perReference,3), ...
    'self',drdm.estimation.pack_observations(context,zeroP,kernels.pairs(:,weak,weak)), ...
    'rdmp_linear',direct,'reference_linear',perReference, ...
    'reference_paths',references,'location',rmfield(candidate,'gain_complex'), ...
    'procedure_id',drdm.statistics.conditional_signature(context), ...
    'model','Conditional finite-ensemble mean; exact quadratic dependence on weak gain.');
end
