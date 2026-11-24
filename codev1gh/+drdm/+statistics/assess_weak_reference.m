function assessment = assess_weak_reference(context,fit,location,samples,gains,requirements,mode)
%ASSESS_WEAK_REFERENCE 联合评价弱径响应、强径重拟合与误差预算。
%
% 输入：context、fit：固定模型和强径拟合；location：候选位置；samples：独立误差；gains：复增益；requirements：要求；mode：可选 joint/rdmp/rdmy。
% 输出：assessment：response、geometry、budget 及局部可辨识性。
%
% 分支比较各自截取原协方差并重新白化，保留相应观测相关性。

if nargin<7
    mode = 'joint';
end
assert(ismember(mode,{'joint','rdmp','rdmy'}),'drdm:ObservationMode','Use joint, rdmp or rdmy.');
assert(fit.converged && isfield(fit,'whitened_jacobian'), ...
    'drdm:ReferenceFit','Use a converged fit with a retained Jacobian.');
assert(isequaln(fit.paths,samples.reference_paths) && ...
    strcmp(samples.procedure_id,drdm.statistics.conditional_signature(context)), ...
    'drdm:CalibrationMismatch','Error samples belong to different paths or processing coordinates.');
response = drdm.model.weak_path_response(context,fit.paths,location);
np = 2*numel(context.observations.rdmp.complex_values);
rows = 1:numel(context.observed);
if strcmp(mode,'rdmp')
    rows = 1:np;
end
if strcmp(mode,'rdmy')
    rows = np+1:numel(context.observed);
end
% 子集比较须从原协方差重建权重，不能截取联合白化向量的若干行。
l = context.cholesky;
jacobian = fit.whitened_jacobian;
if ~isempty(l)
    jacobian = l*jacobian;
    covariance = l(rows,:)*l(rows,:).';
    l = chol(covariance,'lower');
end
geometry = drdm.statistics.profile_weak_response(response.linear(rows,:),response.self(rows),jacobian(rows,:),l);
budget = drdm.statistics.weak_reference_budget(geometry,samples.errors(rows,:),gains,requirements);
budget.budget_satisfied = budget.budget_satisfied & fit.identifiable;
assessment = struct('mode',mode,'observation_rows',rows,'reference_paths',fit.paths, ...
    'location',response.location,'response',response,'geometry',geometry,'budget',budget, ...
    'reference_fit_identifiable',fit.identifiable, ...
    'usable_local_model',fit.identifiable && geometry.gain_estimable, ...
    'status','local_weak_reference_diagnostic','procedure_id',samples.procedure_id);
end
