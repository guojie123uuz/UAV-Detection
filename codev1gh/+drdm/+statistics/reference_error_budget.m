function budget = reference_error_budget(context,fit,uncertainty,directions,signalAmplitudes)
%REFERENCE_ERROR_BUDGET 传播路径先验误差并计算重拟合残差预算。
%
% 输入：context、fit、uncertainty；directions：实观测投影；signalAmplitudes：投影设计幅度。
% 输出：budget：先验均方误差比例、总残差及设计信号与误差功率比。
%
% 先验传播项与重拟合后的总误差分别统计，保留二者的相关关系。

assert(strcmp(uncertainty.procedure_id,drdm.statistics.conditional_signature(context)), ...
    'drdm:CalibrationMismatch','Bootstrap and prediction configurations differ.');
assert(isequaln(uncertainty.reference_paths,fit.paths), ...
    'drdm:CalibrationMismatch','Bootstrap belongs to a different fitted path set.');
validateattributes(directions,{'numeric'},{'2d','real','finite','nonempty'});
assert(size(directions,1)==numel(context.observed) && all(vecnorm(directions)>0), ...
    'drdm:BudgetDirections','Directions must use the fixed normalized real observation coordinates.');
validateattributes(signalAmplitudes,{'numeric'},{'vector','real','finite','positive'});
if isscalar(signalAmplitudes)
    signalAmplitudes = repmat(signalAmplitudes,1,size(directions,2));
end
assert(numel(signalAmplitudes)==size(directions,2),'drdm:BudgetSignal','One projected target amplitude per direction.');
% 对缩放参数做中心差分，得到固定观测模型的局部敏感度。
[q,lower,upper] = drdm.estimation.path_vector(context,fit.paths);
jacobian = zeros(numel(context.observed),numel(q));
for k = 1:numel(q)
    step = context.options.finite_difference_step*max(1,abs(q(k)));
    left = q;
    right = q;
    left(k) = max(lower(k),q(k)-step);
    right(k) = min(upper(k),q(k)+step);
    minus = drdm.model.predict_from_context(context,drdm.estimation.paths_from_vector(context,left));
    plus = drdm.model.predict_from_context(context,drdm.estimation.paths_from_vector(context,right));
    jacobian(:,k) = (plus.normalized_vector-minus.normalized_vector)/(right(k)-left(k));
end
projected = directions.'*jacobian;
priorVariance = max(0,diag(projected*uncertainty.scaled_covariance*projected.'));
priorBias = projected*uncertainty.scaled_bias;
% 总残差直接来自每次重拟合后的观测，先验传播误差另行报告。
residuals = directions.'*uncertainty.refitted_residual_samples(:,uncertainty.valid);
totalBias = mean(residuals,2);
totalMse = mean(residuals.^2,2);
signalPower = signalAmplitudes(:).^2;
budget = struct('procedure_id',uncertainty.procedure_id, ...
    'reference_paths',fit.paths, ...
    'valid',all(uncertainty.valid) && uncertainty.tail_resolved && all(isfinite(totalMse)), ...
    'direction_count',size(directions,2),'signal_amplitudes',signalAmplitudes(:), ...
    'prior_variance',priorVariance,'prior_bias',priorBias, ...
    'prior_mse_fraction',(priorVariance+priorBias.^2)./signalPower, ...
    'total_residual_bias',totalBias,'total_residual_mse',totalMse, ...
    'signal_to_total_error',signalPower./max(totalMse,realmin), ...
    'scaled_model_jacobian',jacobian,'directions',directions, ...
    'scope','Local linear prior diagnostic; total error measured AFTER refitting, no prior covariance added again.');
end
