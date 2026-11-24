function diagnostics = local_diagnostics(context,fit,options)
%LOCAL_DIAGNOSTICS 计算联合参数区间和相关重拟合残差预算。
%
% 输入：context、fit：固定权重下的联合拟合；options：不确定度方法、区间系数及精度要求。
% 输出：diagnostics：拟合点、区间中心、协方差、先验比例、信号误差比及判定。
%
% 独立误差模式校正区间中心偏差，精度约束按区间相对于实际拟合点的最大距离计算。

count = numel(fit.paths);
q = drdm.estimation.path_vector(context,fit.paths);
parameters = 4*count;
if ~fit.identifiable || ~fit.converged || ~isfield(fit,'whitened_jacobian')
    diagnostics = struct('eligible',false,'reason','unresolved_or_failed','interval',nan(parameters,2), ...
        'unit_halfwidth',inf(parameters,1),'physical_centre',nan(parameters,1), ...
        'prior_mse_fraction',inf,'signal_to_total_error',0,'covariance',nan(parameters));
    return;
end
scales = [repmat([context.delay_scale_s;context.doppler_scale_hz],count,1);ones(2*count,1)];
estimateCentre = q.*scales;
if strcmp(options.uncertainty.method,'independent_linearized')
    % 在当前拟合处传播独立误差，区分拟合权重与实际误差协方差。
    moments = drdm.statistics.independent_uncertainty(context,fit);
    covariance = moments.corrected_covariance;
    parameterBias = moments.parameter_bias.*scales;
    centre = estimateCentre-parameterBias;
    priorMse = sum(reshape(moments.prior_mse,2,[]),1).';
    residualMse = sum(reshape(moments.residual_mse,2,[]),1).';
    sampleCount = moments.sample_count;
else
    % 信息矩阵近似保留作对照，观测误差协方差在此近似中取拟合权重。
    [u,s,v] = svd(fit.whitened_jacobian,0);
    singular = diag(s);
    covariance = (v./singular.')*(v./singular.').';
    centre = estimateCentre;
    parameterBias = zeros(parameters,1);
    projectedNoise = context.budget_directions.'*context.cholesky;
    projectedPrior = projectedNoise*u;
    bias = context.cholesky\context.background_mean_error;
    priorBias = projectedPrior*(u.'*bias);
    totalBias = projectedNoise*bias-priorBias;
    priorVariance = sum(projectedPrior.^2,2);
    residualVariance = max(0,sum(projectedNoise.^2,2)-priorVariance);
    priorMse = sum(reshape(priorVariance+priorBias.^2,2,[]),1).';
    residualMse = sum(reshape(residualVariance+totalBias.^2,2,[]),1).';
    sampleCount = 0;
end
% 缩放参数协方差恢复为物理单位，再应用训练得到的区间系数。
unitHalfwidth = sqrt(max(0,diag(covariance))).*scales;
halfwidth = options.interval_multiplier*unitHalfwidth;
interval = [centre-halfwidth centre+halfwidth];
% 偏差校正引起的区间移动也计入实际参照点的精度要求。
extent = max(abs(interval-estimateCentre),[],2);
amplitudes = options.budget_signal_amplitudes(:);
if isscalar(amplitudes)
    amplitudes = repmat(amplitudes,numel(priorMse),1);
end
priorFraction = options.interval_multiplier^2*priorMse./amplitudes.^2;
ratio = amplitudes.^2./max(options.interval_multiplier^2*residualMse,realmin);
precision = all(extent(1:2:2*count)<=options.max_delay_halfwidth_s) && ...
    all(extent(2:2:2*count)<=options.max_doppler_halfwidth_hz);
budget = all(priorFraction<=options.max_prior_mse_fraction) && all(ratio>=options.min_signal_to_total_error);
boundary = any(fit.at_boundary);
reason = 'eligible';
if boundary
    reason = 'boundary';
elseif ~precision
    reason = 'prior_uncertain';
elseif ~budget
    reason = 'cancellation_budget_exceeded';
end
diagnostics = struct('eligible',~boundary && precision && budget,'reason',reason, ...
    'interval',interval,'unit_halfwidth',unitHalfwidth, ...
    'physical_centre',centre,'prior_mse_fraction',priorFraction,'signal_to_total_error',ratio, ...
    'estimate_centre',estimateCentre,'parameter_bias',parameterBias,'uncertainty_extent',extent, ...
    'covariance',covariance.*(scales*scales.'),'precision_pass',precision,'budget_pass',budget, ...
    'uncertainty_sample_count',sampleCount,'method',options.uncertainty.method);
end
