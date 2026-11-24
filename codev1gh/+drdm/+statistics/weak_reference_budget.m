function budget = weak_reference_budget(geometry,errors,gains,requirements)
%WEAK_REFERENCE_BUDGET 按候选响应计算相关误差预算及信号保留率。
%
% 输入：geometry：局部投影结果；errors：独立实误差列；gains：待评价复增益；requirements：信号误差比和先验比例要求。
% 输出：budget：预算判定、方向均方误差、交叉矩及增益误差统计。
%
% 平方矩包含均值偏差；观测误差与强径拟合误差以同一实现配对统计。

validateattributes(errors,{'numeric'},{'2d','real','finite'});
validateattributes(gains,{'numeric'},{'vector','finite','nonempty'});
assert(size(errors,1)==size(geometry.whitened_linear,1) && size(errors,2)>=3, ...
    'drdm:UncertaintySamples','Need at least three errors in the fixed observation coordinates.');
assert(all(abs(gains)>0),'drdm:WeakGain','Candidate gains must be nonzero.');
validateattributes(requirements.min_signal_to_error,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(requirements.max_prior_fraction,{'numeric'},{'scalar','real','finite','positive'});
% 同一误差实现同时产生观测项、强径拟合项和残差，保留三者的相关关系。
white = drdm.estimation.whiten(errors,geometry.cholesky);
u = geometry.nuisance_basis;
prior = u*(u.'*white);
residual = white-prior;
fixed = geometry.fixed_basis;
profiled = geometry.profiled_basis;
observedFixed = fixed.'*white;
priorFixed = fixed.'*prior;
residualFixed = observedFixed-priorFixed;
count = size(errors,2);
observationMse = observedFixed*observedFixed.'/count;
priorMse = priorFixed*priorFixed.'/count;
crossMoment = observedFixed*priorFixed.'/count;
residualMse = residualFixed*residualFixed.'/count;
profiledSamples = profiled.'*residual;
profiledMse = profiledSamples*profiledSamples.'/count;

% 信号同时包含强弱径线性交叉项与弱径自项，按候选相位计算方向误差。
gains = gains(:).';
coefficients = [real(gains);imag(gains)];
signal = geometry.whitened_linear*coefficients+geometry.whitened_self*abs(gains).^2;
effective = signal-u*(u.'*signal);
fixedSignal = fixed.'*signal;
effectiveSignal = profiled.'*effective;
signalPower = sum(effectiveSignal.^2,1);
errorPower = zeros(size(gains));
priorFraction = inf(size(gains));
for k = 1:numel(gains)
    if signalPower(k)>0
        direction = effectiveSignal(:,k)/sqrt(signalPower(k));
        errorPower(k) = direction.'*profiledMse*direction;
    end
    fixedPower = sum(fixedSignal(:,k).^2);
    if fixedPower>0
        direction = fixedSignal(:,k)/sqrt(fixedPower);
        priorFraction(k) = max(0,direction.'*priorMse*direction)/fixedPower;
    end
end
ratio = signalPower./max(errorPower,realmin);
retained = sum(effective.^2,1)./max(sum(signal.^2,1),realmin);
% 在两个增益分量均可估时，单独报告增益误差及弱径自项引起的局部偏差。
gainMean = [];
gainCovariance = [];
selfBias = [];
if geometry.gain_estimable
    gainSamples = geometry.gain_operator*residual;
    gainMean = mean(gainSamples,2);
    centred = gainSamples-gainMean;
    gainCovariance = centred*centred.'/(count-1);
    selfBias = (geometry.gain_operator*geometry.whitened_self)*abs(gains).^2;
end
budget = struct('gains',gains,'sample_count',count,'requirements',requirements, ...
    'gain_estimable',geometry.gain_estimable,'profiled_signal_power',signalPower, ...
    'directional_error_mse',errorPower,'signal_to_error',ratio, ...
    'prior_mse_fraction',priorFraction,'retained_signal_fraction',retained, ...
    'budget_satisfied',geometry.gain_estimable & ratio>=requirements.min_signal_to_error & ...
    priorFraction<=requirements.max_prior_fraction, ...
    'fixed_observation_mse',observationMse,'fixed_prior_mse',priorMse, ...
    'observation_prior_cross_moment',crossMoment,'fixed_residual_mse',residualMse, ...
    'profiled_error_mean',mean(profiledSamples,2),'profiled_error_mse',profiledMse, ...
    'gain_error_mean',gainMean,'gain_error_covariance',gainCovariance, ...
    'weak_self_gain_bias',selfBias, ...
    'scope','Conditional local error budget; finite-sample moments, not a calibrated detection decision.');
end
