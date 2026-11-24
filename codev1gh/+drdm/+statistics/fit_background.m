function background = fit_background(draws,options)
%FIT_BACKGROUND 估计实堆叠观测的均值与收缩协方差。
%
% 输入：draws：背景样本；options：收缩系数、方差下限及维数限制。
% 输出：background：均值、协方差、Cholesky 因子和稳定性诊断。
%
% 采用 (1-λ)S+λdiag(S)+εI，并报告有效秩及分组均值差。

x = draws.samples;
validateattributes(x,{'numeric'},{'2d','real','finite','nonempty'});
assert(size(x,2)>=3,'drdm:BackgroundSamples','At least three independent draws are required.');
assert(size(x,1)<=options.max_observation_dimension, ...
    'drdm:CovarianceMemory','Reduce the fixed observation domain before fitting a dense covariance.');
validateattributes(options.shrinkage,{'numeric'},{'scalar','real','>=',0,'<=',1});
validateattributes(options.relative_variance_floor,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(options.absolute_variance_floor,{'numeric'},{'scalar','real','finite','positive'});
% 估计实/虚部联合协方差，并向其对角阵收缩以稳定求解。
average = mean(x,2);
centred = x-average;
raw = (centred*centred.')/(size(x,2)-1);
floorValue = max(options.absolute_variance_floor,options.relative_variance_floor*max(diag(raw)));
covariance = (1-options.shrinkage)*raw+options.shrinkage*diag(diag(raw))+floorValue*eye(size(raw));
covariance = (covariance+covariance.')/2;
[factor,flag] = chol(covariance,'lower');
assert(flag==0,'drdm:CovariancePositive','Regularized background covariance is not positive definite.');
half = floor(size(x,2)/2);
splitDifference = mean(x(:,1:half),2)-mean(x(:,half+1:end),2);
background = struct('mean',average,'covariance',covariance,'cholesky',factor, ...
    'sample_count',size(x,2),'dimension',size(x,1),'shrinkage',options.shrinkage, ...
    'variance_floor',floorValue,'raw_effective_rank',trace(raw)^2/max(sum(raw(:).^2),realmin), ...
    'split_mean_whitened_rms',norm(factor\splitDifference)/sqrt(size(x,1)), ...
    'mean_standard_error',sqrt(diag(raw)/size(x,2)), ...
    'reference_paths',draws.reference_paths,'payload_seeds',draws.payload_seeds, ...
    'noise_seeds',draws.noise_seeds,'procedure_id',draws.procedure_id, ...
    'scope','Fixed reference covariance; finite Monte Carlo; conditional on assumed traffic.');
end
