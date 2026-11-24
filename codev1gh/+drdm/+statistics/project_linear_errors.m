function moments = project_linear_errors(jacobian,cholesky,errors,directions)
%PROJECT_LINEAR_ERRORS 传播独立观测误差并计算参数偏差和三明治协方差。
%
% 输入：jacobian：白化实雅可比；cholesky：拟合协方差因子；errors：误差列；directions：实观测预算方向。
% 输出：moments：缩放参数偏差、协方差、先验投影及重拟合残差矩。
%
% 令 J=L\G，参数误差为 J⁺L⁻¹e；偏差校正均值的有限样本方差通过 C(1+1/B) 计入。

count = size(errors,2);
assert(count>=3 && size(errors,1)==size(jacobian,1), ...
    'drdm:UncertaintySamples','Need at least three observation error columns matching the Jacobian.');
[u,s,v] = svd(jacobian,0);
singular = diag(s);
assert(numel(singular)==size(jacobian,2) && ...
    all(singular>max(size(jacobian))*eps(max([singular;1]))), ...
    'drdm:UncertaintyRank','Parameter Jacobian must have full column rank.');

whitened = cholesky\errors;
coefficients = u.'*whitened;
parameterErrors = (v./singular.')*coefficients;
parameterBias = mean(parameterErrors,2);
centred = parameterErrors-parameterBias;
covariance = (centred*centred.')/(count-1);

% 偏差校正的均值由独立有限样本估计，其方差贡献为 C/B。
correctedCovariance = covariance*(1+1/count);
prior = directions.'*(cholesky*u)*coefficients;
residual = directions.'*errors-prior;
moments = struct('parameter_bias',parameterBias,'parameter_covariance',covariance, ...
    'corrected_covariance',correctedCovariance,'sample_count',count, ...
    'parameter_error_samples',parameterErrors,'prior_samples',prior,'residual_samples',residual, ...
    'prior_mse',second_moment(prior),'residual_mse',second_moment(residual), ...
    'method','Independent simulated errors with fixed-weight linear propagation');
end

function value = second_moment(samples)
% 直接对同一实现的平方取平均；重拟合残差已包含拟合与观测的相关性。
value = mean(samples.^2,2);
end
