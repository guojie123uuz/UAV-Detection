function uncertainty = independent_uncertainty(context,fit)
%INDEPENDENT_UNCERTAINTY 在当前拟合点估计独立条件误差的不确定度。
%
% 输入：context：含 uncertainty_draws 和预算方向；fit：路径及白化雅可比。
% 输出：uncertainty：参数偏差、三明治协方差、相关残差预算及种子。
%
% 独立误差样本与条件均值、拟合权重和观测仿真使用不同随机种子。

draws = context.uncertainty_draws;
samples = drdm.statistics.conditional_error_samples(context,fit.paths,draws);
uncertainty = drdm.statistics.project_linear_errors(fit.whitened_jacobian, ...
    context.cholesky,samples.errors,context.budget_directions);
uncertainty.payload_seeds = draws.payload_seeds;
uncertainty.noise_seeds = draws.noise_seeds;
end
