function context = attach_background(context,background)
%ATTACH_BACKGROUND 将估计背景的白化权重装入上下文。
%
% 输入：context：当前观测模型；background：匹配的背景统计。
% 输出：context：更新权重、白化目标和初始残差代价。
%
% 同轮全部候选共享相同协方差及归一化尺度。

assert(strcmp(background.procedure_id,drdm.statistics.conditional_signature(context)), ...
    'drdm:CalibrationMismatch','Background belongs to a different observation/model configuration.');
context.cholesky = background.cholesky;
np = 2*numel(context.observations.rdmp.complex_values);
context.p_cholesky = chol(background.covariance(1:np,1:np),'lower');
context.background = background;
% 当前 IQ 和双图保持不变，只更新白化目标及零路径模型的代价。
context.target = drdm.estimation.whiten(context.observed,context.cholesky);
emptyPaths = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
prediction = drdm.model.predict_from_context(context,emptyPaths);
residual = drdm.estimation.whiten(prediction.normalized_vector-context.observed,context.cholesky);
context.initial_cost = sum(residual.^2);
end
