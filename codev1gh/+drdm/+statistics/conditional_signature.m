function value = conditional_signature(context)
%CONDITIONAL_SIGNATURE 生成固定条件观测与拟合配置的标识。
%
% 输入：context：条件模型、观察域、权重和归一化信息。
% 输出：value：配置标识。
%
% 标识由处理参数与样本坐标构成，观测 IQ 数值独立保存。

assert(isfield(context,'conditional_model'),'drdm:ConditionalModel','Expected conditional context.');
description = struct('version','drdm-conditional-fit-v1','source_model',context.conditional_model.fingerprint, ...
    'cfg',context.cfg,'estimation',context.options,'mask',context.block.valid_mask, ...
    'fs_hz',context.block.fs_hz,'sample_count',numel(context.block.iq), ...
    'p_scale',context.p_scale,'y_scale',context.y_scale, ...
    'y_indices',context.y_indices,'y_imag_indices',context.y_imag_indices, ...
    'cholesky',context.cholesky);
value = drdm.statistics.fingerprint(description);
end
