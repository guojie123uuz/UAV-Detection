function updated = rebind_observation(context,block)
%REBIND_OBSERVATION 在固定估计上下文中替换接收数据。
%
% 输入：context：模型、权重及坐标；block：同域新数据块。
% 输出：updated：重算双图、实观测向量及初始代价的上下文。
%
% 采样率、样本数和有效掩码保持一致。

block = drdm.frontend.validate_block(block);
assert(block.fs_hz==context.block.fs_hz && numel(block.iq)==numel(context.block.iq) && ...
    isequal(block.valid_mask,context.block.valid_mask), ...
    'drdm:StatisticalDomain','Cannot change sampling/mask within a statistical experiment.');
updated = context;
updated.block = block;
updated.observations = struct('rdmp',drdm.rdm.compute_rdmp_re(block,context.known,context.cfg.rdmp), ...
    'rdmy',drdm.rdm.compute_rdmy(block,context.cfg.rdmy));
updated.observed = drdm.estimation.pack_observations(updated, ...
    updated.observations.rdmp.complex_values/context.p_scale, ...
    updated.observations.rdmy.complex_values/context.y_scale);
updated.target = drdm.estimation.whiten(updated.observed,updated.cholesky);
empty = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
prediction = drdm.model.predict_from_context(updated,empty);
updated.initial_cost = sum(drdm.estimation.whiten(prediction.normalized_vector-updated.observed,updated.cholesky).^2);
end
