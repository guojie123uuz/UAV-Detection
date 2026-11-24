function draws = simulate_background(context,paths,payloadSeeds,noiseSeeds)
%SIMULATE_BACKGROUND 产生独立的双图背景统计样本。
%
% 输入：context、paths：条件模型与参考路径；payloadSeeds、noiseSeeds：配对随机种子。
% 输出：draws：实观测矩阵、参考路径和种子信息。
%
% 每列由同一 IQ 同时生成 RDMp/RDMy，保留分支间相关性。

drdm.statistics.check_seed_sets(context.conditional_model.payload_seeds,payloadSeeds,noiseSeeds);
assert(numel(payloadSeeds)==numel(noiseSeeds),'drdm:SampleCount','Payload/noise draw counts differ.');
samples = zeros(numel(context.observed),numel(payloadSeeds));
for k = 1:numel(payloadSeeds)
    block = drdm.statistics.simulate_conditional_block(context,paths,payloadSeeds(k),noiseSeeds(k));
    draw = drdm.statistics.rebind_observation(context,block);
    samples(:,k) = draw.observed;
end
draws = struct('samples',samples,'payload_seeds',payloadSeeds(:).','noise_seeds',noiseSeeds(:).', ...
    'reference_paths',paths,'procedure_id',drdm.statistics.conditional_signature(context));
end
