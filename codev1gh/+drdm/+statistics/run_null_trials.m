function trials = run_null_trials(context,nullPaths,payloadSeeds,noiseSeeds)
%RUN_NULL_TRIALS 在给定原假设下重复候选产生和加径拟合。
%
% 输入：nullPaths：原假设路径；context：条件背景；payloadSeeds、noiseSeeds：抽样种子。
% 输出：trials：各次统计量、完成状态、候选数量和处理标识。
%
% 样本包含候选选择及原/备择拟合，用于单轮门限训练。

reserved = [];
if isfield(context,'reserved_analysis_seeds')
    reserved = context.reserved_analysis_seeds;
end
drdm.statistics.check_seed_sets(context.conditional_model.payload_seeds, ...
    context.background.payload_seeds,context.background.noise_seeds,reserved,payloadSeeds,noiseSeeds);
assert(numel(payloadSeeds)==numel(noiseSeeds),'drdm:SampleCount','Null draw counts differ.');
statistics = nan(1,numel(payloadSeeds));
complete = false(size(statistics));
candidateCounts = zeros(size(statistics));
procedureId = '';
for k = 1:numel(payloadSeeds)
    block = drdm.statistics.simulate_conditional_block(context,nullPaths,payloadSeeds(k),noiseSeeds(k));
    experiment = drdm.statistics.rebind_observation(context,block);
    scan = drdm.detection.scan_additional_path(experiment,nullPaths);
    statistics(k) = scan.statistic;
    complete(k) = scan.complete;
    candidateCounts(k) = scan.candidate_count;
    procedureId = scan.procedure_id;
end
trials = struct('statistics',statistics,'complete',complete,'candidate_counts',candidateCounts, ...
    'procedure_id',procedureId,'context_id',drdm.statistics.conditional_signature(context), ...
    'payload_seeds',payloadSeeds(:).','noise_seeds',noiseSeeds(:).', ...
    'scope','Conditional single add-one search round at declared H0 parameters and fixed background weights.');
end
