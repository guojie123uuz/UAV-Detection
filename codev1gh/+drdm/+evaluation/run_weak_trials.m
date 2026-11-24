function trials = run_weak_trials(prepared,strongThreshold,weakOptions,profile,seeds,tolerances)
%RUN_WEAK_TRIALS 在声明路径总体上重复完整强弱径搜索。
%
% 输入：prepared、strongThreshold、weakOptions：处理设置；profile：随机路径总体；seeds：每行参数/业务/噪声种子；tolerances：秒/Hz 对应容差。
% 输出：trials：逐例统计量、完成状态、强径错误、搜索记录及真值。
%
% 用于零假设标定时 profile 仅含强径；未完成搜索和错误强径以 Inf 统计量计入。

validateattributes(seeds,{'numeric'},{'2d','nonempty','ncols',3});
validateattributes(tolerances,{'numeric'},{'vector','numel',2,'real','finite','positive'});
drdm.statistics.check_seed_sets(prepared.internal_seeds,seeds(:));
records = repmat(struct('seeds',[],'statistic',0,'complete',false, ...
    'strong_false',false,'strong_missed_count',0,'reason','','trace',[],'truth',[]),1,size(seeds,1));
procedureId = '';
for k = 1:size(seeds,1)
    [block,truth] = drdm.evaluation.simulate_case(prepared,profile,seeds(k,:));
    trace = drdm.pipeline.weak_trace(block,prepared,strongThreshold,weakOptions);
    assignment = drdm.evaluation.assign_paths(trace.strong.paths,truth.paths,tolerances(1),tolerances(2));
    score = trace.statistic;
    % 求解失败及错误强径无法由提高弱径门限纠正，作为独立控制事件保留。
    if ~trace.complete || assignment.false_count>0
        score = inf;
    end
    if isempty(procedureId)
        procedureId = trace.procedure_id;
    end
    assert(strcmp(procedureId,trace.procedure_id),'drdm:CalibrationMismatch','Trial processing changed.');
    records(k) = struct('seeds',seeds(k,:),'statistic',score,'complete',trace.complete, ...
        'strong_false',assignment.false_count>0,'strong_missed_count',assignment.missed_count, ...
        'reason',trace.reason,'trace',trace,'truth',truth);
end
trials = struct('kind','drdm-weak-trials-v1','procedure_id',procedureId, ...
    'prepared_id',prepared.procedure_id,'code',prepared.description.code, ...
    'population',profile,'tolerances',tolerances(:).','seeds',seeds,'records',records);
trials.id = drdm.statistics.fingerprint(rmfield(trials,'records'));
end
