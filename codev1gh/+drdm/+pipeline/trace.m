function trace = trace(block,prepared)
%TRACE 生成与存在性门限无关的完整参照搜索记录。
%
% 输入：block、prepared：当前接收数据及准备对象。
% 输出：trace：上下文、初始搜索、逐轮记录或联合模型阶次比较。
%
% 按 reference_selection.mode 分派顺序筛选与联合模型筛选；两种模式均重新估计当前 CPI 的背景。

block = drdm.frontend.validate_block(block);
options = prepared.options;
[context,pilot] = drdm.pipeline.background_context(block,prepared);
if strcmp(options.reference_selection.mode,'joint_model')
    trace = drdm.pipeline.joint_trace(context,pilot,prepared);
    return;
end
paths = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
current = drdm.estimation.fit_paths(context,paths);
roundTemplate = struct('null_cost',0,'score',0,'complete',false,'kind','', ...
    'fit',[],'diagnostics',[],'candidates',struct([]),'trials',struct([]),'chosen_index',0);
rounds = repmat(roundTemplate,1,options.max_paths);
roundCount = 0;
endReason = 'budget_limited';
% 各轮沿同一候选选择规则推进，最终门限只截取该序列的前缀。
for iteration = 1:options.max_paths
    if current.cost<=options.minimum_threshold
        roundCount = roundCount+1;
        record = roundTemplate;
        record.null_cost = current.cost;
        record.complete = true;
        record.kind = 'below_procedure_floor';
        rounds(roundCount) = record;
        endReason = record.kind;
        break;
    end
    prediction = drdm.model.predict_from_context(context,current.paths);
    p = context.observations.rdmp;
    y = context.observations.rdmy;
    p.complex_values = p.complex_values-prediction.rdmp.complex_values;
    y.complex_values = y.complex_values-prediction.rdmy.complex_values;
    candidates = drdm.estimation.propose_candidates(p,y,current.paths,context);
    record = roundTemplate;
    record.null_cost = current.cost;
    record.candidates = candidates;
    if isempty(candidates)
        record.complete = true;
        record.kind = 'no_candidates';
    else
        [record.fit,record.diagnostics,record.score,record.trials,record.chosen_index] = ...
            fit_candidates(context,current,candidates,options);
        record.complete = numel(record.trials)==numel(candidates) && all([record.trials.converged]);
        if ~record.complete
            record.kind = 'optimization_or_candidate_budget';
        elseif isempty(record.fit)
            record.kind = 'no_qualified_reference';
        else
            record.kind = 'qualified_candidate';
        end
    end
    roundCount = roundCount+1;
    rounds(roundCount) = record;
    if ~strcmp(record.kind,'qualified_candidate')
        endReason = record.kind;
        break;
    end
    current = record.fit;
end
id = drdm.statistics.fingerprint(struct('procedure',prepared.procedure_id, ...
    'fs_hz',block.fs_hz,'samples',numel(block.iq),'valid_mask',block.valid_mask));
trace = struct('version','drdm-sequential-trace-v1','procedure_id',id, ...
    'prepared_id',prepared.procedure_id,'options',options,'context',context,'pilot',pilot, ...
    'rounds',rounds(1:roundCount),'end_reason',endReason, ...
    'scope','Frozen pilot/background/quality rules; threshold only terminates this same prefix.');
end

function [best,diagnostics,score,trials,bestIndex] = fit_candidates(context,current,candidates,options)
[order,~] = drdm.estimation.rank_candidates(context,current.paths,candidates);
attempts = min(numel(order),context.options.refine_candidates);
template = struct('candidate_index',0,'converged',false,'fit',[],'diagnostics',[],'reason','');
trials = repmat(template,1,attempts);
best = [];
diagnostics = [];
score = 0;
bestIndex = 0;
for k = 1:attempts
    candidateIndex = order(k);
    fit = drdm.estimation.fit_paths(context,[current.paths seed(candidates(candidateIndex))]);
    check = drdm.pipeline.local_diagnostics(context,fit,options);
    if drdm.estimation.has_close_paths(fit.paths,context)
        check.eligible = false;
        check.reason = 'unresolved_cluster';
    end
    trials(k) = struct('candidate_index',candidateIndex,'converged',fit.converged, ...
        'fit',fit,'diagnostics',check,'reason',check.reason);
    improvement = current.cost-fit.cost;
    if check.eligible && improvement>score
        best = fit;
        diagnostics = check;
        score = improvement;
        bestIndex = candidateIndex;
    end
end
end

function path = seed(candidate)
path = struct('tau_s',candidate.tau_s,'nu_hz',candidate.nu_hz,'gain_complex',0);
end
