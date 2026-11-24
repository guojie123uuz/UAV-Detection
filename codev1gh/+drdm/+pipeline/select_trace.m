function selection = select_trace(trace,threshold)
%SELECT_TRACE 对固定搜索记录应用存在性门限。
%
% 输入：trace：完整搜索记录；threshold：不低于流程下限的门限。
% 输出：selection：路径、拟合、区间预算、完成状态及停止原因。
%
% 顺序模式选取相同搜索序列的前缀，联合模式整体接纳或拒绝最终模型。

validateattributes(threshold,{'numeric'},{'scalar','real','nonnegative'});
assert(threshold>=trace.options.minimum_threshold,'drdm:ThresholdFloor', ...
    'The trace was computed only for thresholds at or above its fixed floor.');
if isfield(trace,'mode') && strcmp(trace.mode,'joint_model')
    selection = drdm.pipeline.select_joint_model(trace,threshold);
    return;
end
paths = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
fit = [];
diagnostics = [];
selectedRound = 0;
inspected = 0;
complete = false;
reason = 'budget_limited';
for k = 1:numel(trace.rounds)
    record = trace.rounds(k);
    inspected = k;
    % 加径后的残差平方和非负，因此改进量 T≤J0。
    % 当上界已低于门限时，可直接停止，无须继续候选精化。
    if record.null_cost<=threshold
        complete = true;
        reason = 'below_threshold_bound';
        break;
    end
    if ~record.complete
        reason = 'optimization_or_candidate_budget';
        break;
    end
    if ~strcmp(record.kind,'qualified_candidate')
        complete = true;
        reason = record.kind;
        break;
    end
    if record.score<=threshold
        complete = true;
        reason = 'below_threshold';
        break;
    end
    fit = record.fit;
    paths = fit.paths;
    diagnostics = record.diagnostics;
    selectedRound = k;
end
selection = struct('paths',paths,'fit',fit,'diagnostics',diagnostics, ...
    'complete',complete,'stop_reason',reason,'selected_round',selectedRound,'inspected_rounds',inspected, ...
    'threshold',threshold,'procedure_id',trace.procedure_id, ...
    'reference_reliable',false,'scope','Complete means this bounded decision procedure, not exhaustive physical absence.');
end
