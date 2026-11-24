function selection = select_joint_model(trace,threshold)
%SELECT_JOINT_MODEL 对最终联合模型整体应用阶次检验门限。
%
% 输入：trace：joint_trace 结果；threshold：不低于流程下限的门限。
% 输出：selection：接纳路径、诊断、完成状态及停止原因。
%
% 统计量严格大于门限时接纳；空模型代价同时提供可直接停止的上界。

validateattributes(threshold,{'numeric'},{'scalar','real','nonnegative'});
assert(threshold>=trace.options.minimum_threshold,'drdm:ThresholdFloor', ...
    'Threshold is below the configured procedure floor.');
record = trace.rounds;
paths = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
fit = [];
diagnostics = [];
selectedRound = 0;
complete = record.complete;
reason = record.kind;
if record.null_cost<=threshold
    complete = true;
    reason = 'below_threshold_bound';
elseif record.complete && strcmp(record.kind,'qualified_candidate')
    if record.score>threshold
        paths = record.fit.paths;
        fit = record.fit;
        diagnostics = record.diagnostics;
        selectedRound = 1;
        reason = 'joint_model_accepted';
    else
        reason = 'below_model_order_threshold';
    end
end
selection = struct('paths',paths,'fit',fit,'diagnostics',diagnostics, ...
    'complete',complete,'stop_reason',reason,'selected_round',selectedRound,'inspected_rounds',1, ...
    'threshold',threshold,'procedure_id',trace.procedure_id,'reference_reliable',false, ...
    'scope','The complete joint model is accepted or rejected as a group after fixed model and quality checks.');
end
