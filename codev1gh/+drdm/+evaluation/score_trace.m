function score = score_trace(trace,truth,threshold,tolerances)
%SCORE_TRACE 评价完整搜索的路径对应、共同覆盖与残差预算。
%
% 输入：trace、truth：已完成搜索与真值；threshold：门限；tolerances：秒/Hz 对应容差。
% 输出：score：错误事件、恢复及覆盖判定、预算、参数误差和停止原因。
%
% point_error 与 parameter_error 分别为真值减拟合点、真值减区间中心；统计保留未完成搜索。

selected = drdm.pipeline.select_trace(trace,threshold);
assignment = drdm.evaluation.assign_paths(selected.paths,truth, ...
    tolerances.delay_s,tolerances.doppler_hz);
maxFalseScore = 0;
intervalScore = 0;
traceIntervalScore = inf;
if isempty(truth)
    traceIntervalScore = 0;
end
% 扫描全部轮次，计入新增错误路径、已有路径漂移和未完成搜索。
for k = 1:numel(trace.rounds)
    record = trace.rounds(k);
    if ~record.complete
        maxFalseScore = max(maxFalseScore,record.null_cost);
    elseif strcmp(record.kind,'qualified_candidate')
        match = drdm.evaluation.assign_paths(record.fit.paths,truth,tolerances.delay_s,tolerances.doppler_hz);
        if match.false_count>0
            maxFalseScore = max(maxFalseScore,record.score);
        end
        if match.false_count==0 && match.missed_count==0 && ~isempty(truth)
            ordered = truth(match.truth_index);
            values = [[ordered.tau_s];[ordered.nu_hz]];
            values = [values(:);real([ordered.gain_complex]).';imag([ordered.gain_complex]).'];
            traceIntervalScore = max(abs(values-record.diagnostics.physical_centre)./ ...
                max(record.diagnostics.unit_halfwidth,realmin));
        end
    end
end
if strcmp(trace.end_reason,'budget_limited') && ~isempty(trace.rounds)
    maxFalseScore = max(maxFalseScore,trace.rounds(end).score);
end
% 共同覆盖要求路径完整恢复且全部实参数同时落在区间内。
recovered = selected.complete && assignment.false_count==0 && assignment.missed_count==0;
covered = recovered;
budgetPassed = recovered;
budgetRatio = inf;
parameterError = [];
pointError = [];
unitHalfwidth = [];
if ~isempty(truth)
    if recovered && ~isempty(selected.paths)
        matchedTruth = truth(assignment.truth_index);
        values = [[matchedTruth.tau_s];[matchedTruth.nu_hz]];
        values = [values(:);real([matchedTruth.gain_complex]).';imag([matchedTruth.gain_complex]).'];
        check = selected.diagnostics;
        parameterError = values-check.physical_centre;
        fitted = [[selected.paths.tau_s];[selected.paths.nu_hz]];
        fitted = [fitted(:);real([selected.paths.gain_complex]).';imag([selected.paths.gain_complex]).'];
        pointError = values-fitted;
        unitHalfwidth = check.unit_halfwidth;
        intervalScore = max(abs(parameterError)./max(check.unit_halfwidth,realmin));
        covered = all(values>=check.interval(:,1) & values<=check.interval(:,2));

        prediction = drdm.model.predict_from_context(trace.context,selected.paths);
        residual = trace.context.observed-prediction.normalized_vector;
        projected = trace.context.budget_directions.'*residual;
        power = sum(reshape(projected.^2,2,[]),1).';
        amplitude = trace.options.budget_signal_amplitudes(:);
        if isscalar(amplitude)
            amplitude = repmat(amplitude,numel(power),1);
        end
        budgetRatio = max(power./amplitude.^2);
        budgetPassed = budgetRatio<=1/trace.options.min_signal_to_total_error;
    else
        intervalScore = inf;
        covered = false;
        budgetPassed = false;
    end
end
score = struct('procedure_id',trace.procedure_id,'max_false_score',maxFalseScore, ...
    'trace_interval_score',traceIntervalScore, ...
    'interval_score',intervalScore,'false_path_event',~selected.complete || assignment.false_count>0, ...
    'complete',selected.complete,'joint_recovery',recovered, ...
    'joint_interval_coverage',covered,'budget_passed',budgetPassed,'budget_ratio',budgetRatio, ...
    'truth_count',numel(truth),'estimated_count',numel(selected.paths), ...
    'false_count',assignment.false_count,'missed_count',assignment.missed_count, ...
    'stop_reason',selected.stop_reason,'parameter_error',parameterError, ...
    'point_error',pointError,'unit_halfwidth',unitHalfwidth, ...
    'selected_paths',selected.paths,'interval',[],'threshold',threshold);
if ~isempty(selected.diagnostics)
    score.interval = selected.diagnostics.interval;
end
end
