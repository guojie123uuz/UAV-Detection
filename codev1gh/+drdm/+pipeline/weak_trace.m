function trace = weak_trace(block,prepared,strongThreshold,weakOptions)
%WEAK_TRACE 执行强径筛选并生成一轮弱径搜索记录。
%
% 输入：block、prepared：接收数据及联合模型准备对象；strongThreshold：强径门限；weakOptions：可选弱径配置。
% 输出：trace：强径结果、弱径候选精化、统计量、完成状态及流程标识。
%
% 采用 joint_model 和 independent_linearized，最多五条强径参照；弱径门限留待选择阶段应用。

if nargin<4
    weakOptions = struct();
end
weakOptions = drdm.detection.weak_options(weakOptions);
assert(strcmp(prepared.description.code,drdm.pipeline.code_signature()), ...
    'drdm:CalibrationMismatch','Prepared code changed; rebuild the processing configuration.');
assert(strcmp(prepared.options.reference_selection.mode,'joint_model') && ...
    prepared.options.max_paths<=5 && ~isempty(prepared.uncertainty_draws), ...
    'drdm:WeakPipeline','Use joint-model reference selection, at most five references and independent error draws.');
validateattributes(strongThreshold,{'numeric'},{'scalar','real','finite','nonnegative'});
strongTrace = drdm.pipeline.trace(block,prepared);
strong = drdm.pipeline.select_trace(strongTrace,strongThreshold);
description = struct('version','drdm-weak-search-v1','strong_procedure_id',strongTrace.procedure_id, ...
    'strong_threshold',strongThreshold,'weak_options',weakOptions);
trace = struct('kind','drdm-weak-trace-v1','procedure_id',drdm.statistics.fingerprint(description), ...
    'prepared_id',prepared.procedure_id,'code',prepared.description.code,'description',description, ...
    'strong_trace',strongTrace,'strong',strong,'scan',[], ...
    'complete',strong.complete,'statistic',inf,'reason','strong_search_incomplete');
if ~strong.complete
    return;
end
if isempty(strong.paths)
    trace.statistic = 0;
    trace.reason = 'no_selected_reference';
    return;
end
trace.scan = drdm.detection.scan_weak_path(strongTrace.context,strong.paths,weakOptions);
trace.complete = trace.scan.complete;
trace.statistic = trace.scan.statistic;
trace.reason = trace.scan.reason;
end
