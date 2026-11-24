function result = select_weak_trace(trace,threshold)
%SELECT_WEAK_TRACE 对固定弱径搜索应用门限并生成复数残差。
%
% 输入：trace：weak_trace 结果；threshold：非负弱径门限。
% 输出：result：联合路径、弱径、更新参照、接纳标志、预测与双图残差。
%
% 接纳后使用强弱径共同重估结果，否则沿用强径阶段的路径集合。

validateattributes(threshold,{'numeric'},{'scalar','real','nonnegative'});
assert(strcmp(trace.kind,'drdm-weak-trace-v1'),'drdm:WeakTrace','Expected a complete weak-search trace.');
accepted = trace.complete && ~isempty(trace.scan) && ...
    ~isempty(trace.scan.best_fit) && trace.statistic>threshold;
paths = trace.strong.paths;
weak = [];
references = paths;
reason = trace.reason;
if accepted
    paths = trace.scan.best_fit.paths;
    weak = trace.scan.weak_path;
    references = trace.scan.updated_references;
    reason = 'weak_path_accepted';
elseif trace.complete && ~isempty(trace.scan) && ~isempty(trace.scan.best_fit)
    reason = 'below_weak_threshold';
end
context = trace.strong_trace.context;
prediction = drdm.model.predict_from_context(context,paths);
p = context.observations.rdmp;
y = context.observations.rdmy;
p.complex_values = p.complex_values-prediction.rdmp.complex_values;
p.normalized_values = p.complex_values/context.p_scale;
y.complex_values = y.complex_values-prediction.rdmy.complex_values;
y.normalized_values = y.complex_values/context.y_scale;
y.normalization_energy = context.y_scale;
result = struct('status','weak_search_uncalibrated','paths',paths,'weak_path',weak, ...
    'references',references,'weak_accepted',accepted,'complete',trace.complete, ...
    'stop_reason',reason,'statistic',trace.statistic,'threshold',threshold, ...
    'prediction',prediction,'rdmp_residual',p,'rdmy_residual',y, ...
    'procedure_id',trace.procedure_id,'trace',trace,'calibration',[]);
end
