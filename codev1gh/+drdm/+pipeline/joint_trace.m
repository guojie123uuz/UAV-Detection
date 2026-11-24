function trace = joint_trace(context,pilot,prepared)
%JOINT_TRACE 形成最终联合模型及多初值约简比较记录。
%
% 输入：context、pilot、prepared：背景上下文、初始搜索及准备对象。
% 输出：trace：最终模型、质量判据、阶次比较及中间搜索历史。
%
% 低阶模型用于生成候选；只有最终满足质量要求的联合模型参与门限选择。

options = prepared.options;
search = drdm.pipeline.fit_joint_model(context,options);
current = search.fit;
comparison = [];
comparisons = {};
diagnostics = [];
record = struct('null_cost',search.empty_fit.cost,'score',0,'complete',search.complete, ...
    'kind',search.reason,'fit',[],'diagnostics',[]);
if search.complete && ~isempty(current.paths)
    % 多初值约简若找到同样满足残差限值的更低阶模型，则继续向下比较。
    while ~isempty(current.paths)
        comparison = drdm.pipeline.compare_reduced_model(context,current,search.history, ...
            options.reference_selection.reduced_optimizer);
        comparisons{end+1} = comparison; %#ok<AGROW>
        if ~comparison.complete
            record.complete = false;
            record.kind = 'reduced_model_optimization_failed';
            break;
        end
        reduced = comparison.reduced_fit;
        if reduced.cost>search.residual_cost_limit
            break;
        end
        current = reduced;
    end
    if record.complete && ~isempty(current.paths)
        diagnostics = drdm.pipeline.local_diagnostics(context,current,options);
        diagnostics.model_adequate = current.cost<=search.residual_cost_limit;
        diagnostics.model_residual_ratio = current.cost/search.residual_cost_limit;
        if drdm.estimation.has_close_paths(current.paths,context)
            diagnostics.eligible = false;
            diagnostics.reason = 'unresolved_cluster';
        end
        record.kind = 'joint_model_quality_rejected';
        if diagnostics.eligible
            record.kind = 'qualified_candidate';
            record.fit = current;
            record.diagnostics = diagnostics;
            record.score = comparison.statistic;
        end
    end
end
if record.complete && isempty(current.paths)
    record.kind = 'adequate_empty_model';
end
block = context.block;
id = drdm.statistics.fingerprint(struct('procedure',prepared.procedure_id, ...
    'fs_hz',block.fs_hz,'samples',numel(block.iq),'valid_mask',block.valid_mask));
trace = struct('version','drdm-joint-model-trace-v1','mode','joint_model', ...
    'procedure_id',id,'prepared_id',prepared.procedure_id,'options',options, ...
    'context',context,'pilot',pilot,'rounds',record,'end_reason',record.kind, ...
    'search_history',search.history,'model_search',search,'model_fit',current, ...
    'model_diagnostics',diagnostics,'order_comparisons',{comparisons}, ...
    'scope','Provisional models are retained for search; only the final joint model is eligible for reference selection.');
end
