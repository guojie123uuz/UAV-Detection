function search = fit_joint_model(context,options)
%FIT_JOINT_MODEL 逐次扩展暂定路径集合以形成联合模型。
%
% 输入：context：固定观测与权重；options：完整流程配置。
% 输出：search：最终拟合、各阶搜索记录、完成状态及整体残差判据。
%
% 整体残差达到 residual_cost_limit 后进入参数精度与模型阶次检验。

paths = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
emptyFit = drdm.estimation.fit_paths(context,paths);
current = emptyFit;
limit = options.reference_selection.residual_cost_limit;
history = struct('fit',{},'candidates',{},'trials',{},'chosen_index',{});
complete = true;
reason = 'model_residual_satisfied';
% 逐次扩展模型，每次重估全部已有路径；各候选共享固定白化权重。
while current.cost>limit
    if numel(current.paths)>=options.max_paths
        complete = false;
        reason = 'model_path_budget';
        break;
    end
    prediction = drdm.model.predict_from_context(context,current.paths);
    p = context.observations.rdmp;
    y = context.observations.rdmy;
    p.complex_values = p.complex_values-prediction.rdmp.complex_values;
    y.complex_values = y.complex_values-prediction.rdmy.complex_values;
    candidates = drdm.estimation.propose_candidates(p,y,current.paths,context);
    [order,~] = drdm.estimation.rank_candidates(context,current.paths,candidates);
    attempts = min(numel(order),context.options.refine_candidates);
    trials = repmat(struct('candidate_index',0,'fit',[],'numerically_eligible',false),1,attempts);
    best = [];
    bestIndex = 0;
    for k = 1:attempts
        index = order(k);
        candidate = candidates(index);
        seed = struct('tau_s',candidate.tau_s,'nu_hz',candidate.nu_hz,'gain_complex',0);
        fit = drdm.estimation.fit_paths(context,[current.paths seed]);
        eligible = fit.converged && fit.identifiable && ...
            ~drdm.estimation.has_close_paths(fit.paths,context) && fit.cost<current.cost;
        trials(k) = struct('candidate_index',index,'fit',fit,'numerically_eligible',eligible);
        if eligible && (isempty(best) || fit.cost<best.cost)
            best = fit;
            bestIndex = index;
        end
    end
    history(end+1) = struct('fit',best,'candidates',candidates,'trials',trials,'chosen_index',bestIndex); %#ok<AGROW>
    allConverged = all(arrayfun(@(t)t.fit.converged,trials));
    if attempts<numel(candidates) || ~allConverged
        complete = false;
        reason = 'optimization_or_candidate_budget';
        break;
    end
    if isempty(best)
        complete = false;
        reason = 'unresolved_model_residual';
        break;
    end
    current = best;
end
search = struct('fit',current,'empty_fit',emptyFit,'history',history, ...
    'complete',complete,'reason',reason,'model_adequate',current.cost<=limit, ...
    'residual_cost_limit',limit);
end
