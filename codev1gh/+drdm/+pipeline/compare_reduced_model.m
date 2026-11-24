function comparison = compare_reduced_model(context,fullFit,history,optimizer)
%COMPARE_REDUCED_MODEL 比较 K 径联合模型与至多 K−1 径模型。
%
% 输入：context、fullFit：上下文和完整拟合；history：搜索记录；optimizer：可选约简优化配置。
% 输出：comparison：约简最小代价、代价差统计量、全部初值结果及完成状态。
%
% 分别删去各径并加入历史 K−1 径模型作为初值，其余路径自由重估。

if nargin<4
    defaults = drdm.pipeline.default_options();
    optimizer = defaults.reference_selection.reduced_optimizer;
end
% 约简模型含未解释成分，较差初值可能沿高残差区域缓慢移动。
% 采用单独配置的差分方式和预算，保留每个初值的实际收敛状态。
reducedContext = context;
reducedContext.options.finite_difference_type = optimizer.finite_difference_type;
reducedContext.options.max_iterations = optimizer.max_iterations;
reducedContext.options.max_function_evaluations = optimizer.max_function_evaluations;
count = numel(fullFit.paths);
assert(count>=1,'drdm:ModelOrder','Reduced-model comparison requires at least one path.');
seeds = cell(1,count);
for k = 1:count
    seeds{k} = fullFit.paths([1:k-1 k+1:count]);
end
for k = 1:numel(history)
    fit = history(k).fit;
    if ~isempty(fit) && numel(fit.paths)==count-1 && ...
            ~any(cellfun(@(s)isequaln(s,fit.paths),seeds))
        seeds{end+1} = fit.paths; %#ok<AGROW>
    end
end
trials = repmat(struct('initial_paths',[],'fit',[]),1,numel(seeds));
% 零增益属于“至多 K-1 径”零假设，空模型提供一个确定可行的代价上界。
empty = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
best = drdm.estimation.fit_paths(context,empty);
for k = 1:numel(seeds)
    fit = drdm.estimation.fit_paths(reducedContext,seeds{k});
    trials(k) = struct('initial_paths',seeds{k},'fit',fit);
    if fit.converged && fit.cost<best.cost
        best = fit;
    end
end
complete = all(arrayfun(@(t)t.fit.converged,trials));
comparison = struct('full_order',count,'null_max_order',count-1, ...
    'reduced_fit',best,'reduced_cost',best.cost,'full_cost',fullFit.cost, ...
    'statistic',max(0,best.cost-fullFit.cost),'complete',complete,'trials',trials, ...
    'optimizer',optimizer);
end
