function scan = scan_additional_path(context,nullSeeds)
%SCAN_ADDITIONAL_PATH 对固定背景进行一轮加径比较。
%
% 输入：context：条件模型与固定权重；nullSeeds：原假设的初始路径。
% 输出：scan：最大残差下降量、原/备择拟合、候选与完成状态。
%
% 原假设和每个加径备择均重新估计已有路径。

nullFit = drdm.estimation.fit_paths(context,nullSeeds);
prediction = drdm.model.predict_from_context(context,nullFit.paths);
p = context.observations.rdmp;
y = context.observations.rdmy;
p.complex_values = p.complex_values-prediction.rdmp.complex_values;
y.complex_values = y.complex_values-prediction.rdmy.complex_values;
candidates = drdm.estimation.propose_candidates(p,y,nullFit.paths,context);
order = drdm.estimation.rank_candidates(context,nullFit.paths,candidates);
attempted = min(numel(order),context.options.refine_candidates);
trials = repmat(struct('candidate_index',0,'converged',false,'identifiable',false,'cost',inf),1,attempted);
best = [];
maximum = 0;
for k = 1:attempted
    candidate = candidates(order(k));
    seed = struct('tau_s',candidate.tau_s,'nu_hz',candidate.nu_hz,'gain_complex',0);
    fitted = drdm.estimation.fit_paths(context,[nullFit.paths seed]);
    trials(k) = struct('candidate_index',order(k),'converged',fitted.converged, ...
        'identifiable',fitted.identifiable,'cost',fitted.cost);
    improvement = nullFit.cost-fitted.cost;
    if fitted.converged && fitted.identifiable && improvement>maximum
        best = fitted;
        maximum = improvement;
    end
end
complete = nullFit.converged && nullFit.identifiable && attempted==numel(candidates) && all([trials.converged]);
statistic = maximum;
if ~complete
    statistic = nan;
end % 保留失败状态，供门限训练和评价计数。
id = drdm.statistics.fingerprint(struct('algorithm','conditional-add-one-v1', ...
    'context',drdm.statistics.conditional_signature(context),'null_seeds',nullSeeds));
scan = struct('statistic',statistic,'attempted_maximum',maximum,'complete',complete, ...
    'null_fit',nullFit,'best_fit',best,'candidates',candidates,'trials',trials, ...
    'candidate_count',numel(candidates),'attempted_count',attempted,'procedure_id',id, ...
    'context_id',drdm.statistics.conditional_signature(context), ...
    'scope','One conditional add-one round, including null-path refitting; not an entire sequential search.');
end
