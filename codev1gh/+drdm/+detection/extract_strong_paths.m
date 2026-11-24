function search = extract_strong_paths(block,known,observations,cfg)
%EXTRACT_STRONG_PATHS 按残差下降准则逐径搜索并联合更新。
%
% 输入：block、known、observations、cfg：接收块、模板、双图和搜索配置。
% 输出：search：路径、迭代历史、候选、复数预测及残差。
%
% 用于初始路径发现；每轮重估全部已有路径。

context = drdm.estimation.prepare_context(block,known,cfg,observations);
options = context.options;
paths = struct('tau_s',{},'nu_hz',{},'gain_complex',{});
fit = drdm.estimation.fit_paths(context,paths);
historyTemplate = struct('iteration',0,'cost_before',0,'cost_after',0,'improvement',0, ...
    'relative_improvement',0,'proposal_source','','candidate_count',0,'trials',struct([]));
history = repmat(historyTemplate,1,options.max_paths+1);
historyCount = 0;
pending = struct([]);
costHistory = nan(1,options.max_paths+1);
costHistory(1) = fit.cost;
costCount = 1;
reason = '';
iteration = 0;
while true
    prediction = drdm.model.predict_observations(known,block,paths,cfg,options.fit_mode);
    pResidual = residual_map(observations.rdmp,prediction.rdmp);
    yForSearch = observations.rdmy;
    if context.joint
        yForSearch = residual_map(observations.rdmy,prediction.rdmy);
    end
    if fit.cost<=options.residual_fraction_stop*max(context.initial_cost,realmin)
        pending = struct([]);
        reason = 'numerical_residual_floor';
        break;
    end
    pending = drdm.estimation.propose_candidates(pResidual,yForSearch,paths,context);
    if numel(paths)>=options.max_paths
        reason = 'budget_limited';
        break;
    end
    if isempty(pending)
        reason = 'no_candidates';
        break;
    end
    iteration = iteration+1;
    % 在同一观测域内重新求解复增益，按加权残差对候选排序。
    order = drdm.estimation.rank_candidates(context,paths,pending);
    best = [];
    bestIndex = 0;
    anyConverged = false;
    anyIdentifiable = false;
    trialCount = min(options.refine_candidates,numel(order));
    trialTemplate = struct('candidate_index',0,'cost',0,'exitflag',0,'identifiable',false, ...
        'reason','','paths',paths,'jacobian_condition',inf);
    trials = repmat(trialTemplate,1,trialCount);
    for trial = 1:trialCount
        k = order(trial);
        candidateFit = drdm.estimation.fit_paths(context,[paths seed_path(pending(k))]);
        anyConverged = anyConverged || candidateFit.converged;
        anyIdentifiable = anyIdentifiable || (candidateFit.converged && candidateFit.identifiable);
        improvement = fit.cost-candidateFit.cost;
        close = drdm.estimation.has_close_paths(candidateFit.paths,context);
        if ~candidateFit.converged
            trialReason = 'optimization_failed';
        elseif ~candidateFit.identifiable || close
            trialReason = 'unresolved_cluster';
        elseif improvement<=0
            trialReason = 'no_improvement';
        elseif improvement/max(context.initial_cost,realmin)<options.min_total_improvement || ...
                improvement/max(fit.cost,realmin)<options.min_relative_improvement
            trialReason = 'below_development_threshold';
        else
            trialReason = 'eligible';
        end
        trials(trial) = struct('candidate_index',k,'cost',candidateFit.cost, ...
            'exitflag',candidateFit.exitflag,'identifiable',candidateFit.identifiable,'reason',trialReason, ...
            'paths',candidateFit.paths,'jacobian_condition',candidateFit.jacobian_condition);
        if strcmp(trialReason,'eligible') && (isempty(best) || candidateFit.cost<best.cost)
            best = candidateFit;
            bestIndex = k;
        end
    end
    before = fit.cost;
    if isempty(best)
        if trialCount<numel(order)
            reason = 'candidate_budget_limited';
        elseif ~anyConverged
            reason = 'optimization_failed';
        elseif ~anyIdentifiable || all(strcmp({trials.reason},'unresolved_cluster'))
            reason = 'unresolved_cluster';
        else
            reason = 'below_development_threshold';
        end
        historyCount = historyCount+1;
        history(historyCount) = struct('iteration',iteration,'cost_before',before,'cost_after',before, ...
            'improvement',0,'relative_improvement',0,'proposal_source','none', ...
            'candidate_count',numel(pending),'trials',trials);
        break;
    end
    fit = best;
    paths = fit.paths;
    historyCount = historyCount+1;
    history(historyCount) = struct('iteration',iteration,'cost_before',before,'cost_after',fit.cost, ...
        'improvement',before-fit.cost,'relative_improvement',(before-fit.cost)/max(before,realmin), ...
        'proposal_source',pending(bestIndex).source,'candidate_count',numel(pending),'trials',trials);
    costCount = costCount+1;
    costHistory(costCount) = fit.cost;
end
history = history(1:historyCount);
costHistory = costHistory(1:costCount);
% 由原始观测和当前全部路径重新计算预测与残差，避免逐次相减积累误差。
prediction = drdm.model.predict_observations(known,block,paths,cfg,options.fit_mode);
pResidual = residual_map(observations.rdmp,prediction.rdmp);
yResidual = [];
if context.joint
    yResidual = residual_map(observations.rdmy,prediction.rdmy);
end
records = paths;
for k = 1:numel(records)
    records(k).role = 'development_estimate';
    records(k).reference_reliable = false;
    records(k).time_center_s = block.time_start_s+(numel(block.iq)-1)/(2*block.fs_hz);
    records(k).reference_id = block.reference_id;
end
search = struct('status','uncalibrated','fit_mode',options.fit_mode,'estimated_paths',records, ...
    'fit',fit,'prediction',prediction,'rdmp_residual',pResidual,'rdmy_residual',yResidual, ...
    'pending_candidates',pending,'history',history,'cost_history',costHistory,'stop_reason',reason, ...
    'candidate_budget_per_round',options.refine_candidates, ...
    'statistic_kind','fixed-weight residual reduction; not a calibrated GLRT', ...
    'options',options,'resolution_scales',[context.delay_scale_s context.doppler_scale_hz]);
end

function path = seed_path(candidate)
path = struct('tau_s',candidate.tau_s,'nu_hz',candidate.nu_hz,'gain_complex',candidate.gain_complex);
end

function map = residual_map(original,prediction)
map = original;
map.complex_values = original.complex_values-prediction.complex_values;
if map.normalization_energy>0
    map.normalized_values = map.complex_values/map.normalization_energy;
else
    map.normalized_values = complex(nan(size(map.complex_values)));
end
map.metadata.residual_source = 'original complex observation minus current complete model';
end
