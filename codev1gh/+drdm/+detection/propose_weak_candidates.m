function candidates = propose_weak_candidates(context,nullFit,options)
%PROPOSE_WEAK_CANDIDATES 合并残差峰、RDMy 双向差分解释与固定网格候选。
%
% 输入：context：联合观测；nullFit：强径拟合；options：完整弱径搜索配置。
% 输出：candidates：时延、多普勒、来源、参照编号及初始分数。
%
% 复数双图减去强径均值响应后生成候选，网格在固定搜索域内按相同尺度去重。

% 先在复数双图中扣除强径模型，再产生残差峰和双向差分候选。
prediction = drdm.model.predict_from_context(context,nullFit.paths);
p = context.observations.rdmp;
y = context.observations.rdmy;
p.complex_values = p.complex_values-prediction.rdmp.complex_values;
y.complex_values = y.complex_values-prediction.rdmy.complex_values;
proposalContext = context;
proposalContext.options.candidates_per_branch = options.candidates_per_branch;
proposalContext.options.candidate_merge_radius = options.candidate_merge_radius;
candidates = drdm.estimation.propose_candidates(p,y,nullFit.paths,proposalContext);
delay = options.grid_delay_s;
frequency = options.grid_doppler_hz;
if isempty(delay)
    delay = p.delay_s;
end
if isempty(frequency)
    frequency = p.doppler_hz;
end
assert(numel(delay)*numel(frequency)<=options.max_grid_points, ...
    'drdm:WeakGridBudget','Declared weak grid exceeds max_grid_points.');
% 固定网格与图上候选按统一的分辨尺度合并。
scale = [context.delay_scale_s;context.doppler_scale_hz];
existing = [[nullFit.paths.tau_s];[nullFit.paths.nu_hz]]./scale;
chosen = zeros(2,0);
if ~isempty(candidates)
    chosen = [[candidates.tau_s];[candidates.nu_hz]]./scale;
end
for f = frequency(:).'
    for t = delay(:).'
        point = [t;f]./scale;
        assert(all(point>=context.lower & point<=context.upper), ...
            'drdm:PathBounds','Declared weak grid lies outside the fixed observation domain.');
        if any(vecnorm(existing-point)<context.options.duplicate_radius)
            continue;
        end
        if any(vecnorm(chosen-point)<options.candidate_merge_radius)
            continue;
        end
        candidates(end+1) = struct('tau_s',t,'nu_hz',f,'gain_complex',0, ...
            'source','parameter_grid','reference_index',0,'pair_sign',0,'seed_score',0); %#ok<AGROW>
        chosen(:,end+1) = point; %#ok<AGROW>
    end
end
end
