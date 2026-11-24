function [context,pilot] = background_context(block,prepared)
%BACKGROUND_CONTEXT 由当前 CPI 发现初始路径并估计背景协方差。
%
% 输入：block、prepared：接收数据与准备对象。
% 输出：context：固定模型、白化权重、误差样本及预算方向；pilot：初始路径搜索结果。
%
% 背景样本由当前初始路径传播缓存的独立源生成；预算点各保留实、虚两个方向。

cfg = prepared.cfg;
known = prepared.known;
observations = struct('rdmp',drdm.rdm.compute_rdmp_re(block,known,cfg.rdmp), ...
    'rdmy',drdm.rdm.compute_rdmy(block,cfg.rdmy));
pilot = drdm.detection.extract_strong_paths(block,known,observations,cfg);
cfg.estimation = prepared.options.search;
context = drdm.estimation.prepare_context(block,known,cfg,observations,prepared.model);
assert(numel(context.observed)<=prepared.options.background.max_observation_dimension, ...
    'drdm:CovarianceMemory','Choose a smaller fixed observation domain.');
% 在当前 CPI 的初始路径估计处重新估计协方差。
count = size(prepared.moment_sources,2);
samples = zeros(numel(context.observed),count);
source = prepared.model.source;
for draw = 1:count
    source.iq = prepared.moment_sources(:,draw);
    synthetic = block;
    synthetic.iq = prepared.moment_noise(:,draw);
    for k = 1:numel(pilot.fit.paths)
        synthetic.iq = synthetic.iq+drdm.simulation.apply_path(source,pilot.fit.paths(k),cfg.fractional_delay);
    end
    experiment = drdm.statistics.rebind_observation(context,synthetic);
    samples(:,draw) = experiment.observed;
end
bgOptions = prepared.options.background;
draws = struct('samples',samples,'reference_paths',pilot.fit.paths, ...
    'payload_seeds',bgOptions.covariance_payload_seeds,'noise_seeds',bgOptions.covariance_noise_seeds, ...
    'procedure_id',drdm.statistics.conditional_signature(context));
background = drdm.statistics.fit_background(draws,bgOptions);
context = drdm.statistics.attach_background(context,background);
context.retain_jacobian = true;
context.uncertainty_draws = prepared.uncertainty_draws;
prediction = drdm.model.predict_from_context(context,pilot.fit.paths);
context.background_mean_error = background.mean-prediction.normalized_vector;
% 每个复数预算点对应实、虚两个投影方向，随后合并为复数功率。
points = prepared.options.budget_y_points;
indices = zeros(size(points,1),1);
directions = zeros(numel(context.observed),2*size(points,1));
np = 2*numel(context.observations.rdmp.complex_values);
nr = numel(context.y_indices);
ym = context.observations.rdmy;
for k = 1:size(points,1)
    [~,d] = min(abs(ym.delay_s-points(k,1)));
    [~,f] = min(abs(ym.doppler_hz-points(k,2)));
    index = sub2ind(size(ym.complex_values),f,d);
    r = find(context.y_indices==index,1);
    s = find(context.y_imag_indices==index,1);
    assert(~isempty(r) && ~isempty(s),'drdm:BudgetDirections', ...
        'Design points must select retained, non-origin RDMy cells.');
    indices(k) = index;
    directions(np+r,2*k-1) = 1;
    directions(np+nr+s,2*k) = 1;
end
assert(numel(unique(indices))==numel(indices),'drdm:BudgetDirections','Design points map to duplicate cells.');
context.budget_indices = indices;
context.budget_directions = directions;
end
