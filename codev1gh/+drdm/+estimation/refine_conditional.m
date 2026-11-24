function analysis = refine_conditional(block,cellInfo,seeds,cfg,options)
%REFINE_CONDITIONAL 完成给定路径数下的条件双图精化与重采样。
%
% 输入：block、cellInfo、seeds、cfg；options：可选条件统计配置。
% 输出：analysis：联合拟合、背景、bootstrap 区间、预测和复数残差。
%
% 均值、协方差及 bootstrap 使用分开的业务/噪声随机种子。

if nargin<5
    options = drdm.statistics.default_options();
end
defaults = drdm.statistics.default_options();
assert(isstruct(options) && isscalar(options) && all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'drdm:BackgroundOptions','Unknown conditional estimation option.');
names = fieldnames(options);
for k = 1:numel(names)
    defaults.(names{k}) = options.(names{k});
end
options = defaults;
drdm.statistics.check_seed_sets(options.mean_payload_seeds,options.covariance_payload_seeds, ...
    options.covariance_noise_seeds,options.bootstrap_payload_seeds,options.bootstrap_noise_seeds);
assert(numel(options.covariance_payload_seeds)==numel(options.covariance_noise_seeds) && ...
    numel(options.bootstrap_payload_seeds)==numel(options.bootstrap_noise_seeds), ...
    'drdm:SampleCount','Paired payload/noise sample counts differ.');
assert(isempty(cfg.estimation.observation_covariance), ...
    'drdm:BackgroundOptions','Conditional refinement estimates its covariance from the declared independent draws.');
known = drdm.lte.build_known_resources(cellInfo,cfg.known_resources);
model = drdm.model.build_conditional_model(cellInfo,known,options.model, ...
    options.mean_payload_seeds,options.noise_power,options.max_ensemble_bytes);
jointCfg = cfg;
jointCfg.estimation.fit_mode = 'joint_conditional';
context = drdm.estimation.prepare_context(block,known,jointCfg,[],model);
assert(numel(context.observed)<=options.max_observation_dimension, ...
    'drdm:CovarianceMemory','Select a smaller fixed map before conditional fitting.');
% 先在初始路径处估计背景，再以固定权重精化全部路径。
draws = drdm.statistics.simulate_background(context,seeds, ...
    options.covariance_payload_seeds,options.covariance_noise_seeds);
background = drdm.statistics.fit_background(draws,options);
context = drdm.statistics.attach_background(context,background);
referencePrediction = drdm.model.predict_from_context(context,seeds);
meanDifference = background.mean-referencePrediction.normalized_vector;
meanCheck = struct('independent_mean_difference',meanDifference, ...
    'whitened_rms',norm(context.cholesky\meanDifference)/sqrt(numel(meanDifference)), ...
    'mean_standard_error',background.mean_standard_error,'stability_validated',false, ...
    'scope','Finite ensemble/independent moment discrepancy; increase draws and check load mismatch separately.');
fit = drdm.estimation.fit_paths(context,seeds);
prediction = drdm.model.predict_from_context(context,fit.paths);
% 通过数值与秩检查后，对固定路径集合进行参数 bootstrap。
uncertainty = [];
if fit.converged && fit.identifiable && ~isempty(fit.paths)
    uncertainty = drdm.statistics.bootstrap_refit(context,fit, ...
        options.bootstrap_payload_seeds,options.bootstrap_noise_seeds,options);
end
context.reserved_analysis_seeds = [options.bootstrap_payload_seeds(:);options.bootstrap_noise_seeds(:)];
pResidual = context.observations.rdmp;
yResidual = context.observations.rdmy;
pResidual.complex_values = pResidual.complex_values-prediction.rdmp.complex_values;
pResidual.normalized_values = pResidual.complex_values/context.p_scale;
yResidual.complex_values = yResidual.complex_values-prediction.rdmy.complex_values;
yResidual.normalized_values = yResidual.complex_values/context.y_scale;
yResidual.normalization_energy = context.y_scale;
analysis = struct('status','conditional_refinement_unvalidated','fit',fit,'background',background, ...
    'uncertainty',uncertainty,'mean_check',meanCheck,'prediction',prediction,'rdmp_residual',pResidual, ...
    'rdmy_residual',yResidual,'reference_reliable',false,'context',context,'options',options, ...
    'scope','Fixed selected path count; conditional RMC mean/covariance; no full-search reliability certification.');
end
