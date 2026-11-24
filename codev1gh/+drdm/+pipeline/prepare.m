function prepared = prepare(cellInfo,cfg,options)
%PREPARE 准备确认资源、条件模型及独立仿真样本。
%
% 输入：cellInfo、cfg、options：小区、处理参数及完整流程配置。
% 输出：prepared：资源、模型、背景与误差样本缓存、种子及源码配置标识。
%
% 缓存波形和噪声供同配置 CPI 复用；每个 CPI 在其初始路径处重新估计背景。

defaults = drdm.pipeline.default_options();
assert(isstruct(options) && isscalar(options) && isequal(sort(fieldnames(options)),sort(fieldnames(defaults))), ...
    'drdm:PipelineOptions','Use a complete pipeline.default_options structure.');
validateattributes(options.max_paths,{'numeric'},{'scalar','integer','>=',1,'<=',6});
for name = {'interval_multiplier','max_delay_halfwidth_s','max_doppler_halfwidth_hz', ...
        'max_prior_mse_fraction','min_signal_to_total_error'}
    validateattributes(options.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(options.development_threshold,{'numeric'},{'scalar','real','nonnegative'});
validateattributes(options.minimum_threshold,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(options.budget_y_points,{'numeric'},{'2d','real','finite','nonempty','ncols',2});
validateattributes(options.budget_signal_amplitudes,{'numeric'},{'vector','real','finite','positive'});
assert(isscalar(options.budget_signal_amplitudes) || ...
    numel(options.budget_signal_amplitudes)==size(options.budget_y_points,1), ...
    'drdm:BudgetSignal','Supply one signal amplitude per design point or a scalar.');
search = drdm.estimation.check_options(options.search);
assert(strcmp(search.fit_mode,'joint_conditional') && isempty(search.observation_covariance), ...
    'drdm:PipelineOptions','Sequential search uses its own conditional background weights.');
bg = options.background;
selection = options.reference_selection;
assert(ismember(selection.mode,{'sequential','joint_model'}), ...
    'drdm:ReferenceSelection','Use sequential or joint_model reference selection.');
validateattributes(selection.residual_cost_limit,{'numeric'},{'scalar','real','finite','positive'});
reduced = selection.reduced_optimizer;
assert(ismember(reduced.finite_difference_type,{'forward','central'}), ...
    'drdm:FiniteDifference','Use forward or central finite differences.');
validateattributes(reduced.max_iterations,{'numeric'},{'scalar','integer','finite','positive'});
validateattributes(reduced.max_function_evaluations,{'numeric'},{'scalar','integer','finite','positive'});
u = options.uncertainty;
assert(ismember(u.method,{'information','independent_linearized'}), ...
    'drdm:UncertaintyMethod','Use information or independent_linearized uncertainty.');
internalSeeds = [bg.mean_payload_seeds(:);bg.covariance_payload_seeds(:);bg.covariance_noise_seeds(:)];
uncertaintyCount = 0;
if strcmp(u.method,'independent_linearized')
    assert(numel(u.payload_seeds)==numel(u.noise_seeds) && numel(u.payload_seeds)>=3, ...
        'drdm:SampleCount','Need at least three independent uncertainty draws.');
    internalSeeds = [internalSeeds;u.payload_seeds(:);u.noise_seeds(:)];
    uncertaintyCount = numel(u.payload_seeds);
end
drdm.statistics.check_seed_sets(internalSeeds);
assert(numel(bg.covariance_payload_seeds)==numel(bg.covariance_noise_seeds) && ...
    numel(bg.covariance_payload_seeds)>=3,'drdm:SampleCount','Need at least three paired moment draws.');
% 构造一次可复用的确认资源与独立条件业务集合。
known = drdm.lte.build_known_resources(cellInfo,cfg.known_resources);
model = drdm.model.build_conditional_model(cellInfo,known,bg.model,bg.mean_payload_seeds,bg.noise_power,bg.max_ensemble_bytes);
count = numel(bg.covariance_payload_seeds);
assert(16*(count+uncertaintyCount)*(numel(known.iq)+known.cpi_num_samples)<=bg.max_ensemble_bytes, ...
    'drdm:EnsembleMemory','Cached moment draws exceed the configured memory budget.');
moments = drdm.statistics.cache_conditional_draws(cellInfo,known,bg.model,bg.noise_power, ...
    bg.covariance_payload_seeds,bg.covariance_noise_seeds);
uncertaintyDraws = [];
if uncertaintyCount>0
    uncertaintyDraws = drdm.statistics.cache_conditional_draws(cellInfo,known,bg.model,bg.noise_power, ...
        u.payload_seeds,u.noise_seeds);
end
% 初始发现采用 RDMp；完整顺序拟合采用 options.search 指定的条件双图模型。
processing = struct('rdmp',cfg.rdmp,'rdmy',cfg.rdmy,'fractional_delay',cfg.fractional_delay, ...
    'known_resources',{cfg.known_resources},'estimation',drdm.estimation.check_options(cfg.estimation));
processing.estimation.fit_mode = 'rdmp_only';
processing.estimation.enabled = true;
processing.estimation.observation_covariance = [];
description = struct('code',drdm.pipeline.code_signature(),'cell_info',cellInfo,'cfg',processing,'options',options);
prepared = struct('kind','sequential_prepared','known',known,'model',model,'cfg',processing, ...
    'cell_info',cellInfo,'options',options,'moment_sources',moments.sources,'moment_noise',moments.noise, ...
    'uncertainty_draws',uncertaintyDraws,'internal_seeds',internalSeeds(:).', ...
    'procedure_id',drdm.statistics.fingerprint(description),'description',description);
end
