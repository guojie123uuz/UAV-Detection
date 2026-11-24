function context = prepare_context(block,known,cfg,observations,conditionalModel)
%PREPARE_CONTEXT 建立固定观测域的参数估计上下文。
%
% 输入：block、known、cfg；observations：可选双图；conditionalModel：条件拟合所需模型。
% 输出：context：观测向量、尺度、坐标范围、模型和白化矩阵。
%
% 观测顺序为 RDMp 实/虚部，再接 RDMy 去重后的实/虚部。

block = drdm.frontend.validate_block(block);
assert(strcmp(known.kind,'known_only'),'drdm:EstimationSource','Only confirmed-resource templates are allowed.');
options = drdm.estimation.check_options(cfg.estimation);
if nargin<4 || isempty(observations)
    observations = struct('rdmp',drdm.rdm.compute_rdmp_re(block,known,cfg.rdmp), ...
        'rdmy',drdm.rdm.compute_rdmy(block,cfg.rdmy));
end
assert(strcmp(observations.rdmp.operator,'rdmp_ofdm_re') && ...
    strcmp(observations.rdmy.operator,'rdmy_symmetric'), ...
    'drdm:ObservationOperator','Expected the RE RDMp and symmetric RDMy operators.');
assert(isequal(observations.rdmp.delay_s,cfg.rdmp.delay_s(:).') && ...
    isequal(observations.rdmp.doppler_hz,cfg.rdmp.doppler_hz(:)) && ...
    isequal(observations.rdmy.delay_s,cfg.rdmy.delay_samples(:).'/block.fs_hz) && ...
    isequal(observations.rdmy.doppler_hz,cfg.rdmy.doppler_hz(:)), ...
    'drdm:ObservationAxes','Observation axes differ from the forward model.');
assert(numel(observations.rdmp.delay_s)>=2 && numel(observations.rdmp.doppler_hz)>=2, ...
    'drdm:SearchDomain','Continuous fitting needs a nonzero delay and Doppler search interval.');
assert(all(isfinite(observations.rdmp.complex_values(:))) && ...
    all(isfinite(observations.rdmy.complex_values(:))),'drdm:ObservationValues','Nonfinite observations.');
assert(strcmp(observations.rdmp.metadata.symbol_window,cfg.rdmp.symbol_window) && ...
    observations.rdmp.metadata.fft_backoff_samples==floor(min(known.cp_samples)*cfg.rdmp.fft_backoff_fraction) && ...
    strcmp(observations.rdmy.metadata.sample_window,cfg.rdmy.sample_window), ...
    'drdm:ObservationWindow','Observation windows differ from the forward model.');
% 上下文只保留观测和传播配置，模拟输入由独立数据结构管理。
modelCfg = struct('rdmp',cfg.rdmp,'rdmy',cfg.rdmy,'fractional_delay',cfg.fractional_delay);
context = struct('block',block,'known',known,'cfg',modelCfg,'options',options, ...
    'observations',observations,'joint',~strcmp(options.fit_mode,'rdmp_only'));
if strcmp(options.fit_mode,'joint_conditional')
    assert(nargin>=5 && isstruct(conditionalModel) && ...
        strcmp(conditionalModel.kind,'conditional_ensemble'), ...
        'drdm:ConditionalModel','joint_conditional requires an explicit independent source ensemble.');
    assert(conditionalModel.source.fs_hz==known.fs_hz && ...
        conditionalModel.source.cpi_first_sample==known.cpi_first_sample && ...
        conditionalModel.source.cpi_num_samples==known.cpi_num_samples && ...
        size(conditionalModel.iq_ensemble,1)==numel(known.iq), ...
        'drdm:ConditionalAlignment','Conditional source and known template must share their coordinates.');
    assert(strcmp(conditionalModel.known_signature, ...
        drdm.statistics.fingerprint(struct('grid',known.grid,'fs_hz',known.fs_hz))), ...
        'drdm:ConditionalResources','Conditional source was built for different confirmed resources.');
    context.conditional_model = conditionalModel;
    context.conditional_noise_mean = drdm.model.noise_rdmy_mean(block,cfg.rdmy,conditionalModel.noise_power);
    n = numel(block.iq);
    assert(16*(n*numel(cfg.rdmy.doppler_hz)+numel(conditionalModel.iq_ensemble))<=conditionalModel.max_working_bytes, ...
        'drdm:EnsembleMemory','Conditional Doppler mixer exceeds the configured memory budget.');
    t = ((0:n-1)'-(n-1)/2)/block.fs_hz;
    context.conditional_mixer = exp(-2i*pi*t*cfg.rdmy.doppler_hz(:).');
end
% 按有效带宽和 CPI 时长缩放参数，改善非线性优化的数值尺度。
active = any(known.observation_grid~=0,2);
frequency = known.subcarrier_hz(active);
assert(numel(frequency)>1,'drdm:ResourceBandwidth','Insufficient occupied bandwidth for delay estimation.');
context.delay_scale_s = 1/(max(frequency)-min(frequency));
context.doppler_scale_hz = block.fs_hz/numel(block.iq);
context.lower = [min(observations.rdmp.delay_s)/context.delay_scale_s; ...
    min(observations.rdmp.doppler_hz)/context.doppler_scale_hz];
context.upper = [max(observations.rdmp.delay_s)/context.delay_scale_s; ...
    max(observations.rdmp.doppler_hz)/context.doppler_scale_hz];
context.p_scale = observations.rdmp.normalization_energy;
assert(context.p_scale>0,'drdm:NoKnownObservation','No known-resource normalization energy.');
% 剔除共轭镜像重复项及 RDMy 原点的恒零虚部。
context.y_indices = find(observations.rdmy.metadata.independent_half_plane);
[f,d] = ind2sub(size(observations.rdmy.complex_values),context.y_indices);
origin = observations.rdmy.delay_s(d).'==0 & observations.rdmy.doppler_hz(f)==0;
context.y_imag_indices = context.y_indices(~origin);
indices = known.cpi_first_sample+(0:known.cpi_num_samples-1);
w = drdm.rdm.window_vector(cfg.rdmy.sample_window,numel(block.iq)).*block.valid_mask;
context.y_scale = sum(abs(known.iq(indices).*w).^2)/block.fs_hz;
assert(context.y_scale>0,'drdm:NoKnownObservation','No known-waveform energy.');
context.observed = drdm.estimation.pack_observations(context, ...
    observations.rdmp.complex_values/context.p_scale, ...
    observations.rdmy.complex_values/context.y_scale);
context.cholesky = [];
context.p_cholesky = [];
% 可选协方差在固定实观测坐标中白化；各候选使用同一权重。
covariance = options.observation_covariance;
if ~isempty(covariance)
    validateattributes(covariance,{'numeric'},{'2d','real','finite','square'});
    assert(size(covariance,1)==numel(context.observed), ...
        'drdm:CovarianceShape','Covariance must match the fixed real observation vector.');
    assert(norm(covariance-covariance','fro')<=1e-12*max(norm(covariance,'fro'),realmin), ...
        'drdm:CovarianceSymmetry','Covariance must be symmetric.');
    [context.cholesky,flag] = chol(covariance,'lower');
    assert(flag==0,'drdm:CovariancePositive','Covariance must be positive definite.');
    count = 2*numel(observations.rdmp.complex_values);
    context.p_cholesky = chol(covariance(1:count,1:count),'lower');
end
context.target = drdm.estimation.whiten(context.observed,context.cholesky);
context.initial_cost = sum(context.target.^2);
if isfield(context,'conditional_model')
    zeroModel = drdm.estimation.pack_observations(context, ...
        zeros(size(observations.rdmp.complex_values)),context.conditional_noise_mean/context.y_scale);
    context.initial_cost = sum(drdm.estimation.whiten(zeroModel-context.observed,context.cholesky).^2);
end
end
