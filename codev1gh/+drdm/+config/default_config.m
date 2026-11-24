function cfg = default_config()
%DEFAULT_CONFIG 建立 LTE 波形、观测图与参数估计的默认配置。
%
% 输入：无。
% 输出：cfg：LTE、小区时序、双图坐标、传播与估计参数。
%
% 时延单位为秒，多普勒单位为 Hz，复增益相位参考 CPI 中心。

cfg.version = 'drdm-observations-v1';
cfg.lte = struct('NDLRB',6,'NCellID',42,'CellRefP',1, ...
    'CyclicPrefix','Normal','DuplexMode','FDD','NFrame',0,'NSubframe',0, ...
    'Ng','Sixth','PHICHDuration','Normal','Windowing',0);
cfg.num_subframes = 20;
cfg.guard_subframes = 1;
cfg.payload_seed = 20261113;
cfg.noise_seed = 20261114;
cfg.noise_power = 0; % 复噪声功率 E|n|²，采用 LTE IFFT 原始幅度单位
cfg.fractional_delay = struct('radius',64,'kaiser_beta',10);
cfg.known_resources = {'pss','sss','crs','pbch'};
cfg.rdmp = struct('delay_s',(-4:4)/1.92e6,'doppler_hz',(-250:25:250)', ...
    'iq_method','fft','sample_window','rectangular', ...
    'symbol_window','rectangular','fft_backoff_fraction',0.5);
cfg.rdmy = struct('delay_samples',-10:10,'doppler_hz',(-500:25:500)', ...
    'method','fft','sample_window','rectangular');
cfg.estimation = drdm.estimation.default_options();
end
