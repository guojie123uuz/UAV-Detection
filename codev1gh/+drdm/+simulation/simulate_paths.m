function [block,truth] = simulate_paths(source,paths,cfg)
%SIMULATE_PATHS 叠加共享源的多径分量和复高斯白噪声。
%
% 输入：source：源波形；paths：路径结构数组；cfg：传播及噪声配置。
% 输出：block：接收数据块；truth：路径、噪声与功率记录。
%
% 噪声功率定义为 E∣n∣²；独立 RandStream 保持全局随机状态。

assert(isstruct(paths) && all(isfield(paths,{'tau_s','nu_hz','gain_complex'})), ...
    'drdm:Paths','Paths need tau_s, nu_hz and gain_complex.');
validateattributes(cfg.noise_power,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(cfg.noise_seed,{'numeric'},{'scalar','integer','>=',0,'<',2^32});
y = complex(zeros(source.cpi_num_samples,1));
for k = 1:numel(paths)
    y = y+drdm.simulation.apply_path(source,paths(k),cfg.fractional_delay);
end
stream = RandStream('mt19937ar','Seed',cfg.noise_seed);
noise = sqrt(cfg.noise_power/2)*(randn(stream,size(y))+1i*randn(stream,size(y)));
block = drdm.frontend.make_block(y+noise,source.fs_hz, ...
    struct('reference_id','simulation_synchronized','amplitude_units',source.amplitude_units));
truth = struct('paths',paths,'noise_power',cfg.noise_power, ...
    'noise_seed',cfg.noise_seed,'noiseless_power',mean(abs(y).^2), ...
    'phase_reference','CPI centre');
end
