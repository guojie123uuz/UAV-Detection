function map = compute_rdmp_re(block,known,options)
%COMPUTE_RDMP_RE 在已确认 OFDM 资源上计算 RDMp。
%
% 输入：block：接收块；known：已知 RE 与符号坐标；options：搜索网格及 FFT 窗参数。
% 输出：map：复数 RDMp、CP 支持范围、有效符号与归一化信息。
%
% 每个候选多普勒先在采样域补偿，再 FFT；时延通过子载波相位扫描。

block = drdm.frontend.validate_block(block);
assert(block.fs_hz==known.fs_hz && numel(block.iq)==known.cpi_num_samples, ...
    'drdm:TemplateAlignment','Template sampling rate/CPI length mismatch.');
[delay,doppler] = drdm.rdm.check_axes(options.delay_s,options.doppler_hz);
assert(all(abs(doppler)<=block.fs_hz/2),'drdm:DopplerAlias','Doppler exceeds sample Nyquist range.');
validateattributes(options.fft_backoff_fraction,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
% FFT 起点向循环前缀回退，以容纳正、负相对时延。
backoff = floor(min(known.cp_samples)*options.fft_backoff_fraction);
validDelay = [-backoff,min(known.cp_samples)-backoff]/block.fs_hz;
assert(min(delay)>=validDelay(1)-eps && max(delay)<=validDelay(2)+eps, ...
    'drdm:CPRange','RE delay grid is outside [%g,%g] seconds; use IQ reference or change the FFT window.', ...
validDelay(1),validDelay(2));
indices = (0:known.nfft-1)'+(known.useful_start_samples(:).'-backoff);
assert(all(indices(:)>=1 & indices(:)<=numel(block.iq)), ...
    'drdm:FFTBounds','FFT windows exceed the observation block.');
validSymbols = all(block.valid_mask(indices),1);
weights = drdm.rdm.window_vector(options.symbol_window,size(indices,2)).'.*validSymbols;
reference = known.observation_grid;
energy = sum(sum(abs(reference).^2.*weights))/(block.fs_hz*known.nfft);
assert(energy>0,'drdm:NoKnownObservation','No usable known-resource symbols.');
% 子载波相位扫描实现时延匹配，phase 补偿 FFT 窗回退引入的确定相位。
steering = exp(1i*2*pi*delay(:)*known.subcarrier_hz.');
phase = exp(1i*2*pi*known.subcarrier_hz*backoff/block.fs_hz);
t = ((0:numel(block.iq)-1)'-(numel(block.iq)-1)/2)/block.fs_hz;
values = complex(zeros(numel(doppler),numel(delay)));
for f = 1:numel(doppler)
    corrected = block.iq.*exp(-1i*2*pi*doppler(f)*t);
    spectrum = fftshift(fft(corrected(indices),[],1),1);
    received = spectrum(known.fft_indices,:).*phase;
    matched = sum(received.*conj(reference).*weights,2);
    values(f,:) = (steering*matched).'/ (block.fs_hz*known.nfft);
end
metadata = struct('symbol_window',options.symbol_window,'fft_backoff_samples',backoff, ...
    'cp_valid_delay_s',validDelay,'valid_symbols',validSymbols,'symbol_weights',weights, ...
    'symbol_centres_s',((known.useful_start_samples-backoff-1+(known.nfft-1)/2)- ...
    (numel(block.iq)-1)/2)/block.fs_hz, ...
    'doppler_compensation','per-candidate sample mixing before FFT', ...
    'other_paths_cp_support','not inferable from confirmed resources alone', ...
    'missing_data_policy','drop any FFT symbol containing invalid samples');
map = drdm.rdm.pack_map(values,delay,doppler,'rdmp_ofdm_re',energy,metadata);
end
