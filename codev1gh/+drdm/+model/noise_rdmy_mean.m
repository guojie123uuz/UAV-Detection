function values = noise_rdmy_mean(block,options,noisePower)
%NOISE_RDMY_MEAN 计算复白噪声在有限观测窗中的 RDMy 均值。
%
% 输入：block：观测坐标与有效掩码；options：RDMy 坐标及窗；noisePower：E∣n∣²。
% 输出：values：与 RDMy 同尺寸的复数均值。
%
% 零时延切片含窗平方的傅里叶变换，非零时延均值为零。

validateattributes(noisePower,{'numeric'},{'scalar','real','finite','nonnegative'});
values = complex(zeros(numel(options.doppler_hz),numel(options.delay_samples)));
zero = find(options.delay_samples==0);
if isempty(zero) || noisePower==0
    return;
end
w = drdm.rdm.window_vector(options.sample_window,numel(block.iq)).*block.valid_mask;
t = ((0:numel(block.iq)-1)'-(numel(block.iq)-1)/2)/block.fs_hz;
for k = 1:numel(options.doppler_hz)
    values(k,zero) = noisePower/block.fs_hz*sum(w.^2.*exp(-2i*pi*options.doppler_hz(k)*t));
end
end
