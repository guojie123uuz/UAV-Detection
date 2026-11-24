function map = cross_ambiguity(blockA,blockB,options)
%CROSS_AMBIGUITY 计算以乘积中点为时间坐标的对称互模糊函数。
%
% 输入：blockA、blockB：同一时间轴的两个信号块；options：整数延迟、Hz 频率与窗。
% 输出：map：多普勒×时延的复数互模糊图。
%
% A_AB(l,f)=Σ xA[n]conj(xB[n-l])exp(-j2πf(t[n]-l/2fs))/fs。

a = drdm.frontend.validate_block(blockA);
b = drdm.frontend.validate_block(blockB);
assert(numel(a.iq)==numel(b.iq) && a.fs_hz==b.fs_hz && ...
    a.time_start_s==b.time_start_s && a.sample_index0==b.sample_index0 && ...
    strcmp(a.reference_id,b.reference_id),'drdm:PairAlignment','Pair blocks have different time/reference axes.');
[lags,doppler] = drdm.rdm.check_axes(options.delay_samples,options.doppler_hz);
assert(all(lags==round(lags)),'drdm:SampleLag','RDMy needs integer sample lags.');
w = drdm.rdm.window_vector(options.sample_window,numel(a.iq));
xA = a.iq.*w.*a.valid_mask;
xB = b.iq.*w.*b.valid_mask;
values = drdm.rdm.linear_correlation(xA,xB,a.fs_hz,lags,doppler,0,options.method);
values = values.*exp(1i*pi*doppler*(lags/a.fs_hz));
energy = sqrt(sum(abs(xA).^2)*sum(abs(xB).^2))/a.fs_hz;
metadata = struct('sample_window',options.sample_window,'method',options.method, ...
    'boundary','zero-extended linear', ...
    'valid_sample_counts',[sum(a.valid_mask) sum(b.valid_mask)]);
map = drdm.rdm.pack_map(values,lags/a.fs_hz,doppler,'symmetric_cross_ambiguity',energy,metadata);
end
