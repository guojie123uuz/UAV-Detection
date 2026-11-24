function map = compute_rdmp_iq(block,known,options)
%COMPUTE_RDMP_IQ 计算已知时域模板与接收信号的互模糊图。
%
% 输入：block：接收块；known：已知模板；options：图坐标、窗及算法。
% 输出：map：RDMp、归一化能量及边界信息。
%
% 该接口采用整数采样时延坐标，参考模板保留 CPI 前后保护区。

block = drdm.frontend.validate_block(block);
assert(block.fs_hz==known.fs_hz && numel(block.iq)==known.cpi_num_samples, ...
    'drdm:TemplateAlignment','Template sampling rate/CPI length mismatch.');
[delay,doppler] = drdm.rdm.check_axes(options.delay_s,options.doppler_hz);
lags = delay*block.fs_hz;
assert(all(abs(lags-round(lags))<1e-8),'drdm:IQGrid','IQ reference requires integer sample delays.');
w = drdm.rdm.window_vector(options.sample_window,numel(block.iq)).*block.valid_mask;
values = drdm.rdm.linear_correlation(block.iq.*w,known.iq,block.fs_hz, ...
    round(lags),doppler,known.cpi_first_sample-1,options.iq_method);
indices = known.cpi_first_sample+(0:numel(block.iq)-1);
energy = sum(w.*abs(known.iq(indices)).^2)/block.fs_hz;
metadata = struct('sample_window',options.sample_window,'method',options.iq_method, ...
    'valid_sample_count',sum(block.valid_mask),'reference_guard_before',known.cpi_first_sample-1, ...
    'boundary','zero-extended linear template');
map = drdm.rdm.pack_map(values,delay,doppler,'rdmp_iq',energy,metadata);
end
