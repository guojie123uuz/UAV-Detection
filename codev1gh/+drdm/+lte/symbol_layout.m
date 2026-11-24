function layout = symbol_layout(enb,totalSubframes)
%SYMBOL_LAYOUT 计算 OFDM 符号边界与有效子载波索引。
%
% 输入：enb：小区配置；totalSubframes：子帧总数。
% 输出：layout：CP/有效符号起点、FFT 索引、子载波频率和采样率。
%
% 样本索引从 1 开始；子载波索引跳过直流频点。

info = drdm.lte.validate_cell(enb);
validateattributes(totalSubframes,{'numeric'},{'scalar','integer','positive','finite'});
cp = repmat(double(info.CyclicPrefixLengths(:)),totalSubframes,1);
nfft = double(info.Nfft);
starts = [1;1+cumsum(nfft+cp(1:end-1))];
layout = struct('fs_hz',double(info.SamplingRate),'nfft',nfft, ...
    'cp_samples',cp,'cp_start_samples',starts, ...
    'useful_start_samples',starts+cp,'num_samples',sum(cp+nfft));
k = [-enb.NDLRB*6:-1 1:enb.NDLRB*6]';
layout.subcarrier_hz = k*15000;
layout.fft_indices = nfft/2+1+k;
layout.samples_per_subframe = sum(double(info.CyclicPrefixLengths))+14*nfft;
end
