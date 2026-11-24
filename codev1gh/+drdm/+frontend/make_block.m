function block = make_block(iq,fs,metadata)
%MAKE_BLOCK 封装同步后的单通道复数采样数据。
%
% 输入：iq：采样向量；fs：采样率 Hz；metadata：可选时间与校正记录。
% 输出：block：统一的 IQ 数据块，iq 和 valid_mask 均为 N×1。
%
% 采样间隔恒定；无效样本由 valid_mask 标记。

if nargin<3
    metadata = struct();
end
assert(isvector(iq) && ~isempty(iq),'drdm:IQShape','Expected one nonempty IQ vector.');
count = numel(iq);
block = struct('iq',double(iq(:)),'fs_hz',fs,'valid_mask',true(count,1), ...
    'sample_index0',0,'time_start_s',0,'reference_id','synchronized_input', ...
    'amplitude_units','native IQ','corrections',struct(), ...
    'version','drdm-block-v1');
fields = fieldnames(metadata);
allowed = {'valid_mask','sample_index0','time_start_s','reference_id', ...
    'amplitude_units','corrections'};
assert(all(ismember(fields,allowed)),'drdm:BlockMetadata','Unsupported metadata field.');
for k = 1:numel(fields)
    block.(fields{k}) = metadata.(fields{k});
end
block = drdm.frontend.validate_block(block);
end
