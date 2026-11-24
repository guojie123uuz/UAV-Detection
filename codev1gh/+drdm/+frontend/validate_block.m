function block = validate_block(block)
%VALIDATE_BLOCK 检查采样数据契约并计算 CPI 中心时刻。
%
% 输入：block：含 IQ、采样率、时间起点、有效掩码和参考系的数据块。
% 输出：block：双精度 IQ、置零的无效样本及 time_center_s。
%
% sample_index0 为零起始全局索引，时间字段以秒计。

required = {'iq','fs_hz','valid_mask','sample_index0','time_start_s', ...
    'reference_id','amplitude_units','corrections'};
assert(isstruct(block) && isscalar(block) && all(isfield(block,required)), ...
    'drdm:Block','Incomplete observation block.');
validateattributes(block.iq,{'numeric'},{'column','nonempty'});
validateattributes(block.fs_hz,{'numeric'},{'scalar','finite','real','positive'});
validateattributes(block.sample_index0,{'numeric'},{'scalar','finite','integer','nonnegative'});
validateattributes(block.time_start_s,{'numeric'},{'scalar','real','finite'});
assert(islogical(block.valid_mask) && isequal(size(block.valid_mask),size(block.iq)), ...
    'drdm:Mask','valid_mask must be a logical column matching IQ.');
assert(all(isfinite(block.iq(block.valid_mask))),'drdm:InvalidIQ','Valid samples must be finite.');
assert(any(block.valid_mask),'drdm:EmptyObservation','No valid samples.');
block.iq = double(block.iq);
block.iq(~block.valid_mask) = 0;
block.time_center_s = block.time_start_s+(numel(block.iq)-1)/(2*block.fs_hz);
end
