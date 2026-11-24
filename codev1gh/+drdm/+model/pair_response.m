function response = pair_response(source,block,pathI,pathJ,cfg)
%PAIR_RESPONSE 计算有序路径对的单位增益互模糊响应。
%
% 输入：source：已知源或仿真源；block：观测坐标；pathI、pathJ：路径；cfg：算子配置。
% 输出：response：有序路径对响应及差分坐标。
%
% 组装 RDMy 时乘以 a_i·conj(a_j)，保留全部有序交叉项。

assert(isfield(source,'kind') && ismember(source.kind,{'known_only','oracle_full'}), ...
    'drdm:KernelSource','Source must be labelled known_only or oracle_full.');
assert(source.fs_hz==block.fs_hz && source.cpi_num_samples==numel(block.iq), ...
    'drdm:TemplateAlignment','Pair source does not match the observation block.');
pathI.gain_complex = 1;
pathJ.gain_complex = 1;
a = block;
b = block;
a.iq = drdm.simulation.apply_path(source,pathI,cfg.fractional_delay);
b.iq = drdm.simulation.apply_path(source,pathJ,cfg.fractional_delay);
response = drdm.rdm.cross_ambiguity(a,b,cfg.rdmy);
response.metadata.source_kind = source.kind;
response.metadata.pair_centre = [pathI.tau_s-pathJ.tau_s,pathI.nu_hz-pathJ.nu_hz];
response.metadata.coefficient = 'a_i * conj(a_j)';
end
