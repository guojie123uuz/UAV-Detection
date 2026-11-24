function map = compute_rdmy(block,options)
%COMPUTE_RDMY 计算完整接收 IQ 的对称自模糊图 RDMy。
%
% 输入：block：完整接收数据块；options：时延、多普勒、窗和相关算法。
% 输出：map：RDMy 与共轭镜像去重掩码。
%
% 满足 A(-τ,-ν)=conj(A(τ,ν))；掩码仅去除共轭重复坐标。

map = drdm.rdm.cross_ambiguity(block,block,options);
map.operator = 'rdmy_symmetric';
canonical = map.delay_s>0 | (map.delay_s==0 & map.doppler_hz>=0);
mirrorExists = ismember(-map.delay_s,map.delay_s) & ismember(-map.doppler_hz,map.doppler_hz);
map.metadata.independent_half_plane = canonical | ~mirrorExists;
map.metadata.independence_scope = 'Conjugate duplicates only; other statistical correlations remain.';
map.metadata.mirror_relation = 'A(-delay,-doppler)=conj(A(delay,doppler))';
end
