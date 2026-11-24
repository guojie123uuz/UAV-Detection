function paths = paths_from_vector(context,vector)
%PATHS_FROM_VECTOR 将缩放参数向量还原为路径结构。
%
% 输入：context：参数尺度；vector：path_vector 定义的 4K 实向量。
% 输出：paths：含 tau_s、nu_hz、gain_complex 的结构数组。
%
% 时延和多普勒分别恢复为秒和 Hz。

assert(mod(numel(vector),4)==0,'drdm:PathParameters','Expected four real parameters per path.');
count = numel(vector)/4;
paths = repmat(struct('tau_s',0,'nu_hz',0,'gain_complex',0),1,count);
for k = 1:count
    paths(k) = struct('tau_s',vector(2*k-1)*context.delay_scale_s, ...
        'nu_hz',vector(2*k)*context.doppler_scale_hz, ...
        'gain_complex',vector(2*count+k)+1i*vector(3*count+k));
end
end
