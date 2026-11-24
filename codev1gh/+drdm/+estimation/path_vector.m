function [vector,lower,upper] = path_vector(context,paths)
%PATH_VECTOR 将路径结构转换为缩放参数向量。
%
% 输入：context：时延/多普勒尺度和边界；paths：K 条路径。
% 输出：vector：4K×1 参数；lower、upper：参数边界。
%
% 顺序为各径时延/多普勒交替排列，再接全部增益实部及虚部。

coordinates = [[paths.tau_s]/context.delay_scale_s;[paths.nu_hz]/context.doppler_scale_hz];
gains = reshape([paths.gain_complex],[],1);
vector = [coordinates(:);real(gains);imag(gains)];
lower = [repmat(context.lower,numel(paths),1);-inf(2*numel(paths),1)];
upper = [repmat(context.upper,numel(paths),1);inf(2*numel(paths),1)];
end
