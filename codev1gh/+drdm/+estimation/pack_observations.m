function vector = pack_observations(context,p,y)
%PACK_OBSERVATIONS 将双图复数观测转换为固定实向量。
%
% 输入：context：索引与拟合模式；p、y：已归一化的 RDMp/RDMy。
% 输出：vector：实部/虚部堆叠后的列向量。
%
% RDMy 原点虚部恒为零，故在拟合向量中省去。

vector = [real(p(:));imag(p(:))];
if context.joint
    realPart = real(y(context.y_indices));
    imaginaryPart = imag(y(context.y_imag_indices));
    vector = [vector;realPart(:);imaginaryPart(:)];
end
end
