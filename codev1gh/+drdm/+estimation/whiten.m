function result = whiten(values,cholesky)
%WHITEN 通过 Cholesky 因子对白化观测或字典。
%
% 输入：values：观测列向量或矩阵；cholesky：协方差下三角因子，可为空。
% 输出：result：L\values；因子为空时等于输入。
%
% 使用线性方程求解，避免显式计算协方差逆矩阵。

if isempty(cholesky)
    result = values;
else
    result = cholesky\values;
end
end
