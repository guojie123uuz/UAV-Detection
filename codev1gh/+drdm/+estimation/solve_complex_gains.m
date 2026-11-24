function [gains,diagnostics] = solve_complex_gains(dictionary,data,cholesky)
%SOLVE_COMPLEX_GAINS 以加权奇异值分解求解全部路径的复增益。
%
% 输入：dictionary：复数响应字典；data：复观测向量；cholesky：可选实堆叠协方差因子。
% 输出：gains：K×1 复增益；diagnostics：秩、条件数及奇异值。
%
% 实/虚部联合表示保留两分量之间的误差相关性。

if nargin<3
    cholesky = [];
end
count = size(dictionary,2);
matrix = [real(dictionary) -imag(dictionary);imag(dictionary) real(dictionary)];
matrix = drdm.estimation.whiten(matrix,cholesky);
target = drdm.estimation.whiten([real(data(:));imag(data(:))],cholesky);
if count==0
    gains = complex(zeros(0,1));
    diagnostics = struct('rank',0,'condition',1,'singular_values',[]);
    return;
end
[u,s,v] = svd(matrix,0);
singular = diag(s);
threshold = max(size(matrix))*eps(max([singular;1]));
keep = singular>threshold;
parameters = v(:,keep)*((u(:,keep)'*target)./singular(keep));
gains = parameters(1:count)+1i*parameters(count+1:end);
condition = inf;
if sum(keep)==2*count
    condition = singular(1)/singular(end);
end
diagnostics = struct('rank',sum(keep),'condition',condition,'singular_values',singular);
end
