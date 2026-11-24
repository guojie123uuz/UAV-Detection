function geometry = profile_weak_response(linearResponse,selfResponse,nuisanceJacobian,cholesky,relativeTolerance)
%PROFILE_WEAK_RESPONSE 计算强径局部重拟合后的有效弱径响应。
%
% 输入：linearResponse：实观测域 H；selfResponse：q；nuisanceJacobian：未白化强径雅可比；cholesky：协方差因子或空；relativeTolerance：可选秩容差。
% 输出：geometry：投影前后响应、增益估计算子、信息矩阵、秩及响应保留率。
%
% 采用 SVD 构造强径参数子空间，投影后的有效响应为 (I−JJ⁺)L⁻¹H。

if nargin<5
    relativeTolerance = 1e-8;
end
validateattributes(relativeTolerance,{'numeric'},{'scalar','real','finite','>',0,'<',1});
validateattributes(linearResponse,{'numeric'},{'2d','real','finite','ncols',2});
rows = size(linearResponse,1);
assert(size(nuisanceJacobian,1)==rows && numel(selfResponse)==rows, ...
    'drdm:ObservationDimension','Weak response and nuisance Jacobian must share observation rows.');
b = drdm.estimation.whiten(linearResponse,cholesky);
q = drdm.estimation.whiten(selfResponse(:),cholesky);
j = drdm.estimation.whiten(nuisanceJacobian,cholesky);
% 强径参数雅可比的列空间给出可被局部重拟合吸收的观测分量。
[u,s,~] = svd(j,0);
singularJ = diag(s);
rankJ = sum(singularJ>relativeTolerance*max([singularJ;0]));
u = u(:,1:rankJ);
% 将弱径响应投影到该空间的正交补，再判断两个增益分量是否可估。
c = b-u*(u.'*b);
[ub,sb,~] = svd(b,0);
singularB = diag(sb);
cutoff = relativeTolerance*max([singularB;0]);
rankB = sum(singularB>cutoff);
ub = ub(:,1:rankB);
[uc,sc,vc] = svd(c,0);
singularC = diag(sc);
% 以未投影响应的尺度定秩；否则消去后的舍入残差会被误认为独立弱径。
rankC = sum(singularC>cutoff);
uc = uc(:,1:rankC);
gainOperator = [];
retention = nan(2,1);
condition = inf;
if rankB==2
    % 在弱径响应的正交基中计算广义能量比，避免 BᵀB 放大条件数。
    % 对任意增益方向，投影前后的能量比上下界等于下式的特征值。
    effectiveBasis = ub-u*(u.'*ub);
    matrix = effectiveBasis.'*effectiveBasis;
    retention = sort(max(0,min(1,eig((matrix+matrix.')/2))));
end
if rankC==2
    gainOperator = (vc(:,1:2)./singularC(1:2).')*uc.';
    condition = singularC(1)/singularC(2);
end
geometry = struct('whitened_linear',b,'whitened_self',q,'effective_linear',c, ...
    'nuisance_basis',u,'fixed_basis',ub,'profiled_basis',uc,'gain_operator',gainOperator, ...
    'nuisance_rank',rankJ,'nuisance_parameter_count',size(j,2), ...
    'fixed_rank',rankB,'profiled_rank',rankC,'gain_estimable',rankC==2, ...
    'profiled_condition',condition,'information',c.'*c, ...
    'retained_fraction_bounds',retention,'relative_rank_tolerance',relativeTolerance, ...
    'cholesky',cholesky,'scope','Local nuisance projection with a fixed observation covariance.');
end
