function [order,costs] = rank_candidates(context,paths,candidates)
%RANK_CANDIDATES 按固定观测域的加权残差排序候选。
%
% 输入：context、paths、candidates：上下文、已有路径与候选。
% 输出：order、costs：候选次序与对应残差代价。
%
% 固定坐标处联合求解复增益，再按残差升序排列候选。

costs = inf(1,numel(candidates));
pData = context.observations.rdmp.complex_values(:)/context.p_scale;
for k = 1:numel(candidates)
    candidate = candidates(k);
    seed = struct('tau_s',candidate.tau_s,'nu_hz',candidate.nu_hz,'gain_complex',0);
    seeds = [paths seed];
    coordinates = [[seeds.tau_s]/context.delay_scale_s; ...
        [seeds.nu_hz]/context.doppler_scale_hz];
    [dictionary,kernels] = drdm.estimation.unit_dictionary(context,coordinates);
    gains = drdm.estimation.solve_complex_gains(dictionary,pData,context.p_cholesky);
    model = drdm.estimation.model_vector(context,dictionary,kernels,gains);
    residual = drdm.estimation.whiten(model-context.observed,context.cholesky);
    costs(k) = sum(residual.^2);
end
[~,order] = sort(costs);
end
