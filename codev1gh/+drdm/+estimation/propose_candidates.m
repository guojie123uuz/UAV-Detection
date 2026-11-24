function candidates = propose_candidates(pResidual,yMap,paths,context)
%PROPOSE_CANDIDATES 由 RDMp 峰和 RDMy 径对差分产生路径候选。
%
% 输入：pResidual、yMap：残差或候选图；paths：已有路径；context：搜索域与尺度。
% 输出：candidates：坐标、来源、参考路径、差分方向及峰值分数。
%
% 差分同时考虑两个方向，经范围裁剪和尺度化距离去重。

options = context.options;
candidates = repmat(make_candidate(0,0,'',0,0,0),1, ...
    options.candidates_per_branch*(1+2*numel(paths)));
candidateCount = 0;
pPower = abs(pResidual.complex_values);
[tau,nu] = meshgrid(pResidual.delay_s,pResidual.doppler_hz);
for k = 1:numel(paths)
    distance = hypot((tau-paths(k).tau_s)/context.delay_scale_s, ...
        (nu-paths(k).nu_hz)/context.doppler_scale_hz);
    pPower(distance<options.duplicate_radius) = 0;
end
peaks = local_peaks(pPower,options.candidates_per_branch);
for k = 1:size(peaks,1)
    candidateCount = candidateCount+1;
    candidates(candidateCount) = make_candidate(pResidual.delay_s(peaks(k,2)), ...
        pResidual.doppler_hz(peaks(k,1)),'rdmp',0,0,peaks(k,3));
end
% 已有路径提供 RDMy 差分的参考坐标，两种差分方向均生成候选。
if ~isempty(paths)
    yPower = abs(yMap.complex_values);
    yPower(~yMap.metadata.independent_half_plane) = 0;
    [delay,frequency] = meshgrid(yMap.delay_s,yMap.doppler_hz);
    % 候选搜索时屏蔽已有路径对的中心邻域；完整响应在预测步骤计算。
    for i = 1:numel(paths)
        for j = 1:numel(paths)
            distance = hypot((delay-(paths(i).tau_s-paths(j).tau_s))/context.delay_scale_s, ...
                (frequency-(paths(i).nu_hz-paths(j).nu_hz))/context.doppler_scale_hz);
            yPower(distance<options.pair_exclusion_radius) = 0;
        end
    end
    peaks = local_peaks(yPower,options.candidates_per_branch);
    for k = 1:size(peaks,1)
        for reference = 1:numel(paths)
            for direction = [-1 1]
                t = paths(reference).tau_s+direction*yMap.delay_s(peaks(k,2));
                f = paths(reference).nu_hz+direction*yMap.doppler_hz(peaks(k,1));
                candidateCount = candidateCount+1;
                candidates(candidateCount) = make_candidate(t,f,'rdmy_difference',reference,direction,peaks(k,3));
            end
        end
    end
end
% 合并候选并按搜索边界、已有路径距离和候选间距离去重。
candidates = candidates(1:candidateCount);
keep = false(size(candidates));
acceptedCoordinates = zeros(2,candidateCount);
acceptedCount = 0;
existing = zeros(2,numel(paths));
if ~isempty(paths)
    existing = [[paths.tau_s]/context.delay_scale_s;[paths.nu_hz]/context.doppler_scale_hz];
end
for k = 1:numel(candidates)
    c = [candidates(k).tau_s/context.delay_scale_s;candidates(k).nu_hz/context.doppler_scale_hz];
    if any(c<context.lower | c>context.upper)
        continue;
    end
    if any(vecnorm(existing-c,2,1)<options.duplicate_radius)
        continue;
    end
    if any(vecnorm(acceptedCoordinates(:,1:acceptedCount)-c,2,1)<options.candidate_merge_radius)
        continue;
    end
    keep(k) = true;
    acceptedCount = acceptedCount+1;
    acceptedCoordinates(:,acceptedCount) = c;
end
candidates = candidates(keep);
end

function candidate = make_candidate(tau,nu,source,reference,direction,score)
candidate = struct('tau_s',tau,'nu_hz',nu,'gain_complex',0,'source',source, ...
    'reference_index',reference,'pair_sign',direction,'seed_score',score);
end

function peaks = local_peaks(power,limit)
peaks = zeros(numel(power),3);
count = 0;
for row = 1:size(power,1)
    for column = 1:size(power,2)
        value = power(row,column);
        neighbors = power(max(1,row-1):min(end,row+1),max(1,column-1):min(end,column+1));
        if isfinite(value) && value>0 && value>=max(neighbors(:))
            count = count+1;
            peaks(count,:) = [row column value];
        end
    end
end
peaks = peaks(1:count,:);
[~,order] = sort(peaks(:,3),'descend');
peaks = peaks(order(1:min(limit,count)),:);
end
