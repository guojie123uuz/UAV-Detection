function association = match_reference_paths(reference,fitted,scales,radius)
%MATCH_REFERENCE_PATHS 将联合重估的路径与原参照一一对应。
%
% 输入：reference、fitted：原参照及多一条路径的拟合；scales：秒和 Hz 尺度；radius：尺度化对应半径。
% 输出：association：complete、reference_indices、weak_index 及距离代价。
%
% 在对应半径内穷举最小平方距离分配，支持最多六条拟合路径。

n = numel(reference);
m = numel(fitted);
assert(m==n+1 && m<=6,'drdm:PathAssociation','Expected at most six paths with one addition.');
distance = zeros(n,m);
for k = 1:n
    distance(k,:) = hypot(([fitted.tau_s]-reference(k).tau_s)/scales(1), ...
        ([fitted.nu_hz]-reference(k).nu_hz)/scales(2));
end
best = [];
bestCost = inf;
visit(1,zeros(1,n),false(1,m),0);
complete = ~isempty(best);
unmatched = [];
if complete
    unmatched = setdiff(1:m,best);
end
association = struct('complete',complete,'reference_indices',best,'weak_index',unmatched,'cost',bestCost);
    function visit(k,indices,used,cost)
        if k>n
            if cost<bestCost
                best = indices;
                bestCost = cost;
            end
            return;
        end
        for j = 1:m
            if ~used(j) && distance(k,j)<=radius
                next = indices;
                next(k) = j;
                mask = used;
                mask(j) = true;
                visit(k+1,next,mask,cost+distance(k,j)^2);
            end
        end
    end
end
