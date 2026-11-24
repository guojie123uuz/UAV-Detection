function assignment = assign_paths(estimated,truth,delayTolerance,dopplerTolerance)
%ASSIGN_PATHS 求解门控条件下的最优路径一一对应。
%
% 输入：estimated、truth：路径数组；delayTolerance、dopplerTolerance：归一化尺度。
% 输出：assignment：真值索引、匹配/多检/漏检数量及代价。
%
% 椭圆门控内先最大化匹配数，再最小化平方距离；最多各 6 条路径。

validateattributes(delayTolerance,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(dopplerTolerance,{'numeric'},{'scalar','real','finite','positive'});
n = numel(estimated);
m = numel(truth);
assert(n<=6 && m<=6,'drdm:AssignmentSize','Exact benchmark assignment supports up to six paths.');
cost = zeros(n,m);
for i = 1:n
    cost(i,:) = hypot(([truth.tau_s]-estimated(i).tau_s)/delayTolerance, ...
        ([truth.nu_hz]-estimated(i).nu_hz)/dopplerTolerance);
end
bestCount = -1;
bestCost = inf;
best = zeros(1,n);
% 枚举至多六条路径的一一对应，优先保留更多有效匹配。
visit(1,false(1,m),zeros(1,n),0,0);
assignment = struct('truth_index',best,'matched_count',bestCount,'false_count',sum(best==0), ...
    'missed_count',m-bestCount,'cost',bestCost);

    function visit(i,used,current,count,total)
        if i>n
            if count>bestCount || (count==bestCount && total<bestCost)
                bestCount = count;
                bestCost = total;
                best = current;
            end
            return;
        end
        visit(i+1,used,current,count,total);
        for j = 1:m
            if ~used(j) && cost(i,j)<=1
                chosen = current;
                chosen(i) = j;
                taken = used;
                taken(j) = true;
                visit(i+1,taken,chosen,count+1,total+cost(i,j)^2);
            end
        end
    end
end
