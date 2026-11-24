function bound = tolerance_limit(values,tailProbability,confidence)
%TOLERANCE_LIMIT 计算有限样本的单侧统计容忍限。
%
% 输入：values：非负统计量；tailProbability：目标上尾概率；confidence：置信水平。
% 输出：bound：顺序统计量上界、阶次、样本数和状态。
%
% 失败样本可记为 Inf；阶次由正则化不完全 beta 函数确定。

validateattributes(values,{'numeric'},{'vector','real','nonempty','nonnegative'});
assert(~any(isnan(values)),'drdm:CalibrationStatistic','Failures must be explicit +Inf, not silently omitted NaNs.');
validateattributes(tailProbability,{'numeric'},{'scalar','real','>',0,'<',1});
validateattributes(confidence,{'numeric'},{'scalar','real','>',0,'<',1});
n = numel(values);
order = nan;
limit = inf;
status = 'insufficient_samples';
for k = 1:n
    if betainc(1-tailProbability,k,n+1-k)<=1-confidence
        sorted = sort(values);
        limit = sorted(k);
        order = k;
        status = 'finite_tolerance_limit';
        if ~isfinite(limit)
            status = 'failure_or_infinite_tail';
        end
        break;
    end
end
bound = struct('limit',limit,'order',order,'sample_count',n,'tail_probability',tailProbability, ...
    'confidence',confidence,'status',status,'infinite_count',sum(isinf(values)));
end
