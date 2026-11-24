function interval = binomial_interval(successes,total,confidence)
%BINOMIAL_INTERVAL 计算二项比例的 Clopper–Pearson 区间。
%
% 输入：successes：事件次数；total：样本数；confidence：置信水平。
% 输出：interval：点估计、双侧区间和两个单侧界。
%
% 零事件和全事件分别按边界公式计算。

validateattributes(total,{'numeric'},{'scalar','integer','positive'});
validateattributes(successes,{'numeric'},{'scalar','integer','>=',0,'<=',total});
validateattributes(confidence,{'numeric'},{'scalar','real','>',0,'<',1});
alpha = 1-confidence;
low = 0;
high = 1;
oneLow = 0;
oneHigh = 1;
if successes>0
    low = betaincinv(alpha/2,successes,total-successes+1);
    oneLow = betaincinv(alpha,successes,total-successes+1);
end
if successes<total
    high = betaincinv(1-alpha/2,successes+1,total-successes);
    oneHigh = betaincinv(1-alpha,successes+1,total-successes);
end
interval = struct('estimate',successes/total,'two_sided',[low high], ...
    'one_sided_lower',oneLow,'one_sided_upper',oneHigh,'confidence',confidence, ...
    'successes',successes,'total',total);
end
