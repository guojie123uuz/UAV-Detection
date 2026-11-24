function check_seed_sets(varargin)
%CHECK_SEED_SETS 检查各统计用途的随机种子是否互异。
%
% 输入：任意多个种子向量。
% 输出：无；重复或越界时抛出带标识的错误。
%
% 种子范围为 [0,2³²)，同一调用涵盖均值、背景、训练或评价用途。

allSeeds = zeros(sum(cellfun(@numel,varargin)),1);
first = 1;
for k = 1:nargin
    seeds = varargin{k};
    if isempty(seeds)
        continue;
    end
    validateattributes(seeds,{'numeric'},{'vector','real','finite','integer','>=',0,'<',2^32});
    allSeeds(first:first+numel(seeds)-1) = seeds(:);
    first = first+numel(seeds);
end
assert(numel(unique(allSeeds))==numel(allSeeds), ...
    'drdm:SampleReuse','Mean, covariance, bootstrap and test streams must have disjoint seeds.');
end
