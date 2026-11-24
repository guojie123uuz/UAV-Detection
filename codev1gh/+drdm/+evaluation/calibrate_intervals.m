function calibration = calibrate_intervals(dataset,designCoverage)
%CALIBRATE_INTERVALS 由独立区间训练样本确定区间放大系数。
%
% 输入：dataset：interval_training 数据集；designCoverage：可选设计覆盖率。
% 输出：calibration：multiplier、设计覆盖率、各类容忍限及训练标识。
%
% 设计覆盖率可高于 suite 的最终共同覆盖目标，调整后的选择流程由独立样本评价。

assert(strcmp(dataset.role,'interval_training'),'drdm:DatasetRole','Expected interval training data.');
% 可在训练前指定高于最终验收目标的设计覆盖率，为独立评价留出统计余量。
if nargin<2
    designCoverage = dataset.suite.minimum_joint_coverage;
end
validateattributes(designCoverage,{'numeric'}, ...
    {'scalar','real','finite','>=',dataset.suite.minimum_joint_coverage,'<',1});
positive = find(arrayfun(@(p)~isempty(p.path_boxes),dataset.suite.profiles));
assert(~isempty(positive),'drdm:CoveragePopulation','Include a nonempty path population.');
confidence = 1-(1-dataset.suite.family_confidence)/numel(positive);
limits = cell(1,numel(positive));
multiplier = dataset.description.options.interval_multiplier;
for k = 1:numel(positive)
    records = dataset.records([dataset.records.profile_index]==positive(k));
    scores = arrayfun(@(r)r.score.trace_interval_score,records);
    limits{k} = drdm.statistics.tolerance_limit(scores,1-designCoverage,confidence);
    multiplier = max(multiplier,limits{k}.limit);
end
calibration = struct('multiplier',multiplier,'usable',isfinite(multiplier),'limits',{limits}, ...
    'design_coverage',designCoverage, ...
    'base_description',dataset.description,'suite_id',dataset.suite_id,'all_seeds',dataset.all_seeds, ...
    'dataset_id',dataset.dataset_id,'scope','Training adjustment only; changed selection is evaluated independently.');
calibration.calibration_id = drdm.statistics.fingerprint(calibration);
end
