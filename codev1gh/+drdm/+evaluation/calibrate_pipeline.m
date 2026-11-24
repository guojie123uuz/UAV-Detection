function calibration = calibrate_pipeline(dataset,intervalCalibration)
%CALIBRATE_PIPELINE 按完整参照选择流程标定存在性门限。
%
% 输入：dataset：threshold_training 数据；intervalCalibration：已冻结的区间训练结果。
% 输出：calibration：门限、各场景容忍限、配置、总体及训练标识。
%
% 顺序模式评价所有可接纳前缀，联合模式评价最终整体模型。

assert(strcmp(dataset.role,'threshold_training'),'drdm:DatasetRole','Expected threshold training data.');
assert(intervalCalibration.usable && strcmp(intervalCalibration.suite_id,dataset.suite_id), ...
    'drdm:CalibrationMismatch','Interval training is not usable for this population.');
drdm.statistics.check_seed_sets(intervalCalibration.all_seeds,dataset.all_seeds);
expected = intervalCalibration.base_description;
expected.options.interval_multiplier = intervalCalibration.multiplier;
assert(isequaln(expected,dataset.description),'drdm:CalibrationMismatch', ...
    'Only the trained interval multiplier may change before threshold training.');
profiles = dataset.suite.profiles;
confidence = 1-(1-dataset.suite.family_confidence)/numel(profiles);
limits = cell(1,numel(profiles));
threshold = dataset.description.options.minimum_threshold;
for k = 1:numel(profiles)
    records = dataset.records([dataset.records.profile_index]==k);
    scores = arrayfun(@(r)r.score.max_false_score,records);
    limits{k} = drdm.statistics.tolerance_limit(scores,dataset.suite.target_false_path_probability,confidence);
    threshold = max(threshold,limits{k}.limit);
end
calibration = struct('threshold',threshold,'usable',isfinite(threshold),'limits',{limits}, ...
    'prepared_id',dataset.prepared_id,'procedure_id',dataset.procedure_id,'description',dataset.description, ...
    'suite',dataset.suite,'suite_id',dataset.suite_id,'all_seeds',dataset.all_seeds, ...
    'interval_calibration',intervalCalibration,'dataset_id',dataset.dataset_id, ...
    'scope','Per-profile average error for the frozen reference-selection procedure; independent validation follows training.');
calibration.calibration_id = drdm.statistics.fingerprint(calibration);
end
