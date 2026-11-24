function calibration = calibrate_thresholds(trials,falsePathProbability,confidence)
%CALIBRATE_THRESHOLDS 由单轮搜索样本估计存在性门限。
%
% 输入：trials：独立原假设样本；falsePathProbability：目标概率；confidence：置信水平。
% 输出：calibration：门限、顺序阶次、尾概率上界和配置标识。
%
% 采用统计量严格大于门限的接纳规则。

validateattributes(falsePathProbability,{'numeric'},{'scalar','real','>',0,'<',1});
validateattributes(confidence,{'numeric'},{'scalar','real','>',0,'<',1});
drdm.statistics.check_seed_sets(trials.payload_seeds,trials.noise_seeds);
n = numel(trials.statistics);
assert(n>0 && numel(trials.complete)==n && numel(trials.payload_seeds)==n && numel(trials.noise_seeds)==n, ...
    'drdm:CalibrationSamples','Inconsistent calibration records.');
threshold = inf;
order = nan;
tailUpper = 1;
status = 'insufficient_tail_samples';
if ~all(trials.complete) || ~all(isfinite(trials.statistics))
    status = 'incomplete_null_search';
else
    assert(all(trials.statistics>=0),'drdm:CalibrationStatistic','Expected nonnegative improvement maxima.');
    for k = 1:n
        if betainc(1-falsePathProbability,k,n+1-k)<=1-confidence
            sorted = sort(trials.statistics);
            threshold = sorted(k);
            order = k;
            tailUpper = 1-betaincinv(1-confidence,k,n+1-k);
            status = 'conditional_threshold';
            break;
        end
    end
end
calibration = struct('status',status,'threshold',threshold,'order_index',order,'sample_count',n, ...
    'target_false_path_probability',falsePathProbability,'confidence',confidence, ...
    'calibration_tail_upper_bound',tailUpper,'procedure_id',trials.procedure_id,'context_id',trials.context_id, ...
    'payload_seeds',trials.payload_seeds,'noise_seeds',trials.noise_seeds, ...
    'comparison','strictly greater than threshold','scope',trials.scope, ...
    'full_pipeline_calibrated',false);
calibration.calibration_id = drdm.statistics.fingerprint(calibration);
end
