function dataset = run_suite(prepared,suite,role,seedBase,perProfile,threshold,progress)
%RUN_SUITE 在声明路径总体上重复完整参照筛选与评价。
%
% 输入：prepared、suite：准备对象和随机总体；role：样本用途；seedBase、perProfile：种子起点及每类数量；threshold：门限；progress：可选进度函数。
% 输出：dataset：逐例真值、路径评价、区间、预算及配置标识。
%
% role 可为 development、interval_training、threshold_training 或 evaluation；真值用于仿真与事后评价。

if nargin<7
    progress = [];
end
assert(ismember(role,{'development','interval_training','threshold_training','evaluation'}), ...
    'drdm:DatasetRole','Use an explicit training or independent evaluation role.');
assert(strcmp(prepared.description.code,drdm.pipeline.code_signature()), ...
    'drdm:CalibrationMismatch','Package code differs from the prepared procedure.');
validateattributes(perProfile,{'numeric'},{'scalar','integer','positive'});
validateattributes(seedBase,{'numeric'},{'scalar','integer','nonnegative'});
assert(numel(unique({suite.profiles.name}))==numel(suite.profiles), ...
    'drdm:Scenario','Profile names must be unique.');
total = numel(suite.profiles)*perProfile;
seeds = reshape(seedBase+(1:3*total),3,total);
drdm.statistics.check_seed_sets(prepared.internal_seeds,seeds(:));
template = struct('profile_index',0,'profile_name','','seeds',zeros(1,3),'truth',struct([]), ...
    'score',[],'background_id','','pilot_path_count',0);
records = repmat(template,1,total);
procedureId = '';
for k = 1:total
    index = ceil(k/perProfile);
    profile = suite.profiles(index);
    [block,truth] = drdm.evaluation.simulate_case(prepared,profile,seeds(:,k));
    trace = drdm.pipeline.trace(block,prepared);
    score = drdm.evaluation.score_trace(trace,truth.paths,threshold,suite.assignment);
    if isempty(procedureId)
        procedureId = trace.procedure_id;
    end
    assert(strcmp(procedureId,trace.procedure_id),'drdm:StatisticalDomain','Suite masks/observations differ.');
    records(k) = struct('profile_index',index,'profile_name',profile.name,'seeds',seeds(:,k).', ...
        'truth',truth,'score',score,'background_id',drdm.statistics.fingerprint(trace.context.background.covariance), ...
        'pilot_path_count',numel(trace.pilot.fit.paths));
    if ~isempty(progress)
        progress(k,total);
    end
end
dataset = struct('role',role,'prepared_id',prepared.procedure_id,'procedure_id',procedureId, ...
    'description',prepared.description,'suite',suite,'suite_id',drdm.statistics.fingerprint(suite), ...
    'records',records,'all_seeds',seeds(:).','per_profile',perProfile,'threshold',threshold, ...
    'scope','Every record reruns the entire data-dependent pipeline; truth used only in simulation/scoring.');
dataset.dataset_id = drdm.statistics.fingerprint(struct('role',role,'seeds',seeds, ...
    'procedure',procedureId,'suite',dataset.suite_id));
end
