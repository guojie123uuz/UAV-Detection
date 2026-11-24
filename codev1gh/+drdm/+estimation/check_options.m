function options = check_options(options)
%CHECK_OPTIONS 合并估计选项并检查模型及参数范围。
%
% 输入：options：完整或部分估计选项。
% 输出：options：补齐默认字段后的估计配置。
%
% 统一在估计入口处理配置，内部函数使用同一份选项。

defaults = drdm.estimation.default_options();
assert(isstruct(options) && isscalar(options),'drdm:EstimationConfig','Expected scalar options.');
fields = fieldnames(options);
assert(all(ismember(fields,fieldnames(defaults))),'drdm:EstimationConfig','Unknown estimation option.');
for k = 1:numel(fields)
    defaults.(fields{k}) = options.(fields{k});
end
options = defaults;
assert(islogical(options.enabled) && isscalar(options.enabled),'drdm:EstimationConfig','enabled must be logical.');
assert(ismember(options.fit_mode,{'rdmp_only','joint_known_only','joint_conditional'}), ...
    'drdm:FitModel','Expected rdmp_only, joint_known_only or joint_conditional.');
assert(ismember(options.finite_difference_type,{'forward','central'}), ...
    'drdm:FiniteDifference','Use forward or central finite differences.');
integerFields = {'max_paths','candidates_per_branch','refine_candidates','max_iterations','max_function_evaluations'};
for k = 1:numel(integerFields)
    validateattributes(options.(integerFields{k}),{'numeric'},{'scalar','integer','positive','finite'});
end
positiveFields = {'function_tolerance','step_tolerance','optimality_tolerance', ...
    'finite_difference_step','max_dictionary_condition','max_jacobian_condition', ...
    'duplicate_radius','candidate_merge_radius','pair_exclusion_radius'};
for k = 1:numel(positiveFields)
    validateattributes(options.(positiveFields{k}),{'numeric'},{'scalar','real','positive','finite'});
end
for field = {'min_total_improvement','min_relative_improvement','residual_fraction_stop'}
    validateattributes(options.(field{1}),{'numeric'},{'scalar','real','finite','>=',0,'<',1});
end
end
