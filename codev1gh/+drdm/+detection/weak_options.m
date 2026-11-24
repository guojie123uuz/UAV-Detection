function options = weak_options(overrides)
%WEAK_OPTIONS 建立一轮新增弱径搜索的配置。
%
% 输入：overrides：可选顶层字段覆盖结构。
% 输出：options：候选网格、精化数量、优化器及误差预算要求。
%
% 空网格采用 RDMp 坐标；optimizer 与 requirements 若覆盖，应提供该子结构的全部字段。

options = struct('candidates_per_branch',2,'refine_candidates',3, ...
    'grid_delay_s',[],'grid_doppler_hz',[],'max_grid_points',400, ...
    'candidate_merge_radius',0.15,'reference_match_radius',0.5, ...
    'optimizer',struct('max_iterations',120,'max_function_evaluations',2400, ...
    'finite_difference_type','central'), ...
    'require_budget',true,'requirements',struct('min_signal_to_error',10,'max_prior_fraction',0.1));
if nargin>=1
    fields = fieldnames(overrides);
    assert(all(ismember(fields,fieldnames(options))),'drdm:WeakOptions','Unknown weak-search option.');
    for k = 1:numel(fields)
        options.(fields{k}) = overrides.(fields{k});
    end
end
for field = {'candidates_per_branch','refine_candidates','max_grid_points'}
    validateattributes(options.(field{1}),{'numeric'},{'scalar','integer','finite','positive'});
end
for field = {'candidate_merge_radius','reference_match_radius'}
    validateattributes(options.(field{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for field = {'grid_delay_s','grid_doppler_hz'}
    value = options.(field{1});
    if ~isempty(value)
        validateattributes(value,{'numeric'},{'vector','real','finite'});
    end
end
assert(islogical(options.require_budget) && isscalar(options.require_budget), ...
    'drdm:WeakOptions','require_budget must be logical.');
for field = {'max_iterations','max_function_evaluations'}
    validateattributes(options.optimizer.(field{1}),{'numeric'},{'scalar','integer','finite','positive'});
end
assert(ismember(options.optimizer.finite_difference_type,{'central','forward'}), ...
    'drdm:WeakOptions','Unknown finite-difference type.');
for field = {'min_signal_to_error','max_prior_fraction'}
    validateattributes(options.requirements.(field{1}),{'numeric'},{'scalar','real','finite','positive'});
end
end
