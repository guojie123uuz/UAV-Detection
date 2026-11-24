function assessment = assess_reference(fit,uncertainty,budget,existence,limits)
%ASSESS_REFERENCE 汇总局部路径分析的存在性、精度和误差预算。
%
% 输入：fit、uncertainty、budget、existence：四类分析结果；limits：设计要求。
% 输出：assessment：各项判据、原因及联合局部判定。
%
% 局部分析结果供诊断；完整顺序流程的参照标记由 pipeline.run 生成。

for name = {'max_delay_halfwidth_s','max_doppler_halfwidth_hz','max_prior_mse_fraction','min_signal_to_total_error'}
    validateattributes(limits.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
assert(strcmp(uncertainty.procedure_id,budget.procedure_id), ...
    'drdm:CalibrationMismatch','Uncertainty and error budget procedures differ.');
assert(isequaln(uncertainty.reference_paths,fit.paths) && isequaln(budget.reference_paths,fit.paths), ...
    'drdm:CalibrationMismatch','Uncertainty and budget must describe this fitted path set.');
count = numel(fit.paths);
reasons = {};
existencePass = false;
existenceValidated = false;
if ~isempty(existence)
    assert(strcmp(existence.context_id,uncertainty.procedure_id), ...
        'drdm:CalibrationMismatch','Existence evidence and refinement contexts differ.');
    assert(isequaln(existence.tested_paths,fit.paths), ...
        'drdm:CalibrationMismatch','Existence evidence tests a different fitted path set.');
    existenceValidated = existence.conditionally_validated;
    existencePass = existenceValidated && existence.accepted;
end
precisionPass = all(uncertainty.valid) && uncertainty.tail_resolved && ...
    all(uncertainty.halfwidth(1:2:2*count)<=limits.max_delay_halfwidth_s) && ...
    all(uncertainty.halfwidth(2:2:2*count)<=limits.max_doppler_halfwidth_hz);
budgetPass = budget.valid && all(budget.prior_mse_fraction<=limits.max_prior_mse_fraction) && ...
    all(budget.signal_to_total_error>=limits.min_signal_to_total_error);
numericalPass = fit.converged && fit.identifiable && ~any(fit.at_boundary);
if ~numericalPass
    reasons{end+1} = 'numerical_or_boundary_failure';
end
if ~existenceValidated
    reasons{end+1} = 'existence_uncalibrated';
elseif ~existencePass
    reasons{end+1} = 'below_existence_threshold';
end
if ~precisionPass
    reasons{end+1} = 'prior_uncertain';
end
if ~budgetPass
    reasons{end+1} = 'cancellation_budget_exceeded';
end
if ~uncertainty.coverage_validated
    reasons{end+1} = 'interval_coverage_unvalidated';
end
% 此处汇总固定路径数和单轮检验的局部诊断；整体参照状态由顺序流程给出。
reasons{end+1} = 'full_pipeline_unvalidated';
assessment = struct('reference_reliable',false,'existence_pass',existencePass, ...
    'precision_pass',precisionPass,'error_budget_pass',budgetPass,'numerical_pass',numericalPass, ...
    'conditional_gates_pass',existencePass && precisionPass && budgetPass && numericalPass, ...
    'reasons',{reasons},'limits',limits,'status','development_assessment', ...
    'scope','Quantitative gates reported; no promotion before selection coverage and complete search validation.');
end
