function scan = scan_weak_path(context,referencePaths,options)
%SCAN_WEAK_PATH 执行一轮弱径候选生成、粗排及强弱径联合精化。
%
% 输入：context：含独立误差样本的联合上下文；referencePaths：一至五条参照；options：可选弱径配置。
% 输出：scan：原假设拟合、全部候选、精化记录、最佳弱径、统计量及完成状态。
%
% 原假设和备择假设均重估强径；粗排选定的每个候选均须完成优化。

if nargin<3
    options = struct();
end
options = drdm.detection.weak_options(options);
assert(context.joint && ~isempty(referencePaths) && numel(referencePaths)<=5, ...
    'drdm:WeakReferences','Use one to five references in a joint observation context.');
assert(isfield(context,'uncertainty_draws') && ~isempty(context.uncertainty_draws), ...
    'drdm:UncertaintySamples','Weak search requires independent conditional error draws.');
work = context;
fields = fieldnames(options.optimizer);
for k = 1:numel(fields)
    work.options.(fields{k}) = options.optimizer.(fields{k});
end
work.retain_jacobian = true;
% 原假设重估强径，为候选排序和嵌套模型比较提供共同基准。
nullFit = drdm.estimation.fit_paths(work,referencePaths);
scan = struct('kind','weak-add-one-scan-v1','options',options,'null_fit',nullFit, ...
    'candidates',struct([]),'coarse',struct([]),'trials',struct([]),'best_fit',[], ...
    'best_candidate_index',0,'weak_path',[],'updated_references',[], ...
    'best_assessment',[],'statistic',inf,'attempted_maximum',0, ...
    'complete',false,'reason','reference_fit_failed','context_id',drdm.statistics.conditional_signature(work));
if ~nullFit.converged || ~nullFit.identifiable
    return;
end
candidates = drdm.detection.propose_weak_candidates(work,nullFit,options);
prediction = drdm.model.predict_from_context(work,nullFit.paths);
whiteResidual = drdm.estimation.whiten(work.observed-prediction.normalized_vector,work.cholesky);
g = nullFit.whitened_jacobian;
if ~isempty(work.cholesky)
    g = work.cholesky*g;
end
coarse = repmat(struct('score',0,'gain',0,'estimable',false),1,numel(candidates));
for k = 1:numel(candidates)
    response = drdm.model.weak_path_response(work,nullFit.paths,candidates(k));
    geometry = drdm.statistics.profile_weak_response(response.linear,response.self,g,work.cholesky);
    if geometry.gain_estimable
        estimate = geometry.gain_operator*whiteResidual;
        coarse(k) = struct('score',sum((geometry.profiled_basis.'*whiteResidual).^2), ...
            'gain',estimate(1)+1i*estimate(2),'estimable',true);
    end
end
% 粗排只使用强径子空间之外的弱径响应，选定固定数量的精化初值。
eligible = find([coarse.estimable]);
[~,order] = sort([coarse(eligible).score],'descend');
selected = eligible(order(1:min(numel(order),options.refine_candidates)));
template = struct('candidate_index',0,'fit',[],'association',[], ...
    'eligible',false,'reason','','statistic',0,'assessment',[]);
trials = repmat(template,1,numel(selected));
samples = [];
if ~isempty(selected)
    samples = drdm.statistics.conditional_error_samples(work,nullFit.paths,work.uncertainty_draws);
end
% 备择模型联合重估全部强弱径，并检查路径对应和独立误差预算。
maximum = 0;
for k = 1:numel(selected)
    index = selected(k);
    c = candidates(index);
    seed = struct('tau_s',c.tau_s,'nu_hz',c.nu_hz,'gain_complex',coarse(index).gain);
    fitted = drdm.estimation.fit_paths(work,[nullFit.paths seed], ...
        [[nullFit.paths.gain_complex] seed.gain_complex]);
    trial = template;
    trial.candidate_index = index;
    trial.fit = fitted;
    trial.statistic = max(0,nullFit.cost-fitted.cost);
    trial.reason = 'optimization_failed';
    if fitted.converged
        trial.reason = 'unresolved_or_boundary';
        if fitted.identifiable && ~any(fitted.at_boundary) && ~drdm.estimation.has_close_paths(fitted.paths,work)
            association = drdm.estimation.match_reference_paths(nullFit.paths,fitted.paths, ...
                [work.delay_scale_s work.doppler_scale_hz],options.reference_match_radius);
            trial.association = association;
            trial.reason = 'reference_association_failed';
            if association.complete && abs(fitted.paths(association.weak_index).gain_complex)>0
                weak = fitted.paths(association.weak_index);
                trial.assessment = drdm.statistics.assess_weak_reference(work,nullFit,weak,samples, ...
                    weak.gain_complex,options.requirements,'joint');
                trial.reason = 'reference_budget_failed';
                trial.eligible = trial.assessment.usable_local_model && ...
                    (~options.require_budget || trial.assessment.budget.budget_satisfied);
                if trial.eligible
                    trial.reason = 'eligible';
                end
            end
        end
    end
    trials(k) = trial;
    if trial.eligible && trial.statistic>maximum
        maximum = trial.statistic;
        scan.best_fit = fitted;
        scan.best_candidate_index = index;
        scan.weak_path = fitted.paths(trial.association.weak_index);
        scan.updated_references = fitted.paths(trial.association.reference_indices);
        scan.best_assessment = trial.assessment;
    end
end
scan.candidates = candidates;
scan.coarse = coarse;
scan.trials = trials;
scan.attempted_maximum = maximum;
% 选定候选中任一优化未完成时，保留未完成状态供统计标定计数。
scan.complete = all(arrayfun(@(t)t.fit.converged,trials));
scan.statistic = maximum;
scan.reason = 'no_eligible_weak_path';
if ~scan.complete
    scan.statistic = inf;
    scan.reason = 'weak_optimization_incomplete';
elseif ~isempty(scan.best_fit)
    scan.reason = 'weak_candidate_available';
end
end
