function uncertainty = bootstrap_refit(context,fit,payloadSeeds,noiseSeeds,options)
%BOOTSTRAP_REFIT 以参数 bootstrap 估计联合路径不确定度。
%
% 输入：context、fit：条件拟合；payloadSeeds、noiseSeeds：独立重采样；options：区间概率和对应门限。
% 输出：uncertainty：参数样本、协方差、偏差、基本 bootstrap 区间及残差。
%
% 每次重估全部路径，并按一一对应统一参数顺序。

drdm.statistics.check_seed_sets(context.conditional_model.payload_seeds, ...
    context.background.payload_seeds,context.background.noise_seeds,payloadSeeds,noiseSeeds);
assert(numel(payloadSeeds)==numel(noiseSeeds) && numel(payloadSeeds)>=2, ...
    'drdm:BootstrapSamples','Need at least two paired bootstrap draws.');
validateattributes(options.interval_probability,{'numeric'},{'scalar','real','>',0,'<',1});
validateattributes(options.bootstrap_match_radius,{'numeric'},{'scalar','real','positive','finite'});
count = numel(fit.paths);
draws = numel(payloadSeeds);
assert(count>=1 && count<=6,'drdm:BootstrapPathCount','Current exact label assignment supports 1 to 6 paths.');
assert(fit.converged && fit.identifiable,'drdm:BootstrapFit','Reference fit must converge and be identifiable.');
centre = drdm.estimation.path_vector(context,fit.paths);
parameters = nan(4*count,draws);
residuals = nan(numel(context.observed),draws);
valid = false(1,draws);
reasons = repmat({'not_run'},1,draws);
% 每份重采样都重新拟合全部路径，再按路径坐标统一参数次序。
for draw = 1:draws
    block = drdm.statistics.simulate_conditional_block(context,fit.paths,payloadSeeds(draw),noiseSeeds(draw));
    experiment = drdm.statistics.rebind_observation(context,block);
    repeated = drdm.estimation.fit_paths(experiment,fit.paths);
    if ~repeated.converged
        reasons{draw} = 'optimization_failed';
        continue;
    end
    if ~repeated.identifiable
        reasons{draw} = 'unresolved_cluster';
        continue;
    end
    if any(repeated.at_boundary)
        reasons{draw} = 'boundary';
        continue;
    end
    [ordered,matched] = match_paths(fit.paths,repeated.paths,context,options.bootstrap_match_radius);
    if ~matched
        reasons{draw} = 'assignment_ambiguous';
        continue;
    end
    parameters(:,draw) = drdm.estimation.path_vector(context,ordered);
    prediction = drdm.model.predict_from_context(experiment,ordered);
    residuals(:,draw) = experiment.observed-prediction.normalized_vector;
    valid(draw) = true;
    reasons{draw} = 'valid';
end
covariance = nan(4*count);
interval = nan(4*count,2);
bias = nan(4*count,1);
% 基本 bootstrap 区间将重复估计的分位数关于当前估计中心反射。
if sum(valid)>=2
    x = parameters(:,valid);
    bias = mean(x,2)-centre;
    dx = x-mean(x,2);
    covariance = dx*dx.'/(size(x,2)-1);
    tail = (1-options.interval_probability)/2;
    interval = [2*centre-empirical_quantile(x,1-tail),2*centre-empirical_quantile(x,tail)];
end
scale = [repmat([context.delay_scale_s;context.doppler_scale_hz],count,1);ones(2*count,1)];
tailResolved = sum(valid)*(1-options.interval_probability)/2>=1-1e-12;
status = 'conditional_local_intervals';
if ~all(valid)
    status = 'bootstrap_failures';
elseif ~tailResolved
    status = 'insufficient_tail_samples';
end
uncertainty = struct('status',status,'interval_probability',options.interval_probability, ...
    'reference_paths',fit.paths, ...
    'scaled_parameters',parameters,'scaled_centre',centre,'scaled_covariance',covariance, ...
    'scaled_bias',bias,'scaled_interval',interval,'physical_interval',interval.*scale, ...
    'physical_covariance',covariance.*(scale*scale.'),'physical_bias',bias.*scale, ...
    'halfwidth',max(abs(interval-centre),[],2).*scale,'valid',valid,'failure_reasons',{reasons}, ...
    'sample_count',draws,'valid_count',sum(valid),'tail_resolved',tailResolved, ...
    'refitted_residual_samples',residuals,'payload_seeds',payloadSeeds(:).','noise_seeds',noiseSeeds(:).', ...
    'procedure_id',drdm.statistics.conditional_signature(context),'coverage_validated',false, ...
    'scope','Basic parametric bootstrap, fixed path count/local basin/weights; selection coverage unvalidated.');
end

function q = empirical_quantile(x,p)
x = sort(x,2);
position = 1+(size(x,2)-1)*p;
left = floor(position);
right = ceil(position);
q = x(:,left)+(position-left)*(x(:,right)-x(:,left));
end

function [ordered,matched] = match_paths(reference,candidate,context,gate)
count = numel(reference);
orders = perms(1:count);
costs = zeros(size(orders,1),1);
distances = zeros(size(orders));
for k = 1:size(orders,1)
    paths = candidate(orders(k,:));
    distances(k,:) = hypot(([paths.tau_s]-[reference.tau_s])/context.delay_scale_s, ...
        ([paths.nu_hz]-[reference.nu_hz])/context.doppler_scale_hz);
    costs(k) = sum(distances(k,:).^2);
end
[sorted,index] = sort(costs);
chosen = index(1);
matched = all(distances(chosen,:)<gate) && (numel(sorted)==1 || sorted(2)-sorted(1)>1e-8);
ordered = candidate(orders(chosen,:));
end
