function reports = coverage_diagnostics(dataset)
%COVERAGE_DIAGNOSTICS 汇总路径恢复、参数误差及区间覆盖。
%
% 输入：dataset：run_suite 输出。
% 输出：reports：各类总数、恢复与覆盖数、点误差及区间中心误差统计。
%
% 误差均值和均方根基于完整恢复案例，同时报告全部案例及恢复案例数量。

profiles = dataset.suite.profiles;
template = struct('name','','sample_count',0,'recovered_count',0,'covered_count',0, ...
    'incomplete_count',0,'point_bias',[],'interval_centre_bias',[], ...
    'point_rmse',[],'interval_centre_rmse',[],'standardized_rmse',[], ...
    'component_noncoverage',[],'maximum_standardized_error',[]);
reports = repmat(template,1,numel(profiles));
for p = 1:numel(profiles)
    records = dataset.records([dataset.records.profile_index]==p);
    recovered = arrayfun(@(r)r.score.joint_recovery,records);
    reports(p).name = profiles(p).name;
    reports(p).sample_count = numel(records);
    reports(p).recovered_count = sum(recovered);
    reports(p).covered_count = sum(arrayfun(@(r)r.score.joint_interval_coverage,records));
    reports(p).incomplete_count = sum(arrayfun(@(r)~r.score.complete,records));
    if isempty(profiles(p).path_boxes) || ~any(recovered)
        continue;
    end
    scores = [records(recovered).score];
    pointErrors = [scores.point_error];
    centreErrors = [scores.parameter_error];
    scales = [scores.unit_halfwidth];
    standardized = centreErrors./max(scales,realmin);
    reports(p).point_bias = mean(pointErrors,2);
    reports(p).interval_centre_bias = mean(centreErrors,2);
    reports(p).point_rmse = sqrt(mean(pointErrors.^2,2));
    reports(p).interval_centre_rmse = sqrt(mean(centreErrors.^2,2));
    reports(p).standardized_rmse = sqrt(mean(standardized.^2,2));
    reports(p).component_noncoverage = sum(abs(standardized)>dataset.description.options.interval_multiplier,2);
    reports(p).maximum_standardized_error = max(abs(standardized),[],1);
end
end
