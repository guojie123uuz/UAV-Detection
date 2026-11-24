function result = run_weak(block,prepared,strongThreshold,weakOptions,calibration)
%RUN_WEAK 运行强径筛选及一轮新增弱径检验。
%
% 输入：block、prepared、strongThreshold；weakOptions：可选配置；calibration：可选完整弱径搜索标定。
% 输出：result：路径、接纳与完成状态、复数残差、搜索记录和标定结果。
%
% 每个 CPI 至多新增一条弱径；省略 calibration 时用零门限输出开发评价结果。

if nargin<4
    weakOptions = struct();
end
if nargin<5
    calibration = [];
end
trace = drdm.pipeline.weak_trace(block,prepared,strongThreshold,weakOptions);
threshold = 0;
if ~isempty(calibration)
    assert(strcmp(calibration.kind,'drdm-weak-search-calibration-v1') && ...
        strcmp(calibration.procedure_id,trace.procedure_id) && ...
        strcmp(calibration.code,trace.code) && ...
        strcmp(calibration.id,drdm.statistics.fingerprint(rmfield(calibration,'id'))), ...
        'drdm:CalibrationMismatch','Weak-search calibration, processing code or observation mask differs.');
    threshold = calibration.threshold;
end
result = drdm.pipeline.select_weak_trace(trace,threshold);
result.calibration = calibration;
if ~isempty(calibration)
    result.status = 'weak_search_calibration_failed';
    if calibration.validated
        result.status = 'weak_search_validated_scope';
    end
end
end
