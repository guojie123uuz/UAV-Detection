function result = run(block,prepared,certificate)
%RUN 从同步 IQ 输出强径参照、参数区间和双图复数残差。
%
% 输入：block、prepared：接收数据及准备对象；certificate：可选完整流程统计验收结果。
% 输出：result：paths、selection、trace、prediction、残差及 reference_reliable。
%
% 统一入口支持 sequential 与 joint_model，标定对象与源码、配置及观测掩码绑定。

if nargin<3
    certificate = [];
end
assert(strcmp(prepared.description.code,drdm.pipeline.code_signature()), ...
    'drdm:CalibrationMismatch','Prepared processing code changed; rebuild and revalidate.');
if ~isempty(certificate)
    assert(isfield(certificate,'kind') && strcmp(certificate.kind,'drdm-complete-pipeline-certificate-v1'), ...
        'drdm:CalibrationMismatch','Expected a complete-pipeline certificate, not a single-round threshold.');
    assert(strcmp(certificate.prepared_id,prepared.procedure_id) && ...
        strcmp(certificate.certificate_id,drdm.statistics.fingerprint(rmfield(certificate,'certificate_id'))), ...
        'drdm:CalibrationMismatch','Certificate configuration or contents differ.');
end
trace = drdm.pipeline.trace(block,prepared);
threshold = max(prepared.options.minimum_threshold,prepared.options.development_threshold);
validated = false;
if ~isempty(certificate)
    % 源码和证书内容已在入口校验；此处核对当前采样域及掩码。
    assert(strcmp(certificate.procedure_id,trace.procedure_id), ...
        'drdm:CalibrationMismatch','Certificate does not match the processing domain or mask.');
    threshold = certificate.threshold;
    validated = certificate.validated;
end
selection = drdm.pipeline.select_trace(trace,threshold);
prediction = drdm.model.predict_from_context(trace.context,selection.paths);
residualP = trace.context.observations.rdmp;
residualY = trace.context.observations.rdmy;
residualP.complex_values = residualP.complex_values-prediction.rdmp.complex_values;
residualP.normalized_values = residualP.complex_values/trace.context.p_scale;
residualY.complex_values = residualY.complex_values-prediction.rdmy.complex_values;
residualY.normalized_values = residualY.complex_values/trace.context.y_scale;
residualY.normalization_energy = trace.context.y_scale;
reliable = validated && selection.complete && ~isempty(selection.paths);
selection.reference_reliable = reliable;
records = selection.paths;
for k = 1:numel(records)
    records(k).reference_reliable = reliable;
    records(k).role = 'conditional_development_reference';
    if reliable
        records(k).role = 'validated_conditional_reference';
    end
    records(k).reference_id = block.reference_id;
end
result = struct('status','sequential_conditional_unvalidated','selection',selection, ...
    'paths',records,'prediction',prediction,'rdmp_residual',residualP,'rdmy_residual',residualY, ...
    'trace',trace,'certificate',certificate,'reference_reliable',reliable, ...
    'scope','Any validation is restricted to the declared simulation population and frozen processing configuration.');
if validated
    result.status = 'sequential_conditional_validated_scope';
end
if strcmp(prepared.options.reference_selection.mode,'joint_model')
    result.status = 'joint_model_conditional_unvalidated';
    if validated
        result.status = 'joint_model_conditional_validated_scope';
    end
end
end
