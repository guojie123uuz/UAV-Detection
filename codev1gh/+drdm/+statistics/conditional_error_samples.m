function samples = conditional_error_samples(context,paths,draws)
%CONDITIONAL_ERROR_SAMPLES 生成指定路径下双图观测相对于条件均值的误差。
%
% 输入：context、paths：观测模型与路径；draws：缓存的独立业务和噪声。
% 输出：samples：实误差矩阵 errors、路径、样本数及处理标识。
%
% 每列为一次实现；沿用 context 中的归一化、共轭去重和预测模型。

count = size(draws.sources,2);
assert(count>=3 && size(draws.noise,2)==count && ...
    size(draws.sources,1)==numel(context.conditional_model.source.iq) && ...
    size(draws.noise,1)==numel(context.block.iq), ...
    'drdm:UncertaintySamples','Source and noise columns must match the CPI and source support.');
drdm.statistics.check_seed_sets(context.conditional_model.payload_seeds, ...
    draws.payload_seeds,draws.noise_seeds);
errors = zeros(numel(context.observed),count);
prediction = drdm.model.predict_from_context(context,paths);
source = context.conditional_model.source;
for k = 1:count
    source.iq = draws.sources(:,k);
    synthetic = context.block;
    synthetic.iq = draws.noise(:,k);
    for p = 1:numel(paths)
        synthetic.iq = synthetic.iq+drdm.simulation.apply_path(source,paths(p),context.cfg.fractional_delay);
    end
    rdmp = drdm.rdm.compute_rdmp_re(synthetic,context.known,context.cfg.rdmp);
    rdmy = drdm.rdm.compute_rdmy(synthetic,context.cfg.rdmy);
    observation = drdm.estimation.pack_observations(context, ...
        rdmp.complex_values/context.p_scale,rdmy.complex_values/context.y_scale);
    errors(:,k) = observation-prediction.normalized_vector;
end
samples = struct('errors',errors,'reference_paths',paths,'sample_count',count, ...
    'payload_seeds',draws.payload_seeds,'noise_seeds',draws.noise_seeds, ...
    'procedure_id',drdm.statistics.conditional_signature(context));
end
