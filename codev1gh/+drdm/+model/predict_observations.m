function prediction = predict_observations(known,block,paths,cfg,mode)
%PREDICT_OBSERVATIONS 按当前路径集合重建已知波形模型的双图响应。
%
% 输入：known、block、paths、cfg：模板、观测块、路径及配置；mode：模型名称。
% 输出：prediction：rdmp、可选 rdmy 和模型标识。
%
% rdmp_only 返回 RDMp；joint_known_only 对整个已知源生成双图预测。

if nargin<5
    mode = 'rdmp_only';
end
assert(strcmp(known.kind,'known_only'),'drdm:EstimationSource','Expected confirmed resources.');
assert(ismember(mode,{'rdmp_only','joint_known_only'}),'drdm:FitModel','Unknown prediction model.');
synthetic = block;
synthetic.iq = complex(zeros(size(block.iq)));
for k = 1:numel(paths)
    synthetic.iq = synthetic.iq+drdm.simulation.apply_path(known,paths(k),cfg.fractional_delay);
end
prediction = struct('rdmp',drdm.rdm.compute_rdmp_re(synthetic,known,cfg.rdmp), ...
    'rdmy',[],'model',mode,'source_kind','known_only');
if strcmp(mode,'joint_known_only')
    prediction.rdmy = drdm.rdm.compute_rdmy(synthetic,cfg.rdmy);
end
end
