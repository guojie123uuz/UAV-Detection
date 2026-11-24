function [vector,p,y] = model_vector(context,dictionary,waveforms,gains)
%MODEL_VECTOR 由单位响应和复增益组装固定观测模型。
%
% 输入：context：拟合模式与索引；dictionary：RDMp 字典；waveforms：波形或条件径对核；gains：复增益。
% 输出：vector：实模型向量；p、y：归一化复数预测。
%
% RDMp 对增益线性，RDMy 包含全部 a_i·conj(a_j) 二次项。

p = dictionary*gains;
y = [];
if context.joint
    if isfield(context,'conditional_model')
        y = waveforms.noise_mean;
        for i = 1:numel(gains)
            for j = 1:numel(gains)
                y = y+reshape(waveforms.pairs(:,i,j),waveforms.map_size)*gains(i)*conj(gains(j));
            end
        end
    else
        synthetic = context.block;
        synthetic.iq = waveforms*gains;
        response = drdm.rdm.compute_rdmy(synthetic,context.cfg.rdmy);
        y = response.complex_values/context.y_scale;
    end
end
vector = drdm.estimation.pack_observations(context,p,y);
end
