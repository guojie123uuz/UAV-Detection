function [dictionary,waveforms] = unit_dictionary(context,coordinates)
%UNIT_DICTIONARY 构造固定观测坐标上的单位路径字典。
%
% 输入：context：拟合上下文；coordinates：2×K 缩放路径坐标。
% 输出：dictionary：RDMp 字典；waveforms：已知源路径波形或条件径对核。
%
% 条件模型分派到 conditional_kernels，保持统一的拟合接口。

if isfield(context,'conditional_model')
    [dictionary,waveforms] = drdm.model.conditional_kernels(context,coordinates);
    return;
end
count = size(coordinates,2);
points = numel(context.observations.rdmp.complex_values);
dictionary = complex(zeros(points,count));
if context.joint
    waveforms = complex(zeros(numel(context.block.iq),count));
else
    waveforms = [];
end
for k = 1:count
    path = struct('tau_s',coordinates(1,k)*context.delay_scale_s, ...
        'nu_hz',coordinates(2,k)*context.doppler_scale_hz,'gain_complex',1);
    synthetic = context.block;
    synthetic.iq = drdm.simulation.apply_path(context.known,path,context.cfg.fractional_delay);
    response = drdm.rdm.compute_rdmp_re(synthetic,context.known,context.cfg.rdmp);
    dictionary(:,k) = response.complex_values(:)/context.p_scale;
    if context.joint
        waveforms(:,k) = synthetic.iq;
    end
end
end
