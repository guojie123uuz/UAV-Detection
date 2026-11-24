function prediction = predict_from_context(context,paths)
%PREDICT_FROM_CONTEXT 使用拟合时的算子和尺度生成当前模型预测。
%
% 输入：context：固定观测域、归一化和模型；paths：路径集合。
% 输出：prediction：复数双图及 normalized_vector。
%
% 条件 RDMy 同时包括路径自项、交叉项和噪声均值。

coordinates = [[paths.tau_s]/context.delay_scale_s;[paths.nu_hz]/context.doppler_scale_hz];
[dictionary,kernels] = drdm.estimation.unit_dictionary(context,coordinates);
[vector,p,y] = drdm.estimation.model_vector(context,dictionary,kernels,reshape([paths.gain_complex],[],1));
pMap = context.observations.rdmp;
pMap.complex_values = reshape(p,size(pMap.complex_values))*context.p_scale;
pMap.normalized_values = pMap.complex_values/pMap.normalization_energy;
yMap = [];
if context.joint
    yMap = context.observations.rdmy;
    yMap.complex_values = y*context.y_scale;
    yMap.normalized_values = yMap.complex_values/context.y_scale;
    yMap.normalization_energy = context.y_scale;
end
prediction = struct('rdmp',pMap,'rdmy',yMap,'normalized_vector',vector,'model',context.options.fit_mode);
end
