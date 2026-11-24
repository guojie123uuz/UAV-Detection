function map = pack_map(values,delay,doppler,operator,energy,metadata)
%PACK_MAP 封装复数观测图、坐标和归一化信息。
%
% 输入：values：复数图；delay、doppler：坐标；operator：算子名；energy：能量；metadata：附加信息。
% 输出：map：保留原始复数值及 normalized_values 的结构。
%
% 归一化能量为零时，归一化值设为 NaN，原始复数值保留。

assert(isfinite(energy) && energy>=0,'drdm:Energy','Invalid normalization energy.');
map = struct('complex_values',values,'delay_s',delay(:).', ...
    'doppler_hz',doppler(:),'operator',operator,'operator_version','drdm-v1', ...
    'normalization_energy',energy,'metadata',metadata, ...
    'phase_reference','CPI centre','correlation_units','IQ squared seconds');
if energy>0
    map.normalized_values = values/energy;
else
    map.normalized_values = complex(nan(size(values)));
end
end
