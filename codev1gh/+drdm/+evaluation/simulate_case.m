function [block,truth] = simulate_case(prepared,profile,seeds)
%SIMULATE_CASE 从声明的参数分布生成独立接收样本。
%
% 输入：prepared：条件模型；profile：路径区间和功率比例；seeds：参数/业务/噪声三个种子。
% 输出：block：接收数据；truth：单独保存的路径真值和分布信息。
%
% 参数在各给定区间内独立均匀抽样。

drdm.statistics.check_seed_sets(prepared.internal_seeds,seeds);
assert(numel(seeds)==3,'drdm:SampleCount','Each case needs parameter, payload and noise seeds.');
boxes = profile.path_boxes;
validateattributes(boxes,{'numeric'},{'2d','real','finite','ncols',8});
assert(size(boxes,1)<=6 && all(all(boxes(:,1:2:end)<=boxes(:,2:2:end))) && ...
    all(boxes(:,5)>=0),'drdm:Scenario','Invalid path population bounds.');
validateattributes(profile.unknown_power_scale,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(profile.noise_power_scale,{'numeric'},{'scalar','real','finite','nonnegative'});
stream = RandStream('mt19937ar','Seed',seeds(1));
parameters = boxes(:,1:2:end)+(boxes(:,2:2:end)-boxes(:,1:2:end)).*rand(stream,size(boxes,1),4);
paths = repmat(struct('tau_s',0,'nu_hz',0,'gain_complex',0),1,size(boxes,1));
for k = 1:numel(paths)
    paths(k) = struct('tau_s',parameters(k,1),'nu_hz',parameters(k,2), ...
        'gain_complex',parameters(k,3)*exp(1i*parameters(k,4)));
end
specification = prepared.model.specification;
specification.unknown_power_scale = specification.unknown_power_scale*profile.unknown_power_scale;
source = drdm.model.draw_conditional_source(prepared.cell_info,prepared.known,specification,seeds(2));
iq = complex(zeros(prepared.known.cpi_num_samples,1));
for k = 1:numel(paths)
    iq = iq+drdm.simulation.apply_path(source,paths(k),prepared.cfg.fractional_delay);
end
stream = RandStream('mt19937ar','Seed',seeds(3));
noisePower = prepared.model.noise_power*profile.noise_power_scale;
iq = iq+sqrt(noisePower/2)*(randn(stream,size(iq))+1i*randn(stream,size(iq)));
block = drdm.frontend.make_block(iq,source.fs_hz,struct('reference_id','conditional_suite'));
truth = struct('paths',paths,'profile',profile.name,'seeds',seeds(:).', ...
    'source_assumption',specification,'noise_power',noisePower);
end
