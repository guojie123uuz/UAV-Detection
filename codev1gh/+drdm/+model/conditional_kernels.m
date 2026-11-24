function [dictionary,kernels] = conditional_kernels(context,coordinates)
%CONDITIONAL_KERNELS 计算共享业务条件下的单径均值和有序径对均值。
%
% 输入：context：条件估计上下文；coordinates：2×K 缩放时延/多普勒。
% 输出：dictionary：RDMp 单径字典；kernels：RDMy 径对核和噪声均值。
%
% 先对共享业务的样本乘积求平均，再进行多普勒积分。

model = context.conditional_model;
count = size(coordinates,2);
points = numel(context.observations.rdmp.complex_values);
ySize = size(context.observations.rdmy.complex_values);
dictionary = complex(zeros(points,count));
pairs = complex(zeros(prod(ySize),count,count));
if count==0
    kernels = struct('pairs',pairs,'noise_mean',context.conditional_noise_mean/context.y_scale,'map_size',ySize);
    return;
end
% RDMp 为线性算子，可以直接传播平均源波形。
source = model.source;
blocks = repmat(context.block,1,count);
for k = 1:count
    path = struct('tau_s',coordinates(1,k)*context.delay_scale_s, ...
        'nu_hz',coordinates(2,k)*context.doppler_scale_hz,'gain_complex',1);
    blocks(k).iq = drdm.simulation.apply_path(source,path,context.cfg.fractional_delay);
    response = drdm.rdm.compute_rdmp_re(blocks(k),context.known,context.cfg.rdmp);
    dictionary(:,k) = response.complex_values(:)/context.p_scale;
end
% RDMy 为二次算子，须保留各份业务及全部有序路径对。
draws = size(model.iq_ensemble,2);
n = numel(context.block.iq);
frequency = context.cfg.rdmy.doppler_hz(:);
assert(16*(n*count*draws+n*numel(frequency)+numel(model.iq_ensemble))<=model.max_working_bytes, ...
    'drdm:EnsembleMemory','Conditional pair-kernel workspace exceeds its configured memory budget.');
waveforms = complex(zeros(n,count,draws));
w = drdm.rdm.window_vector(context.cfg.rdmy.sample_window,n).*context.block.valid_mask;
for draw = 1:draws
    source.iq = model.iq_ensemble(:,draw);
    for k = 1:count
        path = struct('tau_s',coordinates(1,k)*context.delay_scale_s, ...
            'nu_hz',coordinates(2,k)*context.doppler_scale_hz,'gain_complex',1);
        waveforms(:,k,draw) = drdm.simulation.apply_path(source,path,context.cfg.fractional_delay).*w;
    end
end
% 先对同一业务实现中的路径样本乘积求平均，再做多普勒积分。
% 此顺序与逐份计算径对图后平均等价，可复用固定的多普勒混频矩阵。
mixer = context.conditional_mixer;
lags = context.cfg.rdmy.delay_samples;
for i = 1:count
    a = reshape(waveforms(:,i,:),n,draws);
    for j = 1:count
        b = reshape(waveforms(:,j,:),n,draws);
        values = complex(zeros(ySize));
        for lagIndex = 1:numel(lags)
            lag = lags(lagIndex);
            indices = (max(1,1+lag):min(n,n+lag))';
            product = mean(a(indices,:).*conj(b(indices-lag,:)),2);
            values(:,lagIndex) = (mixer(indices,:).'*product).* ...
                exp(1i*pi*frequency*lag/context.block.fs_hz)/(context.block.fs_hz*context.y_scale);
        end
        pairs(:,i,j) = values(:);
    end
end
kernels = struct('pairs',pairs,'noise_mean',context.conditional_noise_mean/context.y_scale,'map_size',ySize);
end
