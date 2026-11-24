function draws = cache_conditional_draws(cellInfo,known,specification,noisePower,payloadSeeds,noiseSeeds)
%CACHE_CONDITIONAL_DRAWS 缓存独立业务波形及接收噪声。
%
% 输入：cellInfo、known、specification：小区、确认资源与业务模型；noisePower：E∣n∣²；payloadSeeds、noiseSeeds：配对种子。
% 输出：draws：带保护区的 sources、CPI 内 noise 及种子记录。
%
% 矩阵每列对应一次独立实现；同次实现的所有路径共享源业务。

validateattributes(noisePower,{'numeric'},{'scalar','real','finite','nonnegative'});
drdm.statistics.check_seed_sets(payloadSeeds,noiseSeeds);
assert(numel(payloadSeeds)==numel(noiseSeeds) && numel(payloadSeeds)>=3, ...
    'drdm:SampleCount','Need at least three paired source/noise draws.');
count = numel(payloadSeeds);
sources = complex(zeros(numel(known.iq),count));
noise = complex(zeros(known.cpi_num_samples,count));
for k = 1:count
    source = drdm.model.draw_conditional_source(cellInfo,known,specification,payloadSeeds(k));
    sources(:,k) = source.iq;
    stream = RandStream('mt19937ar','Seed',noiseSeeds(k));
    noise(:,k) = sqrt(noisePower/2)*(randn(stream,known.cpi_num_samples,1)+ ...
        1i*randn(stream,known.cpi_num_samples,1));
end
draws = struct('sources',sources,'noise',noise, ...
    'payload_seeds',payloadSeeds(:).','noise_seeds',noiseSeeds(:).');
end
