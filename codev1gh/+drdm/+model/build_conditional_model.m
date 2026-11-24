function model = build_conditional_model(cellInfo,known,specification,payloadSeeds,noisePower,maxBytes)
%BUILD_CONDITIONAL_MODEL 建立独立业务随机集合及条件均值模型。
%
% 输入：cellInfo、known：确认资源；specification：RMC/功率假设；payloadSeeds：种子；noisePower：噪声功率；maxBytes：可选内存限额。
% 输出：model：共享业务集合、平均波形及模型标识。
%
% 固定随机集合使连续参数拟合使用同一组蒙特卡洛积分点。

if nargin<6
    maxBytes = 256*1024^2;
end
drdm.statistics.check_seed_sets(payloadSeeds);
validateattributes(noisePower,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(maxBytes,{'numeric'},{'scalar','real','finite','positive'});
assert(numel(payloadSeeds)>=2,'drdm:EnsembleSize','Use at least two independent source draws.');
assert(16*numel(known.iq)*numel(payloadSeeds)<=maxBytes, ...
    'drdm:EnsembleMemory','Conditional waveform ensemble exceeds its explicit memory budget.');
iq = complex(zeros(numel(known.iq),numel(payloadSeeds)));
for k = 1:numel(payloadSeeds)
    source = drdm.model.draw_conditional_source(cellInfo,known,specification,payloadSeeds(k));
    iq(:,k) = source.iq;
end
description = struct('version','drdm-conditional-rmc-v1','cell_info',cellInfo, ...
    'known_resources',{known.resources},'specification',specification, ...
    'payload_seeds',payloadSeeds(:).','noise_power',noisePower);
model = struct('kind','conditional_ensemble','cell_info',cellInfo,'specification',specification, ...
    'source',rmfield(source,'payload_seed'),'iq_ensemble',iq,'mean_iq',mean(iq,2), ...
    'payload_seeds',payloadSeeds(:).','noise_power',noisePower,'max_working_bytes',maxBytes, ...
    'fingerprint',drdm.statistics.fingerprint(description), ...
    'scope','Finite Monte Carlo mean under declared traffic assumptions; not the transmitted realization.');
model.source.iq = model.mean_iq;
model.known_signature = drdm.statistics.fingerprint(struct('grid',known.grid,'fs_hz',known.fs_hz));
end
