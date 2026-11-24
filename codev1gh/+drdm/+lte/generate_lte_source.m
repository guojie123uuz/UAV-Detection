function [source,cellInfo] = generate_lte_source(cfg)
%GENERATE_LTE_SOURCE 生成标准 LTE 参考测量信道波形。
%
% 输入：cfg：含 LTE 参数、CPI 长度、保护子帧和 payload_seed 的配置。
% 输出：source：完整波形与资源栅格；cellInfo：资源重建所需小区信息。
%
% 同一信号的各条路径共享该源波形；功率保持工具箱原始尺度。

drdm.lte.validate_cell(cfg.lte);
validateattributes(cfg.num_subframes,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.guard_subframes,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.payload_seed,{'numeric'},{'scalar','integer','>=',0,'<',2^32});
total = cfg.num_subframes+2*cfg.guard_subframes;
assert(total<10240,'drdm:Duration','Source must span less than an SFN cycle.');
first = mod(cfg.lte.NFrame*10+cfg.lte.NSubframe-cfg.guard_subframes,10240);
enb = cfg.lte;
enb.NFrame = floor(first/10);
enb.NSubframe = mod(first,10);
layout = drdm.lte.symbol_layout(enb,total);
cellInfo = struct('version','drdm-cell-v1','enb',enb,'total_subframes',total, ...
    'cpi_first_subframe',cfg.guard_subframes+1,'cpi_num_subframes',cfg.num_subframes, ...
    'confirmation_source','simulation_configuration');
numbers = mod(first+(0:total-1),10240);
frames = unique(floor(numbers/10),'stable');
records = repmat(struct('frame_number',0,'frame_phase',0,'mib_bits',[], ...
    'crc_ok',true),numel(frames),1);
for k = 1:numel(frames)
    e = enb;
    e.NFrame = frames(k);
    records(k) = struct('frame_number',frames(k),'frame_phase',mod(frames(k),4), ...
        'mib_bits',lteMIB(e),'crc_ok',true);
end
cellInfo.pbch_records = records;
% 按标准带宽选择 RMC，使用独立业务比特产生连续资源栅格。
bandwidths = [6 15 25 50 75 100];
rmcs = {'R.4','R.5','R.6','R.7','R.8','R.9'};
rmc = lteRMCDL(rmcs{find(bandwidths==enb.NDLRB,1)});
fields = fieldnames(enb);
for k = 1:numel(fields)
    rmc.(fields{k}) = enb.(fields{k});
end
rmc.PDSCH.TxScheme = 'Port0';
rmc.PDSCH.NLayers = 1;
rmc.OCNGPDSCHEnable = 'Off';
rmc.OCNGPDCCHEnable = 'Off';
stream = RandStream('mt19937ar','Seed',cfg.payload_seed);
grid = complex(zeros(12*enb.NDLRB,14*total));
done = 0;
payloadBits = 0;
while done<total
    start = mod(first+done,10240);
    count = min(total-done,10240-start); % 每次工具箱调用的帧号保持在 SFN 周期内
    rmc.NFrame = floor(start/10);
    rmc.NSubframe = mod(start,10);
    rmc.TotSubframes = count;
    nbits = max(1,max(rmc.PDSCH.TrBlkSizes(:))*count);
    bits = randi(stream,[0 1],nbits,1);
    payloadBits = payloadBits+nbits;
    [~,part] = lteRMCDLTool(rmc,bits);
    grid(:,14*done+(1:14*count)) = part;
    done = done+count;
end
iq = lteOFDMModulate(enb,grid,0);
source = struct('iq',iq,'grid',grid,'fs_hz',layout.fs_hz, ...
    'cpi_first_sample',cfg.guard_subframes*layout.samples_per_subframe+1, ...
    'cpi_num_samples',cfg.num_subframes*layout.samples_per_subframe, ...
    'payload_seed',cfg.payload_seed,'payload_bit_count',payloadBits, ...
    'kind','oracle_full','amplitude_units','LTE Toolbox IFFT 1/Nfft');
end
