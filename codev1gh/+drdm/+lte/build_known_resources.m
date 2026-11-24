function known = build_known_resources(cellInfo,resourceNames)
%BUILD_KNOWN_RESOURCES 由已确认信息重建 PSS、SSS、CRS 和 PBCH。
%
% 输入：cellInfo：小区及帧时序信息；resourceNames：可选资源名称元胞数组。
% 输出：known：已知资源栅格、标签、时域模板及符号坐标。
%
% PBCH 按 CRC 通过的 MIB 和四帧相位重建，其他 RE 取零。

if nargin<2
    resourceNames = {'pss','sss','crs','pbch'};
end
assert(iscellstr(resourceNames) && all(ismember(resourceNames,{'pss','sss','crs','pbch'})), ...
    'drdm:Resources','Unknown resource name.');
assert(numel(unique(resourceNames))==numel(resourceNames),'drdm:Resources','Duplicate resource name.');
enb = cellInfo.enb;
total = cellInfo.total_subframes;
layout = drdm.lte.symbol_layout(enb,total);
validateattributes(cellInfo.cpi_first_subframe,{'numeric'},{'scalar','integer','positive'});
validateattributes(cellInfo.cpi_num_subframes,{'numeric'},{'scalar','integer','positive'});
assert(cellInfo.cpi_first_subframe+cellInfo.cpi_num_subframes-1<=total, ...
    'drdm:CPIBounds','CPI exceeds resource span.');
% 逐子帧映射已确认资源，并以标签记录其类型。
grid = complex(zeros(enb.NDLRB*12,total*14));
labels = zeros(size(grid),'uint8');
first = enb.NFrame*10+enb.NSubframe;
for sf = 1:total
    number = mod(first+sf-1,10240);
    e = enb;
    e.NFrame = floor(number/10);
    e.NSubframe = mod(number,10);
    g = complex(zeros(enb.NDLRB*12,14));
    label = zeros(size(g),'uint8');
    if ismember('pss',resourceNames)
        idx = ltePSSIndices(e);
        g(idx) = ltePSS(e);
        label(idx) = 1;
    end
    if ismember('sss',resourceNames)
        idx = lteSSSIndices(e);
        g(idx) = lteSSS(e);
        label(idx) = 2;
    end
    if ismember('crs',resourceNames)
        idx = lteCellRSIndices(e);
        g(idx) = lteCellRS(e);
        label(idx) = 3;
    end
    if ismember('pbch',resourceNames) && e.NSubframe==0 && isfield(cellInfo,'pbch_records')
        records = cellInfo.pbch_records;
        selected = find([records.frame_number]==e.NFrame & [records.crc_ok]);
        assert(numel(selected)<=1,'drdm:PBCHRecord','Conflicting confirmed PBCH records.');
        if ~isempty(selected)
            record = records(selected);
            bits = record.mib_bits(:);
            assert(numel(bits)==24 && all(ismember(bits,[0 1])), ...
                'drdm:PBCHBits','A confirmed MIB must contain 24 bits.');
            decoded = lteMIB(bits);
            assert(record.frame_phase==mod(e.NFrame,4) && ...
                decoded.NFrame==4*floor(e.NFrame/4) && decoded.NDLRB==e.NDLRB && ...
                strcmp(decoded.Ng,e.Ng) && strcmp(decoded.PHICHDuration,e.PHICHDuration), ...
                'drdm:PBCHMetadata','Confirmed MIB or frame phase conflicts with cell metadata.');
            symbols = ltePBCH(e,lteBCH(e,bits));
            idx = ltePBCHIndices(e);
            rows = record.frame_phase*size(idx,1)+(1:size(idx,1));
            assert(~any(label(idx)),'drdm:ResourceCollision','PBCH overlaps another known RE.');
            g(idx) = symbols(rows,:);
            label(idx) = 4;
        end
    end
    columns = (sf-1)*14+(1:14);
    grid(:,columns) = g;
    labels(:,columns) = label;
end
% 截取 CPI 对应符号，保留模板前后保护区供有符号时延使用。
firstSample = (cellInfo.cpi_first_subframe-1)*layout.samples_per_subframe+1;
count = cellInfo.cpi_num_subframes*layout.samples_per_subframe;
symbols = (cellInfo.cpi_first_subframe-1)*14+(1:14*cellInfo.cpi_num_subframes);
known = struct('iq',lteOFDMModulate(enb,grid,0),'grid',grid, ...
    'resource_labels',labels,'mask',labels>0,'fs_hz',layout.fs_hz, ...
    'cpi_first_sample',firstSample,'cpi_num_samples',count, ...
    'observation_grid',grid(:,symbols),'observation_labels',labels(:,symbols), ...
    'useful_start_samples',layout.useful_start_samples(symbols)-firstSample+1, ...
    'cp_samples',layout.cp_samples(symbols),'subcarrier_hz',layout.subcarrier_hz, ...
    'fft_indices',layout.fft_indices,'nfft',layout.nfft,'kind','known_only', ...
    'resources',{resourceNames},'version','drdm-known-v1');
end
