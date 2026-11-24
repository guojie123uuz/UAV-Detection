function info = validate_cell(enb)
%VALIDATE_CELL 检查小区参数并取得 OFDM 调制参数。
%
% 输入：enb：LTE Toolbox 小区配置。
% 输出：info：FFT 长度、采样率及循环前缀长度。
%
% 采用 FDD、普通循环前缀、单发射端口及零发射窗重叠。

required = {'NDLRB','NCellID','CellRefP','CyclicPrefix','DuplexMode', ...
    'NFrame','NSubframe','Ng','PHICHDuration','Windowing'};
assert(isstruct(enb) && isscalar(enb) && all(isfield(enb,required)), ...
    'drdm:CellConfig','Incomplete LTE configuration.');
validateattributes(enb.NDLRB,{'numeric'},{'scalar','finite','integer'});
assert(ismember(enb.NDLRB,[6 15 25 50 75 100]),'drdm:Bandwidth','Unsupported bandwidth.');
validateattributes(enb.NCellID,{'numeric'},{'scalar','integer','>=',0,'<=',503});
validateattributes(enb.NFrame,{'numeric'},{'scalar','integer','>=',0,'<=',1023});
validateattributes(enb.NSubframe,{'numeric'},{'scalar','integer','>=',0,'<=',9});
assert(isequal(enb.CellRefP,1),'drdm:Ports','This release requires one transmit port.');
assert(strcmp(enb.DuplexMode,'FDD') && strcmp(enb.CyclicPrefix,'Normal'), ...
    'drdm:FrameType','This release requires FDD and normal CP.');
assert(isequal(enb.Windowing,0),'drdm:TxWindow','This release requires Windowing=0.');
assert(exist('lteOFDMInfo','file')==2,'drdm:LTEToolbox','LTE Toolbox is required.');
info = lteOFDMInfo(enb);
end
