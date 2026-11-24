function source = draw_conditional_source(cellInfo,known,specification,payloadSeed)
%DRAW_CONDITIONAL_SOURCE 按指定 RMC 模型生成一次独立业务实现。
%
% 输入：cellInfo、known：确认资源；specification：业务模型；payloadSeed：种子。
% 输出：source：含保护样本的条件源波形。
%
% 已确认 RE 保持固定，unknown_power_scale 缩放其余 RE 的功率。

assert(strcmp(known.kind,'known_only'),'drdm:EstimationSource','Expected confirmed resources.');
assert(strcmp(specification.family,'assumed_lte_rmc'),'drdm:ConditionalFamily','Unsupported traffic family.');
validateattributes(specification.unknown_power_scale,{'numeric'},{'scalar','real','finite','nonnegative'});
assert(ischar(specification.condition_source) && ~isempty(specification.condition_source), ...
    'drdm:ConditionalAssumption','Document the source of the traffic assumption.');
guard = cellInfo.cpi_first_subframe-1;
assert(guard>=1 && cellInfo.total_subframes==cellInfo.cpi_num_subframes+2*guard, ...
    'drdm:ConditionalSpan','RMC conditional generator currently requires equal guard subframes.');
cfg = drdm.config.default_config();
cfg.lte = cellInfo.enb;
first = mod(cellInfo.enb.NFrame*10+cellInfo.enb.NSubframe+guard,10240);
cfg.lte.NFrame = floor(first/10);
cfg.lte.NSubframe = mod(first,10);
cfg.num_subframes = cellInfo.cpi_num_subframes;
cfg.guard_subframes = guard;
cfg.payload_seed = payloadSeed;
[sample,~] = drdm.lte.generate_lte_source(cfg);
assert(isequal(size(sample.grid),size(known.grid)) && ...
    max(abs(sample.grid(known.mask)-known.grid(known.mask)))<1e-12, ...
    'drdm:ConditionalResources','Assumed source conflicts with confirmed resources.');
% 仅调整未纳入已确认模板的资源功率。
grid = known.grid+sqrt(specification.unknown_power_scale)*(sample.grid-known.grid);
source = struct('iq',lteOFDMModulate(cellInfo.enb,grid,0),'fs_hz',known.fs_hz, ...
    'cpi_first_sample',known.cpi_first_sample,'cpi_num_samples',known.cpi_num_samples, ...
    'kind','conditional_sample','payload_seed',payloadSeed);
end
