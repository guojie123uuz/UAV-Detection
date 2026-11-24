function result = process_cpi(block,cellInfo,cfg)
%PROCESS_CPI 构造双 RDM 并按配置执行初始路径搜索。
%
% 输入：block：同步 IQ；cellInfo：确认小区信息；cfg：观测与搜索配置。
% 输出：result：复数双图、可选 strong_search、参考系和处理时间。
%
% cfg.estimation.enabled 控制是否在双图之后执行初始搜索。

block = drdm.frontend.validate_block(block);
known = drdm.lte.build_known_resources(cellInfo,cfg.known_resources);
token = tic;
rdmp = drdm.rdm.compute_rdmp_re(block,known,cfg.rdmp);
rdmy = drdm.rdm.compute_rdmy(block,cfg.rdmy);
% 结果记录处理参数，业务和噪声种子保存在模拟输入中。
processingConfig = struct('version',cfg.version,'known_resources',{cfg.known_resources}, ...
    'rdmp',cfg.rdmp,'rdmy',cfg.rdmy);
result = struct('status','observations_only','version',cfg.version, ...
    'rdmp',rdmp,'rdmy',rdmy,'config',processingConfig,'cell_info',cellInfo, ...
    'time_center_s',block.time_center_s,'reference_id',block.reference_id, ...
    'sample_index0',block.sample_index0,'processing_seconds',toc(token));
if isfield(cfg,'estimation') && cfg.estimation.enabled
    result.strong_search = drdm.detection.extract_strong_paths(block,known,result,cfg);
    result.status = 'strong_search_uncalibrated';
    result.config.estimation = result.strong_search.options;
    result.config.fractional_delay = cfg.fractional_delay;
    result.processing_seconds = toc(token);
end
end
