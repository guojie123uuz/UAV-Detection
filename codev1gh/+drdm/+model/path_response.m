function response = path_response(known,block,path,cfg,mode)
%PATH_RESPONSE 计算单位增益路径的已知资源完整响应。
%
% 输入：known、block：模板及观测坐标；path：路径参数；cfg：传播/图配置；mode：ofdm_re 或 iq_reference。
% 输出：response：该路径在指定 RDMp 算子下的复数响应。
%
% 直接传播模板并调用观测算子，保留窗口、边界和相位影响。

if nargin<5
    mode = 'ofdm_re';
end
assert(strcmp(known.kind,'known_only'),'drdm:KernelSource','RDMp model requires a known-only template.');
path.gain_complex = 1;
synthetic = block;
synthetic.iq = drdm.simulation.apply_path(known,path,cfg.fractional_delay);
switch mode
    case 'ofdm_re'
        response = drdm.rdm.compute_rdmp_re(synthetic,known,cfg.rdmp);
    case 'iq_reference'
        response = drdm.rdm.compute_rdmp_iq(synthetic,known,cfg.rdmp);
    otherwise
        error('drdm:Mode','Unknown RDMp mode.');
end
response.metadata.unit_path = struct('tau_s',path.tau_s,'nu_hz',path.nu_hz);
response.metadata.source_kind = 'known_only';
end
