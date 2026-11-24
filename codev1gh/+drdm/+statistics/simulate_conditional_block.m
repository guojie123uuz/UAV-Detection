function block = simulate_conditional_block(context,paths,payloadSeed,noiseSeed)
%SIMULATE_CONDITIONAL_BLOCK 生成指定路径集合下的条件接收样本。
%
% 输入：context：条件模型与观测坐标；paths：路径；payloadSeed、noiseSeed：随机种子。
% 输出：block：同一观测域的新接收数据块。
%
% 一次抽样内所有路径共享一份源业务。

model = context.conditional_model;
drdm.statistics.check_seed_sets(model.payload_seeds,payloadSeed,noiseSeed);
source = drdm.model.draw_conditional_source(model.cell_info,context.known,model.specification,payloadSeed);
iq = complex(zeros(size(context.block.iq)));
for k = 1:numel(paths)
    iq = iq+drdm.simulation.apply_path(source,paths(k),context.cfg.fractional_delay);
end
stream = RandStream('mt19937ar','Seed',noiseSeed);
iq = iq+sqrt(model.noise_power/2)*(randn(stream,size(iq))+1i*randn(stream,size(iq)));
block = context.block;
block.iq = iq;
end
