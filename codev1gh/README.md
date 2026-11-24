# 函数库

本目录包含 98 个 MATLAB 库函数，采用 `drdm` 包命名空间。函数名使用小写单词及下划线；各文件的中文帮助给出用途、输入输出、单位及关键计算约定。

## 目录

```text
codev1/
  +drdm/
    +config/       基础配置                         1
    +frontend/     IQ 数据封装与坐标检查            2
    +lte/          LTE 标准源及已确认资源重建        4
    +simulation/   多径传播与接收样本生成            2
    +rdm/          RDMp、RDMy 及相关计算             8
    +model/        路径、径对、条件均值与弱径增量    9
    +estimation/   候选生成、联合估计和路径对应     16
    +detection/    强径搜索与一轮新增弱径检验        5
    +statistics/   背景、区间、相关误差与参照预算   25
    +pipeline/     强径参照筛选与强弱径联合处理     15
    +evaluation/   路径评价、覆盖诊断和统计标定     10
    process_cpi.m  双图计算与可选初始搜索
```

主要入口：

- `drdm.process_cpi(block,cellInfo,cfg)`：构造复数双 RDM，可选执行初始路径搜索。
- `drdm.estimation.refine_conditional(block,cellInfo,seeds,cfg,options)`：固定路径数的条件双图精化与参数 bootstrap。
- `drdm.pipeline.prepare(cellInfo,cfg,options)` 与 `drdm.pipeline.run(block,prepared,certificate)`：顺序筛选或联合模型参照筛选；`certificate` 可省略。
- `drdm.pipeline.run_weak(block,prepared,strongThreshold,weakOptions,calibration)`：强径筛选及一轮新增弱径检验，返回联合重估路径和复数残差；后两个参数可省略。
- `drdm.statistics.assess_weak_reference`：根据候选实际响应、强径参数敏感度和独立误差样本计算局部参照预算。
- `drdm.evaluation.run_suite`、`calibrate_intervals`、`calibrate_pipeline`、`validate_pipeline`：完整流程的统计标定与独立评价。
- `drdm.evaluation.run_weak_trials` 与 `calibrate_weak_search`：在声明的强径零假设总体上训练弱径门限并独立评价。

例如，`help drdm.rdm.compute_rdmp_re` 可查看中文接口说明。

## 运行约定

运行环境为 MATLAB R2021a、LTE Toolbox；非线性参数估计使用 Optimization Toolbox。配置标识计算使用 MATLAB 自带 JVM。

当前 LTE 实现采用 FDD、普通循环前缀、单发射端口、单接收通道和零发射窗重叠；覆盖 6/15/25/50/75/100 RB。前端提供同步校正后的连续 IQ 及确认的小区信息。时间、频率分别采用秒和 Hz，路径复增益参考 CPI 中心；RDM 采用“多普勒行、时延列”。

条件统计配置明确给出 RMC 业务假设、噪声功率和独立随机种子。`options.reference_selection.mode` 支持 `sequential` 与 `joint_model`；`options.uncertainty.method` 支持信息矩阵近似 `information` 与独立误差传播 `independent_linearized`。弱径流程采用后两者中的 `joint_model` 和 `independent_linearized`，每个 CPI 最多使用五条强径参照、新增一条弱径。

`options` 由 `drdm.pipeline.default_options()` 生成，再按实际图域、精度与预算要求修改。弱径配置由 `drdm.detection.weak_options()` 生成。RDMy 的观测域应覆盖待研究的径对差分坐标；RDMp 的时延域须满足循环前缀与 FFT 窗支持范围。

`weak_accepted` 给出当前门限下的规则判决，`complete` 表示本次规定搜索是否完成，`status` 给出标定状态。省略弱径标定时采用零门限，供开发评价。配置与源码版本共同确定标定对象，本版通过重新 `prepare`、训练和独立评价生成对应结果。

通用条件仿真与统计标定函数属于库接口，随机场景由调用方通过 `suite`、`profile` 指定。
