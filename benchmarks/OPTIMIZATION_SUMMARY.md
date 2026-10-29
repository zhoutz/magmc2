# 单线程优化全部尝试与加速倍数总结

记录日期：2026-10-30（Asia/Shanghai）。本文根据已落盘的测速、统计结果和阶段补丁整理；不补写未完成的验收，也不把不同批次的耗时混算为已验证加速比。

## 1. 当前结论

**尚未形成同时满足“至少 10 倍加速”和全部最终统计验收的版本。**

- `release` 在同一批交替测速中：原版 **6.769887 秒**，优化版 **0.652632 秒**，加速 **10.373 倍**。
- 该版本的独立 30 万粒子比较中，只有 V（圆偏振）的分布检验未通过；其他检验及关键总体量容差通过。
- 后续 `try_couple` 将绝热切换参数收紧到 0.005，诊断用的两组统计比较均通过，但尚未用新的独立样本、配对原版测速完成最终验收。
- 最后一批相对测速已落盘，见第 5 节。它没有同批未修改原版，不能据此宣称新候选已达到 10 倍。
- 工作区末尾还包含 `try_optical_build` 类改动；该试验没有测出净加速，也未完成统计验收。当前工作区不能直接视为已验收发布版本。

## 2. 测量口径与基线

- 基线提交：`70b7c2218d9520011f0bf02314d353a97837204e`；`original` 为未修改源码，`stage00` 为加入 CLI、统计和可关闭诊断的参考程序。
- 默认物理参数、10000 粒子、起始种子 42；默认粒子文本输出重定向到 `/dev/null`。耗时包含进程启动、读表、模拟、输出格式化与退出刷新。
- 每个版本预热一次，再测五次，取墙钟时间中位数；按版本列表正序、逆序交替运行。模拟逐个串行执行，不使用多线程、并发模拟或 GPU。
- GCC 16.2.0、C++23、`-O3 -march=native`，未使用 `-ffast-math`；macOS 15.7.3 arm64、GSL 2.8。
- **本步加速比 = 同批对照版本中位数 ÷ 试验版本中位数**。大于 1 表示更快，小于 1 表示更慢；例如 0.974× 表示速度降低，耗时增加约 2.7%。
- 不同批次受系统负载影响，不能直接相除或把各批的加速比连续相乘。第 3 节保留各次探索记录，第 4 节单列统一复测的累计收益。
- 初步基线三次为 6.520、6.534、6.496 秒，中位数 6.520 秒；它不是后续所有批次统一使用的分母。

环境和输入表哈希见 [environment.json](results/environment.json)。`original` 与 `stage00` 的默认文本输出 SHA-256 相同，见 [cli-checks.json](results/cli-checks.json)。

## 3. 按实施过程记录的全部主要尝试

“进入 release”仅描述代码归属，不代表最终统计验收通过。每个阶段名均对应 `stages/<阶段>.patch`；补丁相对于固定基线完整保存，不依赖前一个阶段。

| 阶段 | 对照 | 尝试内容 | 对照 / 秒 | 试验 / 秒 | 本步加速比 | 处理及说明 | 原始记录 |
|---|---|---|---:|---:|---:|---|---|
| [stage01](stages/stage01.patch) | stage00 | 保守磁场上界，跳过远处不可能发生共振的区间 | 6.474211 | 5.840680 | **1.108×** | 进入 release；默认 10000 粒子输出与基线逐字节一致 | [step01](results/step01.json) |
| [stage02](stages/stage02.patch) | stage01 | 轨道段内自适应四次插值，共享几何及旋转系数 | 5.869595 | 3.673543 | **1.598×** | 进入 release；不插值散射率本身 | [step02](results/step02.json) |
| [stage03](stages/stage03.patch) | stage02 | 轨道 atol/rtol 改为 1e-8，最大步长改为 0.5r | 3.815792 | 2.719819 | **1.403×** | 进入后续参数调优主线 | [step03](results/step03.json) |
| [stage04](stages/stage04.patch) | stage03 | 跨连续区间继承偏振步长；Magnus atol/rtol 改为 1e-4 | 2.677463 | 2.075169 | **1.290×** | 进入 release | [step04](results/step04.json) |
| [stage05](stages/stage05.patch) | stage04 | 绝热耦合切换参数 η 从 1e-3 改为 1e-2 | 2.049235 | 1.922090 | **1.066×** | 进入 release，但后续 V 分布失败促使重新检查此参数 | [step05](results/step05.json) |
| [stage06](stages/stage06.patch) | stage05 | 累计光学深度缓存及反演；正弦变量变换处理端点奇性 | 1.983841 | 2.167996 | **0.915×** | 变慢，撤回 | [step06](results/step06.json) |
| [stage07](stages/stage07.patch) | stage05 | 累计光学深度方案改用多项式变量变换 | 1.932308 | 1.960570 | **0.986×** | 仍比未采用该方案慢，撤回 | [step07](results/step07.json) |
| [stage08](stages/stage08.patch) | stage05 | 几何插值检查容差从 1e-5 改为 1e-3 | 1.895549 | 1.680078 | **1.128×** | 进入 release | [step08-10](results/step08-10.json) |
| [stage09](stages/stage09.patch) | stage08 | 轨道 atol/rtol 改为 1e-6，最大步长改为 0.8r | 1.680078 | 1.191827 | **1.410×** | 进入后续参数调优主线 | [step08-10](results/step08-10.json) |
| [stage10](stages/stage10.patch) | stage09 | 保守磁场下界，跳过近星强磁场无共振区间 | 1.191827 | 1.024662 | **1.163×** | 进入 release | [step08-10](results/step08-10.json) |
| [stage11](stages/stage11.patch) | stage10 | 复用 PhotonEvolution、积分及缓存工作区；每个粒子重设原种子 | 1.016289 | 0.971668 | **1.046×** | 进入 release；增加严格精度 CLI | [step11](results/step11.json) |
| [stage12](stages/stage12.patch) | stage11 | 预计算角向 twist；仅在相对插值误差界 ≤1e-8 时查表 | 0.969978 | 0.935869 | **1.036×** | 进入 release；否则直接计算幂函数 | [step12-13](results/step12-13.json) |
| [stage13](stages/stage13.patch) | stage12 | 弱旋转区直接求值，避免构建不划算的几何缓存 | 0.935869 | 0.924987 | **1.012×** | 进入 release | [step12-13](results/step12-13.json) |
| [stage14](stages/stage14.patch) | stage13 | 利用轴对称性，直接投影计算输运系数，简化完整基底 | 0.916611 | 0.926591 | **0.989×** | 变慢，撤回 | [step14](results/step14.json) |
| [stage15](stages/stage15.patch) | stage14 | 在含 stage14 的试验分支上，把轨道 atol/rtol 改为 1e-5 | 0.942899 | 0.872368 | **1.081×** | 参数有收益；随后移除 stage14，以 stage15a 独立重测 | [step15](results/step15.json) |
| [stage16](stages/stage16.patch) | stage15 | 径向幂函数使用 IEEE 尾数分段 Hermite 插值与指数表 | 0.853039 | 0.847989 | **1.006×** | 约 0.6% 的收益低于波动，撤回 | [step16](results/step16.json) |
| [stage15a](stages/stage15a.patch) | stage13 | 不含 stage14，仅把轨道 atol/rtol 改为 1e-5 | 0.924112 | 0.849643 | **1.088×** | 进入 release；这是保留主线上的独立比较 | [step15a-17](results/step15a-17.json) |
| [stage17](stages/stage17.patch) | stage15a | 几何缓存细分时复用父区间已有采样点 | 0.849643 | 0.862909 | **0.985×** | 额外复制与管理开销抵消节省，撤回 | [step15a-17](results/step15a-17.json) |
| [stage18](stages/stage18.patch) | stage15a | 复用共振扫描缓冲区、预计算边界常量、合并相同荷质比的重复边界 | 0.836948 | 0.785492 | **1.066×** | 进入 release | [step18](results/step18.json) |
| [stage19](stages/stage19.patch) | stage18 | 每组分速度切点数从 4 减为 2，保留速度端点及根合并边界 | 0.786002 | 0.752167 | **1.045×** | 进入 release；不删除物理支持边界 | [step19](results/step19.json) |
| [stage20](stages/stage20.patch) | stage19 | 光学深度积分 atol/rtol 从 1e-6 改为 1e-4 | 0.748636 | 0.729325 | **1.026×** | 进入 release | [step20](results/step20.json) |
| [stage21](stages/stage21.patch) | stage20 | PowerLaw2 分布与共振 Jacobian 解析消项，避免重复平方根 | 0.732465 | 0.688987 | **1.063×** | 进入 release；其他分布保留通用路径 | [step21](results/step21.json) |
| [stage22](stages/stage22.patch) | stage21 | 轨道 atol/rtol 进一步放宽至 1e-4 | 0.688751 | 0.684750 | **1.006×** | 精度降低而收益约 0.6%，撤回 | [step22](results/step22.json) |
| [stage23](stages/stage23.patch) | stage21 | 相邻、相同荷质比的组分共享解析共振根 | 0.772434 | 0.704881 | **1.096×** | 进入 release；分布权重和偏振重叠仍分别计算 | [step23](results/step23.json) |
| [stage24](stages/stage24.patch) | stage23 | 用 from_chars 读取磁场表格，保留输入校验 | 0.758588 | 0.742262 | **1.022×** | 进入 release | [step24-and-original](results/step24-and-original.json) |
| [stage25](stages/stage25.patch) | stage24 | 单位向量横向范数使用平方和开方；极小值仍回退 hypot | 0.733811 | 0.725637 | **1.011×** | 进入 release | [step25](results/step25.json) |
| [stage26](stages/stage26.patch) | stage25 | 把共振磁场上下界换算为每段飞行的径向窗口，避免每步计算幂函数 | 0.728388 | 0.719837 | **1.012×** | 初测略快，统一复测变慢约 0.5%；最终撤回 | [step26-and-original](results/step26-and-original.json) |
| [stage27](stages/stage27.patch) | stage26 | 先尝试内联 GK15 小区间积分；误差不合格时回退 QAGS | 0.722408 | 0.670801 | **1.077×** | 算法保留；剔除 stage26 后以 release 独立重测 | [step27](results/step27.json) |
| [release](stages/release.patch) | stage25 | 只在保留主线 stage25 上加入 GK15，并整理代码和零宽区间处理 | 0.696786 | 0.652632 | **1.068×** | 性能超过 10 倍，但独立大样本 V 分布检验失败 | [release-performance](results/release-performance.json) |

补充阶段与实现记录：

- `stage00` 只建立基准工具，不作为性能优化计入。可重建补丁为 [stage00.patch](stages/stage00.patch)。
- `final` 是 `stage27` 后的合并/格式整理快照，仍含后来撤回的 `stage26`，另加入零宽积分区间处理；不是新的独立优化。统一复测中为 0.638274 秒，对同批原版 6.685226 秒为 10.474×。该名称不表示验收完成。补丁见 [final.patch](stages/final.patch)。
- `release` 移除了 `stage26`，其 GK15 净收益重新相对 `stage25` 测量，因此没有把撤回的优化计入最终发布候选。
- `try_optical` 和 `try_optical_build` 的初次编译遇到依赖模板调用缺少 `template` 关键字的问题；修正后重新构建再测速，没有把编译失败的版本计入性能结果。

## 4. 统一复测的累计加速比

以下全部数字来自同一个 [final-progression.json](results/final-progression.json)，累计分母为该批未修改原版的中位数。它用于观察优化过程，不替代后来的统计验收或 `release` 的独立测速。

| 阶段 | 中位数 / 秒 | 相对同批原版 | 备注 |
|---|---:|---:|---|
| original | 6.685226 | 1.000× | 未修改原版 |
| stage00 | 6.711090 | 0.996× | 基准工具版本 |
| stage01 | 6.030298 | 1.109× | 当时的保留主线阶段 |
| stage02 | 3.829793 | 1.746× | 当时的保留主线阶段 |
| stage03 | 2.732004 | 2.447× | 当时的保留主线阶段 |
| stage04 | 2.117660 | 3.157× | 当时的保留主线阶段 |
| stage05 | 1.994583 | 3.352× | 当时的保留主线阶段 |
| stage08 | 1.733264 | 3.857× | 当时的保留主线阶段 |
| stage09 | 1.225336 | 5.456× | 当时的保留主线阶段 |
| stage10 | 1.056582 | 6.327× | 当时的保留主线阶段 |
| stage11 | 1.015723 | 6.582× | 当时的保留主线阶段 |
| stage12 | 0.981901 | 6.808× | 当时的保留主线阶段 |
| stage13 | 0.958397 | 6.975× | 当时的保留主线阶段 |
| stage15a | 0.882307 | 7.577× | 当时的保留主线阶段 |
| stage18 | 0.826435 | 8.089× | 当时的保留主线阶段 |
| stage19 | 0.783756 | 8.530× | 当时的保留主线阶段 |
| stage20 | 0.765023 | 8.739× | 当时的保留主线阶段 |
| stage21 | 0.725043 | 9.220× | 当时的保留主线阶段 |
| stage23 | 0.705968 | 9.470× | 当时的保留主线阶段 |
| stage24 | 0.695202 | 9.616× | 当时的保留主线阶段 |
| stage25 | 0.686392 | 9.740× | 当时的保留主线阶段 |
| stage26 | 0.689803 | 9.692× | 后来撤回 |
| final | 0.638274 | 10.474× | 旧合并快照，仍含 stage26 |

其中 `stage26` 在此批次相对 `stage25` 为 **0.995×**，与初测的轻微收益不一致，因此最终撤回。其余小于几个百分点的收益也应结合波动理解。

## 5. V 分布未通过后的精度诊断与新优化试验

最后一批测速已在中断前后完成并写入 [accuracy-and-optical.json](results/accuracy-and-optical.json)。所有数字均使用该批 `release` 或对应直接父版本作为对照，不套用其他批次的原版耗时。

| 试验 | 相对 release 的修改 | 中位数 / 秒 | 相对同批 release 的速度 | 相对直接父版本的速度 | 统计与结论 |
|---|---|---:|---:|---:|---|
| release | 同批对照 | 0.639294 | 1.000× | — | 已知独立样本 V_KS 未通过 |
| [try_magnus](stages/try_magnus.patch) | Magnus atol/rtol：1e-4 → 1e-5 | 0.660034 | **0.969×** | **0.969×**（release） | 原版对照仍因 V_KS 失败；同种子严格模式对照通过 |
| [try_couple](stages/try_couple.patch) | η_couple：0.01 → 0.005 | 0.656485 | **0.974×** | **0.974×**（release） | 两组诊断比较均通过；新独立最终验收尚未完成 |
| [try_cache](stages/try_cache.patch) | 几何插值容差：1e-3 → 1e-4 | 0.783276 | **0.816×** | **0.816×**（release） | 原版对照仍因 V_KS 失败；同种子严格模式对照通过 |
| [try_optical](stages/try_optical.patch) | 光学深度查询只插值四个需要的系数，不计算三个旋转系数 | 0.647837 | **0.987×** | **0.987×**（release） | 未测出加速；没有完成统计验证 |
| [try_optical_build](stages/try_optical_build.patch) | 在 try_optical 基础上，绝热状态下不以旋转系数误差触发缓存细分 | 0.644899 | **0.991×** | **1.005×**（try_optical） | 比 try_optical 略快，但仍慢于 release；没有完成统计验证 |

这些精度试验分别改变一个参数，没有把三个精度修正同时组合。当前记录也没有“收紧精度 + 新缓存改动”的组合版本完成验收的数据。

### 5.1 已完成的统计对照

- 原版诊断样本：`4000042..4300041`，30 万粒子。
- `release` 及三个精度试验：`5000042..5300041`，各 30 万粒子。
- 严格精度诊断也使用 `5000042..5300041`；与精度试验是同种子对照，不是另一组独立样本。
- 原计划另用 `6000042..6300041` 做最终严格模式验证，但首次独立验收失败后该流程中止；不能声称它已经完成。
- 这批数据随后用于选择修正方向，因此新的最终候选仍需保留新的独立验证样本。

采用总体量容差及分布检验：平均能量相对差 ≤1%；逃逸/吸收比例、平均 muk 和平均 Q/U/V 绝对差 ≤0.01；全部分布/均值检验使用 Holm 多重比较校正，族显著性水平 1%。

| 对照记录 | 总体结论 | V_KS 校正后 p 值 | 说明 |
|---|---|---:|---|
| [stats11-300k](results/stats11-300k.json) | 通过 | 1 | 较早 stage11 的独立样本检查 |
| [stats-holdout-attempt1](results/stats-holdout-attempt1.json) | 未通过 | 0.0056450505 | 首次 release 独立验收；stats-holdout.json 为同一结果，不是重复独立试验 |
| [stats-diagnose-strict](results/stats-diagnose-strict.json) | 通过 | 1 | 原版 4000042 起与严格模式 5000042 起 |
| [stats-try_magnus](results/stats-try_magnus.json) | 未通过 | 0.00730781014 | 原版与收紧 Magnus 精度 |
| [stats-try_couple](results/stats-try_couple.json) | 通过 | 0.27259357 | 原版与收紧绝热切换参数 |
| [stats-try_cache](results/stats-try_cache.json) | 未通过 | 0.00689209534 | 原版与收紧几何插值精度 |
| [stats-paired-try_magnus](results/stats-paired-try_magnus.json) | 通过 | 0.17623863 | 同种子严格模式与收紧 Magnus 精度 |
| [stats-paired-try_couple](results/stats-paired-try_couple.json) | 通过 | 1 | 同种子严格模式与收紧绝热切换参数 |
| [stats-paired-try_cache](results/stats-paired-try_cache.json) | 通过 | 0.155410834 | 同种子严格模式与收紧几何插值精度 |

首次失败的 V_KS 距离为 **0.005891266**，未校正 p 值为 **6.34275337e-05**，Holm 校正后为 **0.0056450505 < 0.01**。不能仅凭其他指标通过就忽略这项失败。

其余已完成的统计内容包括能谱、出射方向、Q/U/V 分布、散射次数、能量—方向联合直方图、分箱平均偏振、分位数，以及有限值和 Stokes 范数检查。多数保留阶段的 10000 粒子快速筛选也已通过，逐项结果在 `results/stats*-quick.json`；快速筛选不等同于最终大样本验收。

## 6. 初始性能诊断与验证范围

初始诊断构建的 10000 粒子结果如下。诊断计时包含计数/计时器开销，不用来计算正式加速比；数据见 [diagnostics-original.json](results/diagnostics-original.json)。

| 指标 | 初始值 |
|---|---:|
| 完整基底求值 | 53,420,274 |
| 磁场求值 | 53,420,274 |
| 共振区间扫描 | 1,931,110 |
| 一维求根调用 | 152,547 |
| 积分采样 | 12,761,595 |
| Magnus 接受步 | 3,004,805 |
| Magnus 拒绝步 | 854,273 |
| 轨道步 | 1,921,457 |
| 光学深度部分累计时间 | 3.349466 秒 |
| 偏振推进部分累计时间 | 2.711714 秒 |

已经运行通过的数值测试覆盖磁场界、解析共振根、窄速度支持、根合并、零电流、解析权重消项、GK15 精度和奇性回退、Stokes 旋转与范数、轨道守恒量、几何插值、再次散射、磁极与向内吸收、表格损坏输入及 CLI。另一个入口 `main.1` 已通过编译检查，没有运行其多线程模拟。

测试在此前 release 阶段通过；末尾新增的混合荷质比/通用分布测试，以及最新光学缓存改动，尚未完成新一轮完整测试。不能将旧测试结果直接当成当前整个工作区的验收。

## 7. 记录与复现

- [run.py](run.py)：从固定 Git 基线与完整阶段补丁重建，进行串行测速。
- [stages/](stages/)：包括撤回尝试在内的所有阶段补丁；`final`、`release` 是历史快照名，不是质量认证。
- [results/](results/)：每次重复计时、各批中位数、统计检验与环境元数据。
- [compare.py](compare.py)：采样与统计比较；NumPy/SciPy 版本在 [requirements.txt](requirements.txt)。
- [test_numerics.cpp](test_numerics.cpp)、[check_cli.py](check_cli.py)：数值与接口检查。
- 粒子样本位于 `build/perf/*.npz`，属于未纳入 Git 的构建产物；统计摘要保存在 `results/`。
- `report.py` 是此前准备的验收报告生成器，需要尚未完成的成功验收数据；本总结不依赖它，也没有使用它生成“全部通过”的报告。

例如复现一次已保存的性能比较（需要仓库现有磁场表、GCC 和 GSL 环境）：

```sh
python3 benchmarks/run.py build original
python3 benchmarks/run.py build stage25
python3 benchmarks/run.py build release
python3 benchmarks/run.py measure original stage25 release --name release-recheck
```

复现最后一批某项精度或缓存尝试：

```sh
python3 benchmarks/run.py build try_couple
python3 benchmarks/run.py measure release try_couple --name couple-recheck
```

最终仍待完成：选定修正后的代码与精度配置、验证实际工作区、使用新的独立样本通过全部统计门槛，并在同批未修改原版对照下重新确认 ≥10 倍加速。本文仅总结已完成工作，不继续实施优化。
