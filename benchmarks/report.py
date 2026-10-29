#!/usr/bin/env python3
"""Regenerate PERFORMANCE.md from the retained measurements."""
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "benchmarks/results"
def read(name):
    return json.loads((RESULTS / (name + ".json")).read_text())

progress = read("final-progression")["medians"]
release = read("release-performance")["medians"]
simulation = read("release-simulation")
holdout = read("stats-holdout")
convergence = read("stats-convergence")
assert holdout["passed"] and convergence["passed"]
assert release["original"] / release["release"] >= 10

steps = [
    ("stage01", "保守磁场上界：跳过远处无共振区间"),
    ("stage02", "轨道段几何/旋转系数自适应四次插值"),
    ("stage03", "轨道误差 1e-8、步长上限 0.5r"),
    ("stage04", "继承偏振步长，Magnus 误差 1e-4"),
    ("stage05", "绝热耦合切换参数 η=1e-2"),
    ("stage08", "几何插值检查容差 1e-3"),
    ("stage09", "轨道误差 1e-6、步长上限 0.8r"),
    ("stage10", "保守磁场下界：跳过近星无共振区间"),
    ("stage11", "跨粒子复用工作区，逐粒子重设原种子"),
    ("stage12", "角向 twist 查表，插值误差界 1e-8"),
    ("stage13", "弱旋转区按需直接求值，避免构建多余缓存"),
    ("stage15a", "轨道误差调整到最终 1e-5"),
    ("stage18", "复用扫描缓冲区、预计算常量、去除重复边界"),
    ("stage19", "只保留速度端点及根合并处的积分切点"),
    ("stage20", "光学深度积分误差 1e-4"),
    ("stage21", "PowerLaw2 分布与共振 Jacobian 的解析消项"),
    ("stage23", "相同荷质比组分共享解析共振根"),
    ("stage24", "用 from_chars 解析表格，保留输入校验"),
    ("stage25", "单位向量横向范数使用平方和；极小值仍用 hypot"),
]
rows = []
previous = "stage00"
for stage, description in steps:
    before, after = progress[previous], progress[stage]
    rows.append(f"| {stage} | {description} | {before:.6f} | {after:.6f} | {before/after:.3f}× | {progress['original']/after:.3f}× |")
    previous = stage
rows.append(f"| release | GK15 小区间快速积分，误差不合格回退 QAGS | {release['stage25']:.6f} | {release['release']:.6f} | {release['stage25']/release['release']:.3f}× | {release['original']/release['release']:.3f}× |")

experiments = [
    ("stage06", "step06", "stage05", "stage06", "累计光学深度、正弦变量变换；更慢，撤回"),
    ("stage07", "step07", "stage05", "stage07", "累计光学深度、多项式变量变换；更慢，撤回"),
    ("stage14", "step14", "stage13", "stage14", "直接投影简化基底；更慢，撤回"),
    ("stage15", "step15", "stage14", "stage15", "含 stage14 的临时分支；参数修改以 stage15a 独立重测后保留"),
    ("stage16", "step16", "stage15", "stage16", "径向幂函数 Hermite 查表；收益低于波动，撤回"),
    ("stage17", "step15a-17", "stage15a", "stage17", "细分缓存复用父节点；额外复制成本使总时间增加，撤回"),
    ("stage22", "step22", "stage21", "stage22", "轨道误差 1e-4；精度更低且收益低于波动，撤回"),
    ("stage26", "step26-and-original", "stage25", "stage26", "预计算共振径向窗口；初测略快，统一复测变为 0.995×，最终撤回"),
    ("stage27 / final", "step27", "stage26", "stage27", "GK15 初始分支；移除 stage26 后，以 release 独立重测并保留"),
]
trial_rows = []
for stage, result, parent, candidate, reason in experiments:
    m = read(result)["medians"]
    trial_rows.append(f"| {stage} | {m[parent]:.6f} | {m[candidate]:.6f} | {m[parent]/m[candidate]:.3f}× | {reason} |")

metrics = {t["test"]: t for t in holdout["tests"]}
metric_rows = []
for field, label in [("escape_fraction","逃逸比例"),("scatter_mean","平均散射次数"),("energy_inf_mean","平均逃逸能量 / keV"),("muk_mean","平均 muk"),("Q_mean","平均 Q"),("U_mean","平均 U"),("V_mean","平均 V")]:
    t = metrics[field]
    metric_rows.append(f"| {label} | {t['baseline']:.8f} | {t['candidate']:.8f} | {t['difference']:+.8f} |")
energy_change = abs(metrics["energy_inf_mean"]["difference"])/abs(metrics["energy_inf_mean"]["baseline"])*100
loop = {stage: statistics.median(row["summary"]["loop_seconds"] for row in records) for stage,records in simulation["records"].items()}
repro_stages = " ".join(["original","stage00"]+[s for s,_ in steps]+["release"])

text = f"""# 单线程优化记录

日期：2026-10-30（Asia/Shanghai）。默认 10000 粒子、种子 42、原物理参数与磁场表。

最终端到端中位数：**{release['original']:.6f} 秒 → {release['release']:.6f} 秒，{release['original']/release['release']:.3f} 倍加速**。
原版与优化版均输出原格式粒子文本到 `/dev/null`；计时包含启动、读表、模拟、格式化与退出刷新。
每版本预热一次、测量五次，交替正序/逆序运行；全部模拟串行，无多线程、多进程并行或 GPU。
同用 GCC 16.2.0、C++23、`-O3 -march=native`，未启用 `-ffast-math`。

关闭输出后，模拟循环中位数为 **{loop['stage00']:.6f} 秒 → {loop['release']:.6f} 秒，{loop['stage00']/loop['release']:.3f} 倍**。
此指标排除读表，包含初始化粒子及其完整传播，证明收益并非来自省略粒子输出。

基线固定为提交 `70b7c2218d9520011f0bf02314d353a97837204e`。
`original` 是未修改的源码；`stage00` 仅加入 CLI、统计与编译期可关闭的诊断。
二者默认文本输出的 SHA-256 相同，见 [CLI/基线检查](benchmarks/results/cli-checks.json)。
环境、GSL 版本和输入磁场表 SHA-256 见 [environment.json](benchmarks/results/environment.json)。

## 每一步的总性能收益

“本步”=完成此前已保留优化的版本耗时 ÷ 加上本次优化的版本耗时；“累计”=同批次未修改原版耗时 ÷ 当前耗时。
前面的阶段使用 [统一复测](benchmarks/results/final-progression.json)，末行使用清理试验分支后的 [最终交替测速](benchmarks/results/release-performance.json)。
每一行的修改前后数据来自同一批测量。毫秒级的小差异会受系统调度影响，不应解读为单项收益保证。

| 阶段 | 优化 | 修改前 / 秒 | 修改后 / 秒 | 本步 | 累计 |
|---|---|---:|---:|---:|---:|
{chr(10).join(rows)}

所有保留阶段都通过默认 10000 粒子的统计初筛；最终版本另通过下面的独立大样本验收。
stage01 的默认输出还与基线逐字节相同。初次单项测量保存在 `benchmarks/results/step*.json`，统一复测没有覆盖这些原始记录。

## 未保留的尝试与试验分支

| 阶段 | 对照 / 秒 | 试验 / 秒 | 比值 | 处理 |
|---|---:|---:|---:|---|
{chr(10).join(trial_rows)}

因此最终实现没有采用累计光学深度反演、径向幂函数查表或预计算共振径向窗口。
它保留原 QAGS/一维求根路径，仅在 GK15 误差估计合格时使用更便宜的小区间积分。
每个尝试都有完整的可重建补丁，包括撤回的版本。

## 统计与数值验收

最终验证使用未参与调参的三组样本，每组 **30 × 10000 = 300000 粒子**：

- 原版：种子 `4000042..4300041`。
- 优化版：种子 `5000042..5300041`。
- 优化实现的严格精度模式：种子 `6000042..6300041`。

每组内及各组之间种子均不重叠；逐粒子检查种子、计数和有限值，模拟异常会使整个验证失败，不会跳过粒子。
完整粒子样本保留在 `build/perf/holdout-*.npz`（构建产物，不纳入 Git）。

原版与优化版的关键统计如下：

| 统计量 | 原版 | 优化版 | 差值（优化−原版） |
|---|---:|---:|---:|
{chr(10).join(metric_rows)}

平均能量相对差 **{energy_change:.4f}%**，低于 1%；逃逸/吸收比例、平均 muk、平均 Q/U/V 的绝对差均低于 0.01。
能量、方向及 Q/U/V 的双样本 KS 检验、散射次数与能量—方向联合直方图检验、25 个能量/方向箱中的条件平均偏振检验均通过。
稀疏箱按预先写定规则合并；对全部统计检验使用 Holm 校正，族显著性水平 1%。
数值有限性、能量正值、方向范围、Stokes 范数检查全部通过。

严格精度模式与优化默认值的独立 300000 粒子比较也通过同一套标准。
完整差值、标准误、分位数、检验 p 值及校正结果：
[独立原版对照](benchmarks/results/stats-holdout.json)、[严格精度收敛对照](benchmarks/results/stats-convergence.json)。

`make test` 覆盖磁场上下界、表格解析与损坏输入、共振根残差、窄速度支持、零电流、根合并、解析权重消项、GK15 多项式精度及奇性回退、已知 Stokes 旋转、轨道守恒量、缓存插值、再次散射、磁极及向内吸收、CLI 参数和输出格式。
`make build/main.1` 编译通过；没有运行其多线程入口。

## 最终数值配置与接口

| 参数 | 默认优化 | `--accuracy strict` |
|---|---:|---:|
| 轨道 atol / rtol | 1e-5 | 1e-10 |
| 轨道最大步长 / r | 0.8 | 0.1 |
| Magnus atol / rtol | 1e-4 | 1e-6 |
| η_couple | 1e-2 | 1e-3 |
| delta_stokes | 1e-3 | 1e-3 |
| 光学深度 atol / rtol | 1e-4 | 1e-8 |
| 速度切点数 / 组分 | 2 | 4 |
| 边界扫描点数 / 无共振采样数 | 4 / 4 | 4 / 4 |
| 几何系数插值 | 容差 1e-3，失败回退直接求值 | 关闭 |
| 小区间 GK15 | 误差合格才使用 | 关闭 |
| twist 插值 | 相对误差界 1e-8 | 关闭 |

保持默认粒子数、种子序列、物理参数和原文本输出格式。新增参数：
`--count N`、`--seed N`、`--output text|csv|none`、`--summary`、`--accuracy fast|strict`。
CSV 增加种子、结局和散射次数，便于统计；摘要写 stderr，不污染粒子输出。

`PhotonEvolution::Options` 集中提供轨道、偏振和几何缓存配置；`ResonantScattering::Options` 提供积分及扫描配置。
缓存属于单个传播对象，不跨散射事件复用旧几何。只插值连续几何系数；散射率仍按共振根和分布支持计算，积分仍在共振边界分段。
散射候选点的接受步骤用完整直接计算的 Basis。非 PowerLaw2 分布保留通用权重路径。
加速与统计结论针对本次默认场景；其他物理参数应重新运行对照，不能据此推断同样的倍数或误差。

## 复现

沿用仓库现有 Homebrew GCC/GSL/HDF5/HighFive 配置和 `table/bfield_t10.txt`。

```sh
make main test build/main.1
python3 -m venv build/perf/venv
build/perf/venv/bin/pip install -r benchmarks/requirements.txt
python3 benchmarks/run.py build original
python3 benchmarks/run.py build stage00
python3 benchmarks/run.py build stage25
python3 benchmarks/run.py build release
python3 benchmarks/run.py measure original stage25 release --name release-performance
```

各补丁是相对于固定基线的完整补丁，不依赖前一个阶段，`build` 会在 `build/perf/<阶段>` 中重建。
任意两步可通过 `python3 benchmarks/run.py measure <前一步> <后一步> --name <记录名>` 交替测量。
统一复测的阶段顺序：`{repro_stages}`（先逐个 build，再交给 measure）。

```sh
build/perf/venv/bin/python3 benchmarks/compare.py collect stage00 holdout-reference --groups 30 --seed 4000042
build/perf/venv/bin/python3 benchmarks/compare.py collect release holdout-fast --groups 30 --seed 5000042
build/perf/venv/bin/python3 benchmarks/compare.py compare holdout-reference holdout-fast --name stats-holdout
build/perf/venv/bin/python3 benchmarks/compare.py collect release holdout-strict --groups 30 --seed 6000042 --accuracy strict
build/perf/venv/bin/python3 benchmarks/compare.py compare holdout-strict holdout-fast --name stats-convergence
python3 benchmarks/run.py measure stage00 release --output none --name release-simulation
python3 benchmarks/check_cli.py build/main --reference-check
make diagnostics
build/diagnostic --output none --summary
python3 benchmarks/report.py
```

诊断计数仅在 `MAGMC_DIAGNOSTICS` 构建中启用，正式测速不包含其开销。
原始测量与阶段补丁分别保存在 `benchmarks/results/` 和 `benchmarks/stages/`。
"""
(ROOT / "PERFORMANCE.md").write_text(text)
print(f"Wrote PERFORMANCE.md: {release['original']/release['release']:.3f}x; both 300k validations passed")
