<div align="center">

# Transient Reverse Brazil Nut Effect in Particle-Laden Inclined Flows

**Liyao Lyu**<sup>1</sup>, **Lingyun Ding**<sup>1,†</sup>, **Yifan Gu**<sup>1</sup>, **Yanming Ba**<sup>2</sup>, **Andrea Bertozzi**<sup>1</sup>

<sup>1</sup>加州大学洛杉矶分校数学系 &nbsp;·&nbsp; <sup>2</sup>哈维穆德学院数学系
<br><sup>†</sup>通讯作者

[English](README.md) | [简体中文](README.zh-CN.md)

</div>

本仓库是论文 *Transient Reverse Brazil Nut Effect in Particle-Laden Inclined
Flows* 的配套代码与数据，包括复现论文数值结果所需的分析脚本、CFD–DEM 算例输入、
处理后的数据，以及论文中使用的实验图像面板。

<p align="center"><img src="paper/figure_assets/figure06.jpg" width="92%"></p>

## 研究内容

在双分散颗粒悬浮液沿斜面流动的实验中，我们观察到瞬态的反巴西果效应（RBNE）：
颗粒竖直分布先进入沉降主导的 RBNE 阶段（小颗粒在上、大颗粒在下），随后经历混合阶段，
最终恢复为常规的巴西果效应（BNE，大颗粒在上）。我们将其解释为差异重力沉降比接触介导的
动力筛分（kinetic sieving）更早发挥作用而导致的瞬态过冲，并用 CFD–DEM 模拟和机制隔离
数值实验验证了这一解释。

竖直分层用平均高度差 Δz̄ = z̄_L − z̄_S 描述：Δz̄ < 0 为 RBNE，Δz̄ > 0 为 BNE。

## 快速开始

只需 Python ≥ 3.9，约一分钟即可重画论文图 6–10，无需安装 CFD–DEM 求解器：

```bash
python -m pip install -r requirements.txt
python scripts/validate_bundle.py      # 数据完整性与一致性检查
python scripts/reproduce_figures.py    # 图 6–10 输出到 outputs/figures/
```

也可使用 `make all`，或 `conda env create -f environment.yml`。

## 图与脚本对应关系

| 图 | 内容 | 脚本 |
| :---: | --- | --- |
| 1 | 斜面流示意图 | `paper/figure_assets/figure01.pdf` |
| 2–5 | 实验侧视与俯视观测 | `data/experiments/`（18 个面板） |
| 6 | 典型算例的瞬态分层与流体速度剖面 | `analysis/figure06.py` |
| 7 | 体积分数 φ 与倾角 α 的影响 | `analysis/figure07.py` |
| 8 | 粒径比的影响 | `analysis/figure08.py` |
| 9 | 两类颗粒的平均高度、接触力与曳力 | `analysis/figure09.py` |
| 10 | 机制隔离对照（仅沉降、仅剪切、等沉降速度） | `analysis/figure10.py` |

## 目录结构

- `analysis/`：作图脚本；`analysis/raw/` 为原始 CFD–DEM 输出的后处理工具
- `cases/`：15 个参数扫描算例与 4 个机制对照算例的输入
- `data/`：平均高度时间序列、粒子快照、流体剖面、机制对照曲线、实验图像
- `simulation/`：求解器补丁、算例准备与运行脚本（见 [simulation/README.md](simulation/README.md)）
- `paper/`：论文 PDF 与排版图；`provenance/`：校验和与数据来源记录

原始的 TB 级 CFD 场、完整粒子轨迹与实验视频不包含在本仓库中。

## 引用

```bibtex
@misc{lyu2026transient,
  title         = {Transient Reverse Brazil Nut Effect in Particle-Laden Inclined Flows},
  author        = {Lyu, Liyao and Ding, Lingyun and Gu, Yifan and Ba, Yanming and Bertozzi, Andrea},
  year          = {2026},
  eprint        = {XXXX.XXXXX},
  archivePrefix = {arXiv},
  note          = {arXiv identifier to be added}
}
```

## 许可

代码采用 [MIT 许可](LICENSE)；数据、图像与论文采用 [CC BY 4.0](LICENSE-DATA)；
`simulation/patches/` 中针对 OpenFOAM、CFDEMcoupling、LIGGGHTS 的补丁沿用上游的 GPL 许可。
