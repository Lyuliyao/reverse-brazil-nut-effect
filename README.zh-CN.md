# 论文代码与必要数据

本库对应 **Transient Reverse Brazil Nut Effect in Particle-Laden Inclined Flows**。
最终版论文保存在 `paper/reverse_Brazil_nut_effect.pdf`，文件未经修改。

安装 Python 依赖后，可直接检查数据并重画论文图 6–10：

```bash
python -m pip install -r requirements.txt
python scripts/validate_bundle.py
python scripts/reproduce_figures.py
```

生成的 PDF 和 PNG 保存在 `outputs/figures/`。图 2–5 的 18 个实验图像面板
保存在 `data/experiments/frames/`，与最终版 PDF 中的面板逐一对应。

库中包括 15 个参数扫描算例、4 个机制对照算例的输入，以及重画论文所需的
平均高度曲线、力数据、流速剖面和粒子快照。原始大型轨迹、全部 CFD 场和
实验视频未放入本库。

各图对应的脚本和数据见 [英文说明](README.md)。重新运行 CFD–DEM 的安装、
算例准备和运行步骤见 [模拟说明](simulation/README.md)。重画已有结果无需
安装 CFD–DEM。仓库准备过程中没有重新运行物理模拟。
