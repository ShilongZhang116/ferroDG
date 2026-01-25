# ferroDG: Ferroptosis in Dentate Gyrus Neurogenesis

## 简介

ferroDG 是一个用于分析小鼠海马齿状回（Dentate Gyrus）神经发生过程中铁死亡（Ferroptosis）活性的R语言分析流程。该项目基于单细胞RNA测序数据，研究不同年龄阶段（Young vs Old）神经发生谱系中铁死亡相关基因的表达模式。

## 主要功能

- **数据预处理**：处理GSE233363数据集的单细胞RNA测序数据
- **细胞类型识别**：识别神经发生谱系细胞类型（qNSC → nIPC → Neuroblast → GC）
- **铁死亡评分**：计算铁死亡促进剂和抑制剂模块评分
- **可视化分析**：生成UMAP图、小提琴图、线图、火山图、热图等
- **轨迹分析**：分析神经发生轨迹上铁死亡活性的变化

## 项目结构

```
ferroDG/
├── README.md                 # 项目说明文档
├── LICENSE                   # 许可证文件
├── .gitignore               # Git忽略文件配置
├── R/                       # R脚本目录
│   ├── 01_data_preprocessing.R    # 数据预处理脚本
│   ├── 02_ferroptosis_analysis.R  # 铁死亡分析脚本
│   └── 03_visualization.R         # 可视化脚本
├── analysis/                # 历史分析脚本与Notebook（分阶段归档）
│   ├── stage-1/             # 第一阶段：Notebook与早期脚本
│   └── stage-2/             # 第二阶段：R脚本（分组与新细胞分组）
├── data/                    # 数据目录
│   └── README.md            # 数据说明文档
├── figures/                 # 生成的图形输出
│   ├── stage-1/             # 第一阶段图形
│   └── stage-2/             # 第二阶段图形
├── docs/                    # 文档目录
│   ├── methodology.md       # 方法说明
│   └── results.md           # 结果说明
└── gene_lists/              # 基因列表
    ├── ferroptosis_promoter.txt    # 铁死亡促进剂基因
    ├── ferroptosis_inhibitor.txt   # 铁死亡抑制剂基因
    └── ferroptosis_regulator.txt  # 铁死亡调节因子基因
```

## 环境要求

### R版本
- R >= 4.0.0

### 主要R包
```r
# 单细胞分析
Seurat >= 4.0.0

# 数据处理
dplyr >= 1.0.0
tidyr >= 1.0.0

# 可视化
ggplot2 >= 3.3.0
ggrepel >= 0.9.0
pheatmap >= 1.0.0

# 统计分析
# 根据需要添加
```

### 安装依赖
```r
# 安装CRAN包
install.packages(c("Seurat", "dplyr", "ggplot2", "ggrepel", "pheatmap"))

# 或者使用remotes安装GitHub开发版
# install.packages("remotes")
# remotes::install_github("satijalab/seurat")
```

## 快速开始

### 1. 克隆仓库
```bash
git clone https://github.com/yourusername/ferroDG.git
cd ferroDG
```

### 2. 准备数据
将GSE233363数据集的Seurat对象文件放置到`data/`目录下：
```
data/
└── GSE233363/
    └── Seurat_combined_with_Celltype_Article.rds
```

### 3. 运行分析流程

#### 方式一：按顺序运行脚本
```r
# 在R中运行
source("R/01_data_preprocessing.R")
source("R/02_ferroptosis_analysis.R")
source("R/03_visualization.R")
```

#### 方式二：使用Rscript命令行
```bash
Rscript R/01_data_preprocessing.R
Rscript R/02_ferroptosis_analysis.R
Rscript R/03_visualization.R
```

### 4. 查看结果
生成的图形文件将保存在`figures/`目录下，包括PNG和PDF格式。

## 分析流程说明

### Stage 1: 初步分析
- 全细胞UMAP可视化
- 神经发生细胞类型UMAP分析
- 细胞组成分析
- 铁死亡评分计算和可视化

### Stage 2: 深度分析
- 铁死亡促进剂/抑制剂模块评分
- 调节因子基因表达分析
- 年龄组间差异表达分析
- 神经发生轨迹上的铁死亡活性变化
- 相关性分析

## 历史分析代码

`analysis/` 中保留了按阶段整理的原始脚本和Notebook，便于追溯分析过程与图形生成逻辑。
主流程建议优先使用 `R/` 目录下的脚本。

## 铁死亡基因集

本项目使用的铁死亡相关基因集包括：

### 促进剂基因 (Promoter)
促进铁死亡发生的基因，包括：
- 铁代谢相关：TFRC, NCOA4, HMOX1, ACO1, IREB2
- 脂质过氧化相关：ACSL4, LPCAT3, ALOX5, ALOX12, ALOX15
- 其他调节因子：PTGS2, CHAC1, NOX1, NOX4等

### 抑制剂基因 (Inhibitor)
抑制铁死亡发生的基因，包括：
- 抗氧化系统：GPX4, SLC7A11, GCLC, GCLM, GSS
- 铁储存：FTH1, NFS1
- Nrf2通路：NFE2L2, KEAP1
- 其他保护因子：AKR1C1, AKR1C2, AKR1C3等

### 调节因子基因 (Regulator)
铁死亡的关键调节因子，包括：
- NQO1, VDAC2, TP53等

## 输出图形说明

### 主要图形类型
1. **UMAP图**：展示细胞聚类和细胞类型分布
2. **小提琴图**：展示不同细胞类型和年龄组的评分分布
3. **线图**：展示评分随年龄或细胞类型的变化趋势
4. **火山图**：展示差异表达基因
5. **热图**：展示基因表达模式的相关性
6. **轨迹图**：展示神经发生轨迹上的铁死亡活性变化

## 引用

如果您在研究中使用了本项目的代码或方法，请引用：

```bibtex
@software{ferroDG,
  title = { ... },
  author = { ... },
  year = {2026},
  url = {https://github.com/shilongzhang/ferroDG},
}
```

## 数据来源

本分析使用的数据集：
- **GSE233363**: 单细胞RNA测序数据，包含不同年龄阶段小鼠海马齿状回的细胞

## 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

## 贡献

欢迎提交Issue和Pull Request！

## 联系方式

如有问题或建议，请通过以下方式联系：
- 提交 GitHub Issue
- 发送邮件至：your.email@example.com

## 更新日志

### Version 1.0.0 (2025-12-25)
- 初始版本发布
- 完成Stage 1和Stage 2的分析流程
- 实现铁死亡评分计算和可视化

## 致谢

感谢所有为铁死亡研究和神经发生研究做出贡献的研究者。

## 注意事项

1. 本项目使用的基因符号已从人类转换为小鼠格式
2. 年龄分组为Young vs Old，其他时间点已被移除
3. 所有图形均以PNG和PDF双格式保存
4. 分析结果仅供参考，请结合实验验证

## 相关资源

- [Ferroptosis Database](http://www.zhounan.org/ferroptosis/)
- [Seurat Documentation](https://satijalab.org/seurat/)
- [GSE233363 on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363)


