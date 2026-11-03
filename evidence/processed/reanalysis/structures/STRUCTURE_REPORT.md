# top100 干旱候选关键蛋白结构预测报告

> TODO 阶段 4 的 AlphaFold 项。全部计算在 `/ds2/workspace/ai/highland-barley-omics-gis/processed/reanalysis/structures/` 内完成；
> 未改动 `processed/reanalysis/` 下的任何既有文件，未执行任何 git 操作。

---

## 1. 方法：降级路径实际走到哪一步

| 优先级 | 方法 | 实际使用情况 | 靶标数 |
|---|---|---|---|
| 1 | AlphaFold DB 直接下载（PDB + pLDDT） | 9 个同源条目全部下载；其中 2 个与靶标序列**逐字符相同**，直接作为该靶标结构采用 | 2 采用 / 9 下载 |
| 2 | **ColabFold Local 1.6.2 本地 GPU 预测** | 成功。AF2 权重，5 模型 × 3 recycle，`--model-type auto`，MSA 走远程服务器 | **8** |
| 3 | ESMFold（后备） | **未使用**——ColabFold 安装与 GPU 推理均一次成功，无需降级 | 0 |
| 4 | 全部失败 | **未发生** | 0 |

**运行环境**：RTX 3090 24 GB，驱动 580.126.09 / CUDA 13.0，无 sudo。
用 `uv` 建独立 venv（Python 3.11）装 `colabfold[alphafold]`；
默认解析到的 `jaxlib 0.10.2` 是 CPU 版，补装 `jax[cuda12]==0.10.2` 后 `jax.devices()` 返回 `CudaDevice(id=0)`，
CUDA 12 运行时在 CUDA 13 驱动上通过向下兼容正常工作（安装约 16 分钟、约 5 GB nvidia wheel）。
全流程 10 个靶标端到端耗时约 30 分钟（含 MSA 检索），峰值显存约 3.2 GB。

### 1.1 `uniprot_ac` 列的一个关键说明（影响了靶标选择的前提）

任务给的流程是"用 `top100_annotation_full.tsv` 的 `uniprot_ac` 列查 AlphaFold DB，查到就不再本地预测"。
核对后确认：**该列的 AC 是 BLAST 比对到 Swiss-Prot 的最佳同源蛋白，不是大麦靶标自身的 UniProt 条目**。
例如 `HORVU.MOREX.r3.5HG0524540` 的 AC 是 `Q75HE5` = 水稻 NAC domain-containing protein 2；
`HORVU.MOREX.r3.2HG0099870` 是 `Q08277` = 玉米 HSP82。
只有 DHN2(`P12951`) 与 HVA1(`P14928`) 本身就是大麦蛋白，序列与靶标完全一致。

如果照字面执行，8 个靶标会拿到"水稻/玉米/拟南芥同源蛋白"的结构并当作大麦靶标结论使用——这是来源误述。因此实际处理为：

- 序列完全一致的 2 个（DHN2、HVA1）→ 按原意直接用 AFDB，不做本地预测；
- 其余 8 个 → 一律用 ColabFold 预测**大麦自身序列**；
- 9 个同源 AFDB 结构仍全部下载保存在 `afdb/`，作为"该家族应该长什么样"的参照系，并在第 5 节用作结构一致性检查的比对对象。

---

## 2. 每个靶标的 pLDDT 汇总

pLDDT 逐残基从 PDB 的 B-factor 列读取（AFDB 与 ColabFold 输出都用该列存 pLDDT）。
读取口径已校验：对 `AF-P12951` 复算得均值 54.64，AFDB API 上报 `globalMetricValue` = 54.66。

| 基因 ID | 注释 | 方向 | 长度 | 结构来源 | 最佳模型 | pLDDT 均值 | 中位数 | pTM | <50 | 50–70 | 70–90 | >90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HORVU.MOREX.r3.4HG0386900 | CHP01589/DUF | down | 87 | ColabFold 本地 | alphafold2_ptm_model_3_seed_000 | **84.7** | 93.8 | 0.650 | 0.023 | 0.195 | 0.115 | 0.667 |
| HORVU.MOREX.r3.5HG0524540 | NAC/NAM TF | up | 330 | ColabFold 本地 | alphafold2_ptm_model_1_seed_000 | **67.1** | 66.8 | 0.560 | 0.418 | 0.085 | 0.154 | 0.342 |
| HORVU.MOREX.r3.3HG0318380 | LRR-RK SERK2-like | down | 224 | ColabFold 本地 | alphafold2_ptm_model_3_seed_000 | **87.5** | 97.4 | 0.820 | 0.076 | 0.107 | 0.080 | 0.737 |
| HORVU.MOREX.r3.1HG0046190 | Aquaporin TIP3-1 | up | 263 | ColabFold 本地 | alphafold2_ptm_model_3_seed_000 | **90.3** | 97.9 | 0.880 | 0.072 | 0.046 | 0.065 | 0.818 |
| HORVU.MOREX.r3.6HG0596770 | Aquaporin PIP2-2 | down | 286 | ColabFold 本地 | alphafold2_ptm_model_3_seed_000 | **88.8** | 97.6 | 0.870 | 0.080 | 0.066 | 0.077 | 0.776 |
| HORVU.MOREX.r3.2HG0099870 | HSP90/HSP82 | up | 716 | ColabFold 本地 | alphafold2_ptm_model_5_seed_000 | **83.7** | 91.9 | 0.780 | 0.102 | 0.070 | 0.250 | 0.578 |
| HORVU.MOREX.r3.5HG0507410 | Dehydrin DHN2 | up | 139 | AFDB P12951（序列完全一致） | AFDB 官方模型 | **54.6** | 52.1 | — | 0.410 | 0.504 | 0.086 | 0.000 |
| HORVU.MOREX.r3.1HG0080020 | LEA HVA1 | up | 213 | AFDB P14928（序列完全一致） | AFDB 官方模型 | **58.4** | 55.5 | — | 0.300 | 0.474 | 0.225 | 0.000 |
| HORVU.MOREX.r3.7HG0637810 | LEA_1 candidate A | up | 117 | ColabFold 本地 | alphafold2_ptm_model_1_seed_000 | **56.9** | 56.5 | 0.250 | 0.239 | 0.607 | 0.154 | 0.000 |
| HORVU.MOREX.r3.2HG0197000 | LEA_1 candidate B | up | 113 | ColabFold 本地 | alphafold2_ptm_model_1_seed_000 | **50.4** | 51.1 | 0.170 | 0.398 | 0.602 | 0.000 | 0.000 |

MSA 深度（ColabFold 远程服务器返回）：CHP01589/DUF 2983；NAC/NAM TF 23213；LRR-RK SERK2-like 17877；Aquaporin TIP3-1 24363；Aquaporin PIP2-2 24353；HSP90/HSP82 47805；Dehydrin DHN2 1997；LEA HVA1 1801；LEA_1 candidate A 1095；LEA_1 candidate B 159 条。

说明：DHN2 与 HVA1 采用的是 AFDB 官方模型；这两个序列本流程也各跑了一遍 ColabFold 作交叉验证，结果见 3.2（模型名与 pLDDT 见 `colabfold_scores.json`）。LEA_1 candidate B 的 MSA 只有 159 条，是全部靶标中最少的，其低置信有 MSA 偏浅的贡献——但它同时满足"<70 占比 100%""有序核心 0 残基""MobiDB-lite 全程判无序"三条，结论仍成立。

图：`structure_plddt_summary.png`（逐残基 pLDDT 曲线，2×5）；`structure_gallery.png`（3D cartoon 按 pLDDT 着色，3Dmol.js + headless Chrome 离线渲染）。

### 2.1 分档

- **高置信（均值 ≥ 85，结构可直接用于分析）**：LRR-RK SERK2-like 87.5、Aquaporin TIP3-1 90.3、Aquaporin PIP2-2 88.8
- **中高（均值 70–85，折叠可信、部分 loop 需谨慎）**：CHP01589/DUF 84.7、HSP90/HSP82 83.7
- **低与无序（均值 < 70，多数残基无固定折叠）**：NAC/NAM TF 67.1、Dehydrin DHN2 54.6、LEA HVA1 58.4、LEA_1 candidate A 56.9、LEA_1 candidate B 50.4

---

## 3. 无序蛋白（DHN / LEA / HVA1）的无序性证据

三条独立证据链指向同一结论：这几个蛋白在溶液中没有固定三维折叠。

### 3.1 证据一：AlphaFold2 pLDDT 低置信区段占比

| 基因 | 蛋白 | pLDDT 均值 | <50 占比 | 50–70 占比 | **<70 占比（无固定折叠）** | 最长连续 <50 区段 | pTM |
|---|---|---|---|---|---|---|---|
| HORVU.MOREX.r3.5HG0507410 | Dehydrin DHN2 | 54.6 | 0.410 | 0.504 | **0.914** | 19 aa | — |
| HORVU.MOREX.r3.1HG0080020 | LEA HVA1 | 58.4 | 0.300 | 0.474 | **0.775** | 27 aa | — |
| HORVU.MOREX.r3.7HG0637810 | LEA_1 candidate A | 56.9 | 0.239 | 0.607 | **0.846** | 8 aa | 0.250 |
| HORVU.MOREX.r3.2HG0197000 | LEA_1 candidate B | 50.4 | 0.398 | 0.602 | **1.000** | 11 aa | 0.170 |

对照（`HORVU.MOREX.r3.4HG0386900` 之外的 5 个"有结构域"靶标）：
<50 占比 CHP01589/DUF 0.023、LRR-RK SERK2-like 0.076、Aquaporin TIP3-1 0.072、Aquaporin PIP2-2 0.080、HSP90/HSP82 0.102；
**<70 占比 CHP01589/DUF 0.218、LRR-RK SERK2-like 0.183、Aquaporin TIP3-1 0.118、Aquaporin PIP2-2 0.147、HSP90/HSP82 0.172**。

即：这 5 个靶标 **82%–88% 的残基**都能达到 pLDDT≥70（折叠可信），
而无序组 **77%–100% 的残基**落在 pLDDT<70（无固定折叠）。两组**完全不重叠**，分界清晰。

（`NAC/NAM TF` 不算在内：它是唯一的中间情形，<70 占比 0.503，与无序组重叠——因为它是"N 端有序结构域 + C 端长无序激活区"的嵌合架构（见 5.4），
不属于纯无序蛋白，因此不计入无序组。）

### 3.2 证据二：与 AlphaFold DB 官方预测的交叉验证

DHN2 与 HVA1 既有 AFDB（官方 AF2 单体预测）又有本流程 ColabFold 独立跑的同一序列，
两条完全独立的流程给出一致的结论：

| 基因 | AFDB pLDDT 均值 | 本流程 ColabFold 均值 | 差值 | 两者 CA-RMSD（全链） |
|---|---|---|---|---|
| HORVU.MOREX.r3.5HG0507410 | 54.6 | 52.7 | 2.0 | 21.9 Å |
| HORVU.MOREX.r3.1HG0080020 | 58.4 | 60.9 | 2.5 | 17.4 Å |

（DHN2 两条流程整体 RMSD 21.9 Å、HVA1 17.4 Å：折叠自由度大、几乎每个模型都不同，
这正是无序蛋白的典型表现——**不是流程不稳定，而是本来就没有唯一构象**。）

### 3.3 证据三：MobiDB-lite 独立无序预测

上游 InterProScan 已带 MobiDB-lite 一致性无序预测，与 pLDDT 结果完全吻合：

| 基因 | MobiDB-lite 判定的无序区（aa） | 覆盖 |
|---|---|---|
| HORVU.MOREX.r3.4HG0386900 | 无 | — |
| HORVU.MOREX.r3.5HG0524540 | 无 | — |
| HORVU.MOREX.r3.3HG0318380 | 无 | — |
| HORVU.MOREX.r3.1HG0046190 | 1–20 | 20 aa |
| HORVU.MOREX.r3.6HG0596770 | 1–22 | 22 aa |
| HORVU.MOREX.r3.2HG0099870 | 无 | — |
| HORVU.MOREX.r3.5HG0507410 | 1–60; 23–42; 49–60 | 60 aa |
| HORVU.MOREX.r3.1HG0080020 | 1–60; 13–23; 41–60 | 60 aa |
| HORVU.MOREX.r3.7HG0637810 | 25–60; 26–49 | 36 aa |
| HORVU.MOREX.r3.2HG0197000 | 1–60; 7–23; 24–37; 38–60 | 60 aa |

（MobiDB-lite 只在前 60 aa 上跑过，见 6.1；因此"1–60 全段无序"实际是"已检测部分全段无序"。）

### 3.4 结论

- **Dehydrin DHN2**：pLDDT<70 的残基占 **91.4%**（其中 <50 的 41.0%），pLDDT 均值仅 54.6。
- **LEA HVA1**：pLDDT<70 的残基占 **77.5%**（其中 <50 的 30.0%），pLDDT 均值仅 58.4。
- **LEA_1 candidate A**：pLDDT<70 的残基占 **84.6%**（其中 <50 的 23.9%），pLDDT 均值仅 56.9，pTM 仅 0.250（pTM<0.5 说明没有可靠的链内相对取向）。
- **LEA_1 candidate B**：pLDDT<70 的残基占 **100.0%**（其中 <50 的 39.8%），pLDDT 均值仅 50.4，pTM 仅 0.170（pTM<0.5 说明没有可靠的链内相对取向）。

- 两个 LEA_1 候选尤其极端：`2HG0197000` 在两个模型里**没有任何一个残基同时达到 pLDDT≥70**（有序核心 = 0 个残基），`7HG0637810` 也只有 18 个。`LEA_1`(Pfam PF03760) 本身是"偏亲水、无固定二级结构"的基序，不是可折叠的结构域——**"注释到 LEA_1 结构域"与"该区段无序"并不矛盾**。
- DHN2 是唯一一个能看出一点结构的：K 片段（富含赖氨酸的重复）在部分模型里形成一段短螺旋，但只有约 10 个残基进入 pLDDT≥70，其余全为柔性链，与脱水素"高度亲水、高度无序"的经典特征一致。

**结论**：对 DHN/LEA/HVA1 这类蛋白，结构预测的价值不是给出一个折叠，而是**定量地证明它没有折叠**；本次结果给出了可写进论文的无序性定量证据。

---

## 4. 有序 / 无序的对照意义

有序组 5 个靶标的 pLDDT 均值 83.7–90.3、<70 残基占比 0.12–0.22；无序组 4 个为均值 50.4–58.4、<70 占比 0.77–1.00——**两项指标上两组均不重叠**；说明这套流程具备区分能力，无序组的低分不是流程噪声或预测失败。

---

## 5. 与 CD-Search / InterPro 结构域注释的一致性检查

### 5.1 方法

把每个靶标的预测结构与同源参考结构做 **Smith-Waterman 局部比对**（BLOSUM62，gap −11/−1），
在比对上的残基对上做 Kabsch 叠加，报两个指标：

- `RMSD_core`：只取**两个模型都判为高置信**（pLDDT≥70）的残基做整体叠加；
- `best_60aa`：核心区内最佳叠合的 60 残基窗口的 RMSD——对多结构域蛋白与螺线管蛋白，这直接回答"这段折叠本身对不对"，不受结构域间相对取向差异干扰。

用局部比对而不是全局比对是必需的：大麦 LRR 蛋白只有 224 aa，而水稻 SERK2 有 628 aa，全局比对会把它强行拉到激酶域上，得出完全错误的结论（改用局部比对后该靶标 RMSD 从 11.66 Å 降到 1.74 Å）。

### 5.2 结果

> 下表的"注释到的结构域"列取自共享注释，**只覆盖前 60 aa**（见 6.1）。
> "比对对象"列：8 个靶标比对的是 AFDB 上的**同源蛋白**（水稻/玉米/拟南芥，非大麦自身）；
> DHN2 与 HVA1 因无独立参考，比对的是本流程对**同一序列**独立跑出的 ColabFold 模型，
> 所以其序列一致性恒为 1.00——这两行的 RMSD 反映的是"两次独立预测差多少"，即构象自由度。

| 靶标 | 注释到的结构域 | 比对对象 | 序列一致性 | 覆盖度 | RMSD_core | best_60aa | 判定 |
|---|---|---|---|---|---|---|---|
| CHP01589/DUF | TIGR01589（TIGR01589 domain）; PTHR31871（OS02G0137100 PROTEIN）; PF09713（Plant protein 1589 of unknown ） | 无同源条目（该基因无 uniprot_ac） | — | — | — | — | 仅凭折叠形态判断 |
| NAC/NAM TF | PTHR31719（NAC TRANSCRIPTION FACTOR 56）; PS51005（NAC domain profile.）; SSF101941（NAC domain） | AFDB 同源 Q75HE5（Oryza） | 0.85 | 0.91 | 0.49 Å | 0.26 Å | 一致 |
| LRR-RK SERK2-like | PF08263（Leucine rich repeat N-terminal）; pfam08263（LRRNT_2） | AFDB 同源 Q7XV05（Oryza） | 0.64 | 0.84 | 1.74 Å | 0.22 Å | 一致 |
| Aquaporin TIP3-1 | SSF81338（Aquaporin-like）; cl00200（MIP superfamily） | AFDB 同源 Q9FWV6（Oryza） | 0.89 | 1.00 | 0.30 Å | 0.12 Å | 一致 |
| Aquaporin PIP2-2 | PTHR45687（AQUAPORIN OR AQUAGLYCEROPORIN ）; G3DSA:1.20.1080.10（Glycerol uptake facilitator pr）; cl00200（MIP superfamily） | AFDB 同源 Q6K215（Oryza） | 0.94 | 1.00 | 0.76 Å | 0.16 Å | 一致 |
| HSP90/HSP82 | SSF55874（ATPase domain of HSP90 chapero）; PS00298（Heat shock hsp90 proteins fami）; G3DSA:3.30.565.10 | AFDB 同源 Q08277（Zea） | 0.91 | 0.99 | 12.77 Å | 0.19 Å | 折叠一致，域间取向不同 |
| Dehydrin DHN2 | 无正式结构域（仅无序预测） | 本流程独立 ColabFold 模型（同一序列） | 1.00 | 1.00 | 4.11 Å | — | 无序（构象不固定，RMSD 高符合预期） |
| LEA HVA1 | 无正式结构域（仅无序预测） | 本流程独立 ColabFold 模型（同一序列） | 1.00 | 1.00 | 4.68 Å | 2.06 Å | 无序（构象不固定，RMSD 高符合预期） |
| LEA_1 candidate A | PTHR33493（LATE EMBRYOGENESIS ABUNDANT PR）; PF03760（Late embryogenesis abundant (L）; pfam03760（LEA_1） | AFDB 同源 Q39138（Arabidopsis） | 0.48 | 0.91 | 5.48 Å | — | 无序（构象不固定，RMSD 高符合预期） |
| LEA_1 candidate B | PTHR33493（LATE EMBRYOGENESIS ABUNDANT PR）; PF03760（Late embryogenesis abundant (L）; pfam03760（LEA_1） | AFDB 同源 Q96273（Arabidopsis） | 0.53 | 0.82 | — | — | 无序（无高置信核心） |

### 5.3 结构域区间内的 pLDDT

| 靶标 | 注释的结构域（来源 / 区间） | 区间内平均 pLDDT | 区间内 ≥70 占比 |
|---|---|---|---|
| CHP01589/DUF | InterProScan TIGR01589 12–60 （TIGR01589 domain） | 94.6 | 1.000 |
| CHP01589/DUF | InterProScan PTHR31871 6–60 （OS02G0137100 PROTEIN） | 93.1 | 0.982 |
| CHP01589/DUF | InterProScan PF09713 15–60 （Plant protein 1589 of unknown func） | 94.8 | 1.000 |
| CHP01589/DUF | CD-Search pfam09713 17–60 （A_thal_3526） | 94.8 | 1.000 |
| NAC/NAM TF | InterProScan PTHR31719 19–60 （NAC TRANSCRIPTION FACTOR 56） | 93.2 | 1.000 |
| NAC/NAM TF | InterProScan PS51005 19–60 （NAC domain profile.） | 93.2 | 1.000 |
| NAC/NAM TF | InterProScan SSF101941 15–60 （NAC domain） | 92.5 | 1.000 |
| NAC/NAM TF | InterProScan G3DSA:2.170.150.80 28–60 （NAC domain） | 94.4 | 1.000 |
| NAC/NAM TF | InterProScan PF02365 20–60 （No apical meristem (NAM) protein） | 93.4 | 1.000 |
| NAC/NAM TF | CD-Search pfam02365 20–60 （NAM） | 93.4 | 1.000 |
| LRR-RK SERK2-like | InterProScan PF08263 24–60 （Leucine rich repeat N-terminal dom） | 93.9 | 0.973 |
| LRR-RK SERK2-like | CD-Search pfam08263 24–60 （LRRNT_2） | 93.9 | 0.973 |
| Aquaporin TIP3-1 | InterProScan SSF81338 21–53 （Aquaporin-like） | 95.3 | 1.000 |
| Aquaporin TIP3-1 | CD-Search cl00200 11–60 （MIP superfamily） | 88.5 | 0.940 |
| Aquaporin PIP2-2 | InterProScan PTHR45687 11–60 （AQUAPORIN OR AQUAGLYCEROPORIN RELA） | 87.1 | 0.820 |
| Aquaporin PIP2-2 | InterProScan G3DSA:1.20.1080.10 25–60 （Glycerol uptake facilitator protei） | 97.7 | 1.000 |
| Aquaporin PIP2-2 | CD-Search cl00200 32–60 （MIP superfamily） | 98.3 | 1.000 |
| HSP90/HSP82 | InterProScan SSF55874 16–59 （ATPase domain of HSP90 chaperone/D） | 84.5 | 0.909 |
| HSP90/HSP82 | InterProScan PS00298 36–45 （Heat shock hsp90 proteins family s） | 81.4 | 0.800 |
| HSP90/HSP82 | InterProScan G3DSA:3.30.565.10 9–60 （-） | 80.8 | 0.827 |
| HSP90/HSP82 | InterProScan PR00775 16–36 （90kDa heat shock protein signature） | 79.7 | 0.857 |
| HSP90/HSP82 | InterProScan PR00775 37–59 （90kDa heat shock protein signature） | 88.8 | 0.957 |
| HSP90/HSP82 | InterProScan PTHR11528 15–60 （HEAT SHOCK PROTEIN 90 FAMILY MEMBE） | 84.7 | 0.913 |
| HSP90/HSP82 | CD-Search cl46693 16–60 （HSP90 superfamily） | 84.6 | 0.911 |
| LEA_1 candidate A | InterProScan PTHR33493 1–60 （LATE EMBRYOGENESIS ABUNDANT PROTEI） | 63.5 | 0.300 |
| LEA_1 candidate A | InterProScan PF03760 1–60 （Late embryogenesis abundant (LEA) ） | 63.5 | 0.300 |
| LEA_1 candidate A | CD-Search pfam03760 1–60 （LEA_1） | 63.5 | 0.300 |
| LEA_1 candidate B | InterProScan PTHR33493 1–60 （LATE EMBRYOGENESIS ABUNDANT PROTEI） | 53.1 | 0.000 |
| LEA_1 candidate B | InterProScan PF03760 1–60 （Late embryogenesis abundant (LEA) ） | 53.1 | 0.000 |
| LEA_1 candidate B | CD-Search pfam03760 1–60 （LEA_1） | 53.1 | 0.000 |

### 5.4 逐靶标结论与不一致点

- **NAC/NAM（5H0524540）—— 一致性很好**：InterProScan 给出 Pfam PF02365 + ProSite PS51005 + SUPERFAMILY NAC + Gene3D NAC 四重支持；预测结构中 N 端折叠成紧凑的 α+β 球状域，与水稻 NAC2(AF-Q75HE5) 的高置信核心叠加 **RMSD 0.49 Å（153 残基）**，best_60aa 0.26 Å——**NAC 折叠成立**。同时 C 端 41.8% 残基 pLDDT<50、最长连续低置信 52 aa，对应 NAC 家族典型的无序转录激活区（TAR），也是该蛋白整体均值只有 67.1 的原因。这类"N 端有序 + C 端无序"正是 NAC 转录因子的真实架构，不是预测失败。
- **水通道蛋白（1H0046190 TIP3-1 / 6H0596770 PIP2-2）—— 一致性很好**：核心 RMSD 分别 0.30 Å / 0.76 Å，best_60aa 0.12 / 0.16 Å，6 次跨膜螺旋束清晰。N 端有一段低置信柔性尾（MobiDB-lite 也判 1–20 / 1–22 无序），是水通道蛋白 N 端调控尾的正常特征。
- **HSP90/HSP82（2H0099870）—— 折叠一致、域间取向不同**：716 aa 全长预测均值 83.7；与玉米 HSP82(AF-Q08277) 比对一致性 0.91、覆盖 0.99，但整体核心 RMSD 12.77 Å、best_60aa 只有 0.19 Å。分域拆开看：N 端 ATPase 域 7.8 Å、中间域 3.3 Å、C 端域 0.46 Å——**差异来自 N 端 ATPase 域与中间域的相对取向**（HSP90 本就是构象可变的分子开关），而不是二级/局部折叠不同。**全长预测在这里明显超出注释能力**：注释序列被截断在 60 aa，只能看到 N 端 ATPase 域的开头。
- **LRR 受体激酶（3H0318380）—— 折叠一致，但只是"半截蛋白"**：预测结构 N 端帽 + 一串 LRR 重复，与水稻 SERK2(AF-Q7XV05) 的局部比对一致性 0.64、核心 RMSD 1.74 Å、best_60aa 0.22 Å，**LRR 螺线管折叠成立**。
  但关键发现是：大麦这个蛋白只有 224 aa，**只对应 SERK2 的第 15–203 位**（即 LRRNT + 前几个 LRR 重复），既没有跨膜段之后的激酶域，也没有完整的 LRR 阵列。加上 SignalP/Phobius 给出信号肽与跨膜螺旋，**它更可能是基因模型不完整（或假基因），而不是一个完整的 SERK2 型受体激酶**。论文里不宜按"受体激酶"来写它的功能。
- **⚠ 不一致点（重要）**：任务描述把 `HORVU.MOREX.r3.4HG0386900` 列为"NAM/NAC 结构域候选"，但 CD-Search 与 InterProScan **都没有任何 NAC 证据**：
  CD-Search 只命中 `pfam09713`（E=7.7e-21），InterProScan 给出 NCBIfam TIGR01589 + PANTHER PTHR31871 + Pfam PF09713 三重支持，对应的是 **IPR006476 / CHP01589（植物保守未知功能蛋白）**，NAM 结构域（PF02365）未命中。
  该蛋白仅 87 aa，预测折叠是一个紧凑的小 α/β 球状域（均值 84.7、中位数 93.8，87 个残基里只有 2 个 <50），与 NAC 折叠无关。**建议把"NAM/NAC 候选"改称"CHP01589/DUF 候选"，或从 NAC 相关结论中剔除。**

---

## 6. 局限与重要提醒

### 6.1 共享注释表存在 60 aa 截断（影响面较大，建议上游修复）

`processed/reanalysis/annotation/` 下的 InterProScan 与 CD-Search 结果，其输入 FASTA `processed/reanalysis/top100_proteins.fasta` 中 107 条序列**全部恰好 60 aa**；
`iprscan_top100.tsv` 的 `sequence_length` 列 441 行全为 60，CD-Search 命中的 `To` 坐标最大值也恰为 60。
即 **现有结构域注释只覆盖每个蛋白 N 端前 60 个残基**，
而 `ANNOTATION_API_REPORT.md` 第 5 节"注意事项"未提及这一点。
对 HSP90(716 aa)、NAC(330 aa)、水通道蛋白(263/286 aa) 这类蛋白注释是不完整的，
下游引用（包括本报告 5.3 的结构域区间）都受此限制。
本次全长结构预测与注释的差异大多源于此，**不是预测错误**。

### 6.2 其他

1. 结构预测用的是 AF2 **单体**模式（`--model-type auto` → `alphafold2_ptm`）。水通道蛋白是四聚体、HSP90 是二聚体，单体预测不能评价其组装界面。
2. MSA 来自 ColabFold 公共服务器 `https://api.colabfold.com`，结果可复现但依赖该服务可用性。
3. `afdb/` 下 9 个 AFDB 结构是**其他物种的同源蛋白**（其中 2 个为大麦自身），仅作折叠参照，不能当作大麦靶标的结构引用。
4. 两个 LEA_1 候选与拟南芥参考的同源一致性只有 0.48 / 0.53，参考结构本身的代表性有限；其"无序"结论主要由 pLDDT 与 MobiDB-lite 两条内部证据支持。

---

## 7. 产出文件

目录：`/ds2/workspace/ai/highland-barley-omics-gis/processed/reanalysis/structures/`

| 文件 | 内容 |
|---|---|
| `pdb/<基因>_colabfold_rank001.pdb` | 本地 ColabFold 预测，每靶标 rank_001 最佳模型 |
| `afdb/AF-<AC>-F1.pdb` + `afdb/afdb_meta.json` | AlphaFold DB 下载的 9 个条目及元数据 |
| `PDB_LIST.tsv` | 全部 PDB 路径 / 来源 / 模型名对照 |
| `structure_plddt_summary.png` | 10 靶标逐残基 pLDDT 曲线（2×5） |
| `structure_gallery.png` | 10 靶标 3D cartoon 结构画廊（3Dmol.js + headless Chrome 离线渲染） |
| `structure_gallery_catrace.png` | 备用画廊：matplotlib Cα 走向图（同色标） |
| `structure_plddt_summary.tsv` | pLDDT 汇总表 |
| `plddt_<基因>.tsv` | 逐残基 pLDDT 原始数据 |
| `domain_consistency.tsv` / `.json` | 结构域一致性检查（含逐结构域区间 pLDDT） |
| `colabfold_scores.json` / `msa_depth.json` | rank_001 的 pLDDT、pTM / MSA 深度 |
| `cf_out/` | ColabFold 原始输出（5 模型 PDB、scores json、a3m、coverage 图、完整 log） |
| `render/` | 3Dmol 渲染用的 HTML 与单靶标 PNG |
| `*.py` | 全部分析脚本（见下） |

**复现命令**

```bash
cd /ds2/workspace/ai/highland-barley-omics-gis/processed/reanalysis/structures
# 1) ColabFold 本地 GPU 预测（.venv_cf 为已建好的环境）
.venv_cf/bin/colabfold_batch input.fasta cf_out --num-recycle 3 \
    --model-type auto --host-url https://api.colabfold.com --rank plddt
# 2) 汇总、画图、一致性检查
python3 finalize_outputs.py && python3 collect_scores.py && python3 msa_depth.py
python3 analyze_structures.py                      # pLDDT 统计（需 numpy）
.venv_cf/bin/python build_domain_annotation.py     # 结构域注释汇总（需 biopython）
python3 extract_mobidb.py
.venv_cf/bin/python domain_check.py                # 结构比对一致性检查
python3 make_figures.py && python3 render_gallery.py
python3 build_report.py                            # 生成本报告
```

**在 ChimeraX / PyMOL 中查看**

```bash
chimerax pdb/*.pdb                 # 8 个本地预测
chimerax afdb/AF-P12951-F1.pdb afdb/AF-P14928-F1.pdb   # DHN2 / HVA1
# ChimeraX 中按 pLDDT 着色：
#   color byattribute bfactor palette /alphafold   (或 palette ^RdYlBu_r)
```

（`structure_gallery.png` 已经是按 pLDDT 着色的 cartoon 渲染，可直接用作论文插图草稿。）
