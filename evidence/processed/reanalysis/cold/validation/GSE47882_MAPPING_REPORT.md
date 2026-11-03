# GSE47882 (GPL17286) 探针→MorexV3 基因映射与干旱共识基因方向一致性验证

数据集：GSE47882（大麦冷驯化时间序列 d0/d1/d7/d21）
平台：GPL17286 = "IPK Hordeum vulgare 12k SEED2"，spotted cDNA 芯片（RZPD/IPK Gatersleben）
目的：为该平台建立探针→MorexV3 基因映射，并在**基因水平**（而非探针水平）重做与干旱共识基因的方向一致性验证。

---

## 1. 结论摘要

1. **克隆序列成功找回**：GPL17286 平台注释只有 3 列（ID/POSITION/CLONE_ID），不含序列。通过 NCBI nuccore 的 dbEST 记录，**12049 个唯一克隆中拿到 12043 个序列（99.95%）**，仅 6 个克隆（对应 12 个探针）缺失。
2. **映射成功**：24172 个探针中 **16330 个（67.6%）映射到唯一 MorexV3 基因**，覆盖 6775 个基因。
3. **验证通过**：三个冷驯化对照（d1/d7/d21 vs d0）与 1197 个干旱共识基因的交集均显示**极显著的同向性**（127/151、149/193、137/191，二项检验 p 分别为 3.97e-18、1.52e-14、1.66e-09）。这与 GSE41513（GPL16130，叶片 6H）此前的结果方向一致，说明冷驯化与干旱胁迫在转录组层面共享大量同向调控基因。
4. 该结论对映射阈值不敏感（见 §7）。

---

## 2. 输入与产出文件

| 文件 | 说明 |
|---|---|
| `/ds2/.../processed/metadata/GPL17286.tsv` | 平台注释（输入） |
| `/ds2/.../processed/reanalysis/cold/GSE47882_{d1,d7,d21}_vs_d0_limma.tsv` | 探针水平 limma 结果（输入，首列为探针 ID） |
| `/ds2/.../processed/reanalysis/drought_consensus_genes.tsv` | 1197 个干旱共识基因（输入） |
| `/ds2/.../reference/MorexV3.cdna.fa` | MorexV3 cDNA 参考（37962 条转录本，含 UTR） |
| **`/ds2/.../processed/mapping/GPL17286_clones.fa`** | **产出**：12043 条克隆序列 |
| **`/ds2/.../processed/mapping/GPL17286_clone_accessions.tsv`** | **产出**：克隆→GenBank 登录号、序列长度 |
| **`/ds2/.../processed/mapping/GPL17286_clone_hits.tsv`** | **产出**：blastn 原始命中 |
| **`/ds2/.../processed/mapping/GPL17286_morex_mapping.tsv`** | **产出**：探针→基因映射表 |
| **`/ds2/.../processed/reanalysis/cold/validation/GSE47882_direction_validation.tsv`** | **产出**：方向一致性验证结果 |

脚本：`/ds1/workspace/ai/highland-barley-omics-gis/fetch_gpl17286_clones.py`、`map_gpl17286.py`、`cold_validation_gpl17286.py`
（未修改 `cold_validation.py` 本体，也未改动 GPL14927/GPL16130 的既有映射文件。）

---

## 3. 序列来源（关键难点）

### 3.1 命中情况

平台注释里的 `CLONE_ID`（如 `HU04A01`、`HA01A05`，共 12049 个，全部为 7 字符）是 EST 克隆名。
先用 20 个样例测试 `esearch db=nuccore term="HU04A01"[All Fields]` → **命中 1 条**：

```
>BQ663773.1 HU04A01u HU Hordeum vulgare subsp. vulgare cDNA clone HU04A01 3-PRIME, mRNA sequence
```

即这些克隆确实以 dbEST 记录存在于 NCBI nuccore，标题格式为
`<ACC> <库名> Hordeum vulgare subsp. vulgare cDNA clone <CLONE_ID> 3'/5', mRNA sequence`，
序列即克隆插入片段。**因此不需要走 ENA / CR-EST / 论文补充材料等备选途径。**

### 3.2 批量抓取策略

采用"按克隆名 OR 批量检索 + 直接 efetch fasta 并从 defline 反解克隆名"的两步法（省去 esummary）：

1. `esearch`（POST，避免超长 URL）：每批 100 个克隆名，`term = "C1"[Title] OR "C2"[Title] OR ...`，`retmax=100000` → 共 **121 批，得 23395 个 uid**（平均每克隆约 1.9 条记录，多数克隆同时有 5' 和 3' 两条 EST）。
2. `efetch`（POST，`rettype=fasta`）：每批 ≤500 个 uid → 共 **47 批**，返回 **23395 条记录**。
3. 本地用正则 `cDNA clone (\S+)` 从 defline 解析克隆名，并只保留 defline 含 `Hordeum vulgare` 的记录（排除同名但来自其他物种的记录，例如高粱的 `hw01a01` genomic survey sequence）。
4. 一个克隆多条记录时取**最长**的一条作为代表序列。

全程无 api_key，请求间 `sleep 0.4s`，符合 NCBI 3 请求/秒的限速要求；运行中遇到 1 次 SSL 中断，由重试逻辑自动恢复。

### 3.3 匹配率

| 指标 | 数值 |
|---|---|
| 平台唯一克隆数 | 12049 |
| **成功取得序列的克隆** | **12043（99.95%）** |
| 缺失克隆 | 6（`HT01M05`、`HV01F02`、`HY02D16`、`HY02F23`、`HY02H06`、`HY02J10`） |
| 受影响探针 | 12 |
| 序列长度 | 均值 521.7 bp，中位数 544 bp，范围 102–882 bp |

缺失的 6 个克隆在 NCBI 中以 `[All Fields]` 也检不到对应 cDNA 记录，属数据库本身缺失，非检索遗漏。

### 3.4 质量控制说明

首轮运行时复用了先前小规模试跑（pilot，300 个克隆）遗留的按**批次序号**命名的 efetch 缓存文件，而该缓存对应的 uid 列表与全量运行不同，导致 2 个批次（约 1000 条记录）被错误跳过，当时匹配率仅 11654/12049（96.7%）。
已改为**按 uid 批次内容哈希命名缓存**并重跑，得到上述 12043/12049。此处保留记录以便数值可追溯：若见到 96.7% 的中间数字，那是修复前的版本。

---

## 4. blastn 映射

```
blastn -query GPL17286_clones.fa -db MorexV3.cdna.fa \
       -num_threads 8 -evalue 1e-10 -max_target_seqs 10 -outfmt 6 <cols>
```

筛选条件：`pident ≥ 95` 且 `qcovhsp ≥ 95`（与 GPL14927/GPL16130 映射一致）。

| 层级 | 结果 |
|---|---|
| MorexV3 转录本 | 37962 |
| 原始命中 | 20984 条（10492 个克隆至少 1 条命中） |
| 通过 95/95 筛选的命中 | 9604 条（8300 个克隆） |
| 完全无命中（evalue 1e-10） | 1551 个克隆 |
| 有命中但未过筛 | 2192 个克隆（绝大多数是 qcovhsp 不足，pident 中位数 99） |

## 5. 映射状态统计

**探针水平**（`GPL17286_morex_mapping.tsv`，共 24172 行）

| status | 探针数 | 占比 |
|---|---|---|
| `unique_gene` | 16330 | 67.6% |
| `ambiguous` | 316 | 1.3% |
| `no_qualified_hit` | 7526 | 31.1% |

**克隆水平**（12049 个克隆）

| status | 克隆数 |
|---|---|
| `unique_gene` | 8142 |
| `ambiguous` | 158 |
| `no_qualified_hit` | 3749（含 6 个无序列克隆） |

16330 个 unique-gene 探针覆盖 **6775 个唯一 MorexV3 基因**（平均 2.4 个探针/基因，与每个克隆约重复点样 2 次相符）。
其中 **1278 个探针命中 531 个干旱共识基因**（占 1197 个共识基因的 44.4%），构成了后续验证的基础。

---

## 6. 基因水平方向一致性验证

方法：探针→基因后，取 `FDR<0.05` 且 `|logFC|≥1` 的基因；同一基因多个探针时取 `|logFC|` 最大的探针为代表；与 1197 个干旱共识基因求交集；方向一致数用 `scipy.stats.binomtest` 对 p=0.5 做双侧检验。

结果文件：`GSE47882_direction_validation.tsv`

| contrast | 显著探针 | 其中已映射 | 显著基因 | 与共识基因交集 | 同向 | 反向 | 同向上调 | 同向下调 | 二项检验 p |
|---|---|---|---|---|---|---|---|---|---|
| d1_vs_d0 | 2539 | 1831 | 1057 | 151 | **127** | 24 | 80 | 47 | **3.97e-18** |
| d7_vs_d0 | 4486 | 3208 | 1824 | 193 | **149** | 44 | 82 | 67 | **1.52e-14** |
| d21_vs_d0 | 3917 | 2773 | 1585 | 191 | **137** | 54 | 79 | 58 | **1.66e-09** |

三个时间点同向率分别为 84.1%、77.2%、71.7%，均极显著高于随机（50%）。
同向上调与同向下调基因数均明显多于反向，说明冷驯化与干旱胁迫的转录响应在**两个方向上**都高度共享，而非单向的胁迫激活偏好。

结论经独立实现复算复核，`overlap`/`same_direction`/`binom_p` 完全一致。

---

## 7. 稳健性检查

放宽 blastn 覆盖度阈值重做 d7 对照（其余条件不变）：

| qcovhsp 阈值 | 已映射探针 | 显著基因 | 交集 | 同向 | p |
|---|---|---|---|---|---|
| ≥95（正式采用） | 16330 | 1824 | 193 | 149 | 1.52e-14 |
| ≥80 | 17628 | 1945 | 203 | 157 | 2.20e-15 |
| ≥60 | 18350 | 1985 | 204 | 157 | 4.82e-15 |

交集与同向数变化很小，**方向一致性的结论不是严格阈值造成的假象**。

---

## 8. 失败与替代方案

### 8.1 未映射部分

- **6 个克隆（12 探针）无序列**：NCBI 无对应 EST 记录。备选途径未再深入，因为其对结果影响可忽略。若日后需要，可尝试 ENA portal API（`result=sequence`）或 IPK/CR-EST 存档（CR-EST 已下线，需经 web.archive.org）。
- **1551 个克隆在 evalue 1e-10 下完全无命中**：怀疑是 MorexV3 注释中不存在（或未注释）的基因、非编码/污染读段、或来自栽培种 Barke 的快速演化区段。
  用更灵敏的 `-task blastn`（非 megablast）重跑这 1551 个克隆：仅 210 个出现任何命中，且只有 3 个能通过 95/95 筛选——说明**这批克隆在 MorexV3 中确实缺乏近缘同源序列**，不是检索灵敏度问题。
- **2192 个克隆有命中但未过 95/95**：pident 中位数高达 99，卡在 `qcovhsp`。已排除 poly-A/T 尾（全库仅 29 条序列末端存在 ≥12 bp 的 A/T 连续串）和参考缺 UTR（MorexV3 cDNA 含 UTR，中位 266 bp）这两个可能原因；更可能是单次测序读段的低质量区、载体/嵌合片段或 3'UTR 在品种间的高度分化。若需提高覆盖度，可放宽 qcovhsp（见 §7，不影响结论）。

### 8.2 与既有脚本的关系

`cold_validation.py` 原本将 GSE47882 归为"无基因映射、只能做探针水平"的平台。本次映射建立后，GSE47882 已可像 GSE41513 一样做基因水平验证；`cold_validation.py` 本体未作改动，新逻辑独立放在 `cold_validation_gpl17286.py`，供后续合并参考。

---

## 9. 复现命令

```bash
python3 /ds1/workspace/ai/highland-barley-omics-gis/fetch_gpl17286_clones.py   # 抓取克隆序列（约 10 分钟，限速）
python3 /ds1/workspace/ai/highland-barley-omics-gis/map_gpl17286.py            # blastn 映射
python3 /ds1/workspace/ai/highland-barley-omics-gis/cold_validation_gpl17286.py # 方向一致性验证
```
