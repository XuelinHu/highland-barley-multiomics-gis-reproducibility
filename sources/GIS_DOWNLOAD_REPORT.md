# 青藏高原青稞项目 — GIS 环境数据下载报告

- **下载完成时间**：2026-09-11
- **存放根目录**：`/ds2/workspace/ai/highland-barley-omics-gis/gis/raw/`
- **研究区域**：青藏高原及周边，经度 60–110°E，纬度 20–45°N
- **总数据量**：约 **33 GB**（5 个来源、5 个子目录）
- **校验结果**：**OK 173 / WARN 0 / FAIL 1**（唯一失败项为 CGIAR 服务器上确实不存在的瓦片，详见第 6 节）

> 本次为纯数据获取任务，未做任何分析；未修改 `DATA_CATALOG.md`、`RESULTS_STATUS.md`、`TODO.md`、`paper/` 等共享文档，也未做任何 git 操作。

---

## 1. 目录结构总览

```
gis/raw/
├── worldclim/    21 GB   WorldClim 2.1 生物气候 19 变量 + 海拔（全球 30″）
├── chelsa/      4.5 GB   CHELSA v2.1 生物气候 19 变量（全球 30″）
├── soilgrids/   2.4 GB   SoilGrids 2.0 土壤属性 8 项 × 2 层（区域子集）
├── elevation/   5.0 GB   SRTM 4.1 90 m 区域瓦片 49 个（zip + 解压 tif）
├── gbif/         38 MB   大麦分布点（全球 + 研究区子集）
├── GIS_DOWNLOAD_REPORT.md   本报告
├── verify_report.json       机器可读校验结果
└── *.sh / *.py              下载与校验脚本（可复现）
```

各来源文件清单与校验明细见下。

---

## 2. WorldClim 2.1

**目录**：`gis/raw/worldclim/`　**大小**：21 GB　**文件数**：22

### 2.1 下载文件

| 文件 | 大小 (bytes) | 状态 |
|---|---:|---|
| `wc2.1_30s_bio.zip` | 10,401,038,857 | ✅ 完整，`unzip -t` 无错误（19 个成员） |
| `wc2.1_30s_elev.zip` | 338,188,707 | ✅ 完整，`unzip -t` 无错误（1 个成员） |

**实际可用 URL（重要，与任务书不同）**：

```
https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip
https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip
```

> ⚠️ 任务书给出的 `https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/...` **已失效**（主机不可达）。
> WorldClim 已迁移至 `geodata.ucdavis.edu`，且路径中新增 `/climate/` 段。

### 2.2 解压后栅格（全球 30 弧秒，EPSG:4326）

`wc2.1_30s_bio_1.tif` … `wc2.1_30s_bio_19.tif`（19 个生物气候变量）+ `wc2.1_30s_elev.tif`（海拔）。

19 个 bio 文件合计约 10.8 GB，海拔文件 342 MB。zip 包已保留，未删除。

### 2.3 许可证（**与任务书假设不符，请务必注意**）

WorldClim 官方 About 页面原文：

> "The data are freely available for **academic use and other non-commercial use**.
> **Redistribution or commercial use is not allowed without prior permission.**
> Using the data to create maps for publishing of academic research articles is allowed.
> … You are allowed (but not required) to publish these articles (and the maps they contain)
> under an open license such as CC-BY…"

**结论：WorldClim 并非 CC BY 4.0，而是「学术/非商业用途免费，未经许可不得再分发或商用」。**
本项目为学术研究，使用与论文配图均被明确允许；但**不建议把原始 WorldClim 栅格再分发给第三方或公开仓库**。

推荐引用：Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. *International Journal of Climatology* 37(12): 4302–4315.

---

## 3. CHELSA v2.1

**目录**：`gis/raw/chelsa/`　**大小**：4.83 GB　**文件数**：19（全部成功）

1981–2010 气候平年值，30 弧秒，全球范围，`43200 × 20880` 像素，EPSG:4326，uint16/float32。

**实际可用 URL（重要，与任务书不同）**：

```
https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio{N}_1981-2010_V.2.1.tif
```

> ⚠️ 任务书给出的 `https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/...` **返回 404**。
> 经查证，CHELSA 已迁移到同一对象存储上的新桶路径 **`/chelsav2/`**（旧 `/envicloud/chelsa/chelsa_V2/` 作废）。
> 另注：该桶的 S3 列表接口不公开 CHELSA 前缀，无法通过列举发现，只能直连已知路径。

### 文件清单

| 文件 | 大小 (bytes) |
|---|---:|
| CHELSA_bio1_1981-2010_V.2.1.tif | 115,490,323 |
| CHELSA_bio2_1981-2010_V.2.1.tif | 65,178,902 |
| CHELSA_bio3_1981-2010_V.2.1.tif | 355,536,492 |
| CHELSA_bio4_1981-2010_V.2.1.tif | 457,053,433 |
| CHELSA_bio5_1981-2010_V.2.1.tif | 126,368,334 |
| CHELSA_bio6_1981-2010_V.2.1.tif | 117,685,874 |
| CHELSA_bio7_1981-2010_V.2.1.tif | 126,959,379 |
| CHELSA_bio8_1981-2010_V.2.1.tif | 129,024,254 |
| CHELSA_bio9_1981-2010_V.2.1.tif | 122,196,972 |
| CHELSA_bio10_1981-2010_V.2.1.tif | 123,145,240 |
| CHELSA_bio11_1981-2010_V.2.1.tif | 114,289,353 |
| CHELSA_bio12_1981-2010_V.2.1.tif | 655,151,886 |
| CHELSA_bio13_1981-2010_V.2.1.tif | 302,056,440 |
| CHELSA_bio14_1981-2010_V.2.1.tif | 178,915,498 |
| CHELSA_bio15_1981-2010_V.2.1.tif | 188,425,245 |
| CHELSA_bio16_1981-2010_V.2.1.tif | 473,777,041 |
| CHELSA_bio17_1981-2010_V.2.1.tif | 345,567,884 |
| CHELSA_bio18_1981-2010_V.2.1.tif | 438,957,395 |
| CHELSA_bio19_1981-2010_V.2.1.tif | 394,609,300 |

**monthly 数据（prec/temp）按任务书说明跳过**——bio 变量已满足第一版需求。如需补齐，URL 结构为同一桶下的 `.../climatologies/1981-2010/{prec,tas,...}/`。

### 许可证

CHELSA 官方数据集页面（chelsa-climate.org/datasets/chelsa_bioclim）明确标注 **"License: Creative Commons Zero - No Rights Reserved (CC0 1.0)"**。

注：部分第三方目录标注为 CC BY 4.0，与官方页面不一致；**以官方 CC0 1.0 为准**（CC0 更宽松，无合规风险）。

推荐引用：Karger, D.N. et al. (2017) *Scientific Data* 4:170122.

---

## 4. SRTM 4.1（90 m）

**目录**：`gis/raw/elevation/`　**大小**：5.0 GB（zip 1.77 GB + 解压 tif 3.53 GB）

### 4.1 CGIAR 服务器状态

任务书担心 CGIAR 已下线——**实测仍然可用**：

- 目录列表 `https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/` → **403**（禁目录索引，属正常）
- 单瓦片下载 → **HTTP 200 正常**

URL 模式：`https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_XX_YY.zip`

### 4.2 瓦片编号与覆盖范围（已实测验证）

编号规则：`XX` 自 180°W 起每 5° 递增（`XX = (lon+180)/5 + 1`），`YY` 自 60°N 起向南每 5° 递增（`YY = (60-lat)/5 + 1`）。

> ⚠️ 任务书给出的 "XX 为 12–22" 与本研究区不符。经用世界文件实测校验
> （`srtm_49_04.tfw` 原点 = 60.0004°E / 44.9996°N），研究区正确编号为：

- **X = 49 … 58**（对应 60–110°E）
- **Y = 4 … 8**（对应 45–20°N）
- 理论共 50 个瓦片

### 4.3 下载结果

| 项 | 数量 | 说明 |
|---|---:|---|
| 成功瓦片 | **49** | 全部 `unzip -t` 通过（每个含 4 个成员：tif/tfw/hdr/readme） |
| 缺失瓦片 | **1** | `srtm_49_08.zip`（60–65°E, 20–25°N）→ 服务器返回 **HTTP 404，该瓦片客观不存在** |

解压目录：`elevation/srtm_tif/`，49 个 `.tif` + 配套 `.tfw`/`.hdr`/`readme.txt`。

**关于缺失瓦片**：`srtm_49_08` 覆盖 60–65°E / 20–25°N（伊朗—阿曼一带），远在青藏高原核心区之外。
该处可由 **WorldClim `wc2.1_30s_elev.tif`**（已下载，全球覆盖）无缝补齐，不影响本研究。zip 与解压 tif 均已保留。

### 许可证

SRTM 为美国政府（NASA/NGA）成果，**属公有领域（public domain）**，CGIAR-CSI 仅做重投影与拼接分发，无使用限制。

---

## 5. SoilGrids 2.0

**目录**：`gis/raw/soilgrids/`　**大小**：2.49 GB　**文件数**：32（全部成功且校验通过）

### 5.1 获取方式

走 **ISRIC WCS 区域子集**（首选方案，成功），未使用 VRT/GDAL 备选路径。

```
https://maps.isric.org/mapserv?map=/map/{prop}.map&SERVICE=WCS&VERSION=2.0.1
  &REQUEST=GetCoverage&COVERAGEID={prop}_{depth}_mean
  &FORMAT=GEOTIFF_INT16
  &SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326
  &SUBSET=X(x0,x1)&SUBSET=Y(20,45)
```

**关键约束**：服务端 `MAXSIZE=16384`，而 60–110°E 在 250 m 分辨率下宽约 20000 像素，**一次性请求会报
`Raster size out of range`**。解决办法是把经度**拆成两段**：

- 段 A：`X(60,85)`　段 B：`X(85,110)`
- 每段输出 **11059 × 10461** 像素，int16，EPSG:4326，恰好覆盖目标范围（已实测边界精确等于 60/85/110°E 与 20/45°N）

### 5.2 属性与土层

下载了 **8 个属性 × 2 个土层 × 2 个经度段 = 32 个文件**：

| 属性 | 含义 | 0–5 cm（必需） | 5–15 cm（可选） |
|---|---|:--:|:--:|
| `clay` | 黏粒含量 | ✅ | ✅ |
| `sand` | 砂粒含量 | ✅ | ✅ |
| `silt` | 粉粒含量 | ✅ | ✅ |
| `soc` | 有机碳含量 | ✅ | ✅ |
| `phh2o` | pH（水） | ✅ | ✅ |
| `nitrogen` | 全氮 | ✅ | ✅ |
| `bdod` | 容重 | ✅ | ✅ |
| `cec` | 阳离子交换量 | ✅ | ✅ |

任务书要求「0–5 cm 必下、5–15 cm 时间够就下」——**两者均已全部完成**。

命名格式：`{prop}_{depth}_mean_X{x0}-{x1}_Y20-45.tif`，例如
`clay_0-5cm_mean_X60-85_Y20-45.tif`、`soc_5-15cm_mean_X85-110_Y20-45.tif`。

各文件大小在 **27 MB – 114 MB** 之间（phh2o 最小，soc/nitrogen 最大），总计 2,492,531,627 bytes。

### 5.3 校验与修复

初次下载后有 **5 个文件损坏**：ISRIC WCS 在传输中断时仍返回 HTTP 200，导致 curl 认为成功，
但 TIFF 数据实际被截断（其中两个文件截断字节数完全相同，均为 20,905,986，是断连的特征）。

已用「**严格校验 + 重下**」流程全部修复：校验方式为打开栅格并**解码最后 256 行像素条带**
（可捕获尾部瓦片缺失），修复后 **32/32 全部通过**。

损坏并已修复的 5 个文件：
`bdod_0-5cm_X60-85`、`phh2o_0-5cm_X85-110`、`sand_5-15cm_X85-110`、`soc_0-5cm_X60-85`、`soc_5-15cm_X85-110`。

### 许可证

ISRIC 官方（isric.org / docs.isric.org）明确标注：

> "The SoilGrids maps are publicly available on the SoilGrids website under the **CC-BY 4.0 License**."

即 **CC BY 4.0**，署名即可自由使用（含商用）。

推荐引用：Poggio, L. et al. (2021) SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty. *SOIL* 7:217–240.

---

## 6. GBIF 分布点

**目录**：`gis/raw/gbif/`　**大小**：38 MB　**文件数**：2

### 6.1 物种匹配

```
GET https://api.gbif.org/v1/species/match?name=Hordeum%20vulgare
→ usageKey / taxonKey = 2706056
  scientificName = "Hordeum vulgare L."   rank = SPECIES
  status = ACCEPTED   confidence = 99   matchType = EXACT
```

### 6.2 输出文件

| 文件 | 数据行数 | 说明 |
|---|---:|---|
| `gbif_hordeum_vulgare.tsv` | **93,536** | 全球所有带坐标记录，**未截断**（API 报告总数恰为 93,536，已全部下完） |
| `gbif_hordeum_vulgare_tibetan_region.tsv` | **9,302** | **额外补充**：仅研究区 60–110°E/20–45°N 内的记录 |

查询条件：`taxonKey=2706056&hasCoordinate=true&hasGeospatialIssue=false`，按 300 条/页翻页，
页间 sleep 0.5 s 礼貌访问。**0 个分页失败**，0 个非法坐标。制表符分隔，UTF-8。

字段（29 列）：
`key, scientificName, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, year, month, day, basisOfRecord, occurrenceStatus, countryCode, country, stateProvince, locality, elevation, elevationAccuracy, recordedBy, institutionCode, collectionCode, catalogNumber, datasetKey, publishingCountry, taxonKey, species, issues, establishmentMeans, occurrenceID, gbifID, license`

### 6.3 ⚠️ 重要发现：全球文件存在严重的地理取样偏差

对全球 93,536 条记录统计，**只有约 1.6% 落在研究区包围盒内**；记录数最多的国家是
**德国、法国、比利时、荷兰、瑞士、英国**——即欧洲。这是因为 GBIF 上 `Hordeum vulgare`
的观测记录绝大部分来自欧洲的农田监测与公民科学项目。

**研究区子集（9,302 条）的国家构成**则合理得多：
PK 2,617｜NP 1,983｜IN 1,607｜CN 1,528｜AF 620｜TJ 357｜UZ 279｜BT 129｜KG 101｜KZ 41｜TM 29｜IR 10

> **建议**：做青藏高原空间分析时**请使用 `gbif_hordeum_vulgare_tibetan_region.tsv`**，
> 直接用全球文件会因欧洲点占绝对多数而使任何环境 niche 分析被欧洲气候主导。
> 另注意中国境内仅 1,528 条，境内大麦分布点仍偏稀疏，必要时需补充国内种质资源或文献点位。

### 6.4 许可证（**混合，需逐条查看**）

GBIF 各数据集许可不同，**不是统一的一种**。本文件内 `license` 列的实际分布：

| 许可 | 全球文件条数 | 研究区子集条数 | 含义 |
|---|---:|---:|---|
| CC BY 4.0 | 74,355 | 7,048 | 署名即可（可商用） |
| **CC BY-NC 4.0** | **11,072** | **1,265** | **署名 + 禁止商用** |
| CC0 1.0 | 8,109 | 989 | 公有领域，无限制 |

> ⚠️ 任务书提示的 **CC BY-NC 4.0 确实存在**（约占 12%）。若本项目涉及商业化应用，
> 需按 `license` 列筛除 NC 记录（`grep -v "by-nc"`）。纯学术发表不受影响，但**建议在论文中注明**。

---

## 7. 校验方法与结果

校验脚本：`verify_downloads.py`（可重复运行），结果存 `verify_report.json`。

校验内容：
1. **体积比对**——本地文件大小 vs 服务器 `Content-Length`（WorldClim / CHELSA / SRTM 逐文件比对）
2. **压缩包完整性**——`unzip -t` 全量测试
3. **栅格完整性**——rasterio 打开并**实际解码像素**，同时核对 CRS、尺寸、边界范围
4. **坐标合法性**——GBIF 经纬度范围检查

### 最终结果

```
OK: 173   WARN: 0   FAIL: 1
```

**唯一失败项**：`elevation/srtm_49_08.zip` — CGIAR 服务器返回 HTTP 404，该瓦片客观不存在（见 4.3 节），
已用 WorldClim 海拔作为兜底，**不影响研究区使用**。

亮点：
- `wc2.1_30s_bio.zip`（10.4 GB，拼接重建）`unzip -t` **19 个成员全部无错误**
- 49 个 SRTM 瓦片全部通过 `unzip -t`
- 32 个 SoilGrids 栅格全部可完整解码，边界精确落在 60/85/110°E 与 20/45°N
- GBIF 93,536 条记录 0 非法坐标、0 分页失败

---

## 8. 过程中遇到的问题与解决办法

| # | 问题 | 解决办法 |
|---|---|---|
| 1 | WorldClim 任务书 URL 主机 `biogeo.ucdavis.edu` 已失效 | 改用 `geodata.ucdavis.edu/climate/worldclim/2_1/base/` |
| 2 | CHELSA 任务书 URL（`/envicloud/chelsa/chelsa_V2/`）404 | 改用迁移后的 `/chelsav2/` 路径 |
| 3 | CGIAR SRTM 单连接被限速至约 **22 KB/s**（50 瓦片需 16 小时以上） | 改为 **6 并发**抓取，实测升至 **4.8 MB/s（约 220×）** |
| 4 | SoilGrids WCS 报 `Raster size out of range, MAXSIZE=16384` | 经度**拆两段**（60–85 / 85–110）分别请求 |
| 5 | WorldClim 10.4 GB 单文件同样被限速（约 1.3 MB/s） | **HTTP Range 分 6 段并发**下载 + 本地按序重组，最终字节数与服务器完全一致 |
| 6 | 部分 Range 连接中途变慢（服务端随机限速） | 把剩余区间**切小（20 MB/块）+ 失败自动换连接重试**，稳定跑完 |
| 7 | 组装脚本漏拼 part1 的子块（一次性 bug） | 定位后 `truncate` 回滚并按正确顺序重组，`unzip -t` 验证无损坏 |
| 8 | **5 个 SoilGrids 文件被静默截断**（WCS 断连仍返回 HTTP 200） | 改用**严格校验**（curl 退出码 + 尾部像素解码）并重下，全部修复 |
| 9 | 系统 `awk` 为 mawk，不支持 `IGNORECASE`，导致 Content-Length 解析为空 | 改用 `grep -i` 解析 |
| 10 | GBIF 顺序翻页较慢 | 改为 6 并发分页抓取后按 offset 顺序合并 |

### 环境备注

- 系统无 GDAL（`gdal_translate` / `gdalinfo` 均不存在），已按预案 `pip install rasterio`（1.5.1）完成栅格校验
- 所有下载经由本机 socks/http 代理（`127.0.0.1:12334`），带宽与其他项目下载任务共享
- 未使用需要认证的 NASA Earthdata（按任务书要求已回避）

---

## 9. 复现脚本

`gis/raw/` 下保留了本次全部可复现脚本：

| 脚本 | 用途 |
|---|---|
| `download_worldclim.sh` / `dl_worldclim_par.sh` / `finish_worldclim2.sh` / `assemble_worldclim.sh` | WorldClim 下载、分段并发、断点续传与重组 |
| `dl_chelsa_par.sh` | CHELSA 19 个 bio 变量并发下载 |
| `download_srtm.sh` / `download_srtm_parallel.sh` | SRTM 瓦片下载（后者为并发版，推荐） |
| `download_soilgrids.sh` / `dl_soilgrids_par.sh` | SoilGrids WCS 区域子集（后者为并发版） |
| `repair_soilgrids_par.py` | SoilGrids 截断文件检测与修复 |
| `download_gbif.py` / `download_gbif_par.py` | GBIF 分布点下载（后者为并发版，推荐） |
| `extract_archives.sh` | 解压 WorldClim 与 SRTM 压缩包 |
| `verify_downloads.py` | 全量校验，输出 `verify_report.json` |

---

## 10. 总结与后续建议

### 完成度

| 来源 | 计划 | 实际 | 状态 |
|---|---|---|---|
| WorldClim 2.1 bio (30″) | 19 变量 | 19 变量，10.4 GB | ✅ 完成 |
| WorldClim 2.1 elev | 1 层 | 1 层，338 MB | ✅ 完成 |
| CHELSA v2.1 bio1–19 | 19 变量 | 19 变量，4.83 GB | ✅ 完成 |
| CHELSA monthly | 可选 | — | ⏭️ 按说明跳过 |
| SRTM 4.1 90 m | 50 瓦片 | **49 瓦片** | ⚠️ 1 片客观不存在，已用 WorldClim elev 兜底 |
| SoilGrids 0–5 cm | 8 属性 | 8 属性 × 2 段 | ✅ 完成 |
| SoilGrids 5–15 cm | 可选 | 8 属性 × 2 段 | ✅ 完成（超出要求） |
| GBIF 分布点 | ≤10 万条 | 93,536 条（全量）+ 研究区 9,302 条 | ✅ 完成，无截断 |

### 后续建议

1. **裁剪到研究区**：WorldClim / CHELSA / SRTM 目前为全球或大区域范围，
   进入分析前建议统一裁剪到 60–110°E / 20–45°N（SoilGrids 已是区域子集）。
2. **统一坐标系与分辨率**：WorldClim 与 CHELSA 均为 EPSG:4326 30″（约 1 km），
   SRTM 为 90 m，SoilGrids 为 250 m。建模前需重采样到统一网格。
3. **GBIF 用研究区子集**：务必使用 `gbif_hordeum_vulgare_tibetan_region.tsv`，
   避免全球文件的欧洲偏差（见 6.3）。
4. **注意 NC 许可**：若涉及商业化，需按 `license` 列剔除 CC BY-NC 记录。
5. **WorldClim 不要再分发**：仅学术/非商业可用，且不得未经许可再分发原始栅格。
6. **CHELSA monthly** 如后续需要（如计算生长季指标），可在同一 `/chelsav2/` 桶下按
   `.../climatologies/1981-2010/{prec,tas,tasmax,...}/` 路径补充。

---

*报告生成于 2026-09-11；数据经逐文件体积比对、压缩包完整性测试与栅格像素级解码校验。*
