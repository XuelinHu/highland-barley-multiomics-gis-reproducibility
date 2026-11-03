#!/usr/bin/env python3
"""Correlation screening + dataset build report.

Outputs
-------
processed/correlation_matrix_full.tsv      (all env variables, cleaned points)
processed/correlation_matrix_plateau.tsv   (all env variables, qingke candidates)
processed/correlation_high.tsv             (pairs with |r| >= 0.9)
processed/DATASET_BUILD_REPORT.md          (Chinese report)
"""
import json
import os
import time

import numpy as np
import pandas as pd

PROC = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed')
RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
RTHRESH = 0.9

TERRAIN = ["elev", "slope", "aspect"]
WC_VARS = [f"wc_bio{i}" for i in range(1, 20)]
CH_VARS = [f"ch_bio{i}" for i in range(1, 20)]
SOIL_PROPS = ["clay", "sand", "silt", "soc", "phh2o", "nitrogen", "bdod", "cec"]
SOIL_VARS = [f"sg_{p}_{d}" for p in SOIL_PROPS for d in ("0_5", "5_15")]
ENV = TERRAIN + WC_VARS + CH_VARS + SOIL_VARS


def corr_pairs(df, cols, min_n=30):
    """Pairwise-complete Pearson correlations; returns long DataFrame."""
    sub = df[cols]
    r = sub.corr(method="pearson", min_periods=min_n)
    n = sub.notna().astype(np.int32).T.dot(sub.notna().astype(np.int32))
    rows = []
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            rows.append((a, b, r.loc[a, b], int(n.loc[a, b])))
    return pd.DataFrame(rows, columns=["var1", "var2", "r", "n"])


def main():
    t0 = time.time()
    full = pd.read_csv(os.path.join(PROC, "barley_env_matrix.tsv.gz"),
                       sep="\t", low_memory=False)
    plat = pd.read_csv(os.path.join(PROC, "qingke_plateau_subset.tsv.gz"),
                       sep="\t", low_memory=False)
    print(f"loaded full={len(full)} plateau={len(plat)}", flush=True)

    pf = corr_pairs(full, ENV)
    pp = corr_pairs(plat, ENV)
    pf.to_csv(os.path.join(PROC, "correlation_matrix_full.tsv"), sep="\t", index=False)
    pp.to_csv(os.path.join(PROC, "correlation_matrix_plateau.tsv"), sep="\t", index=False)

    for df in (pf, pp):
        df["abs_r"] = df.r.abs()
    hi = pf.merge(pp, on=["var1", "var2"], suffixes=("_full", "_plateau"))
    mask = (hi.abs_r_full >= RTHRESH) | (hi.abs_r_plateau >= RTHRESH)
    hi = hi[mask].copy()
    hi["high_in"] = np.where(
        (hi.abs_r_full >= RTHRESH) & (hi.abs_r_plateau >= RTHRESH), "both",
        np.where(hi.abs_r_full >= RTHRESH, "all_points", "plateau_only"))
    hi["note"] = np.where(
        hi.var1.eq("aspect") | hi.var2.eq("aspect"),
        "aspect is circular (0-360 deg) and undefined on flat ground; Pearson r is not interpretable",
        "")
    hi = hi.sort_values("abs_r_full", ascending=False)
    hi.to_csv(os.path.join(PROC, "correlation_high.tsv"), sep="\t", index=False)
    print(f"high pairs |r|>={RTHRESH}: {len(hi)}", flush=True)

    # a square matrix for convenience / downstream plotting
    full[ENV].corr(min_periods=30).to_csv(
        os.path.join(PROC, "correlation_matrix_full_square.tsv"), sep="\t")

    # ---------------- report -------------------------------------------
    clean = json.load(open(os.path.join(PROC, "cleaning_counts.json")))
    covj = json.load(open(os.path.join(PROC, "env_coverage.json")))
    recon = json.load(open(os.path.join(PROC, "recon_rasters.json")))

    def tsv_table(df, cols=None, floatfmt="%.4f"):
        cols = cols or list(df.columns)
        out = ["| " + " | ".join(cols) + " |",
               "|" + "|".join(["---"] * len(cols)) + "|"]
        for _, r in df.iterrows():
            cells = []
            for c in cols:
                v = r[c]
                if isinstance(v, float):
                    cells.append("" if pd.isna(v) else floatfmt % v)
                else:
                    cells.append(str(v))
            out.append("| " + " | ".join(cells) + " |")
        return "\n".join(out)

    L = []
    A = L.append
    A("# 青稞/大麦环境数据集构建报告 (DATASET_BUILD_REPORT)\n")
    A(f"生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}  ")
    A(f"脚本: `gis_build/build_rasters.py`, `clean_gbif.py`, `extract_env.py`, `finalize_dataset.py`  ")
    A(f"输出目录: `{PROC}`\n")

    A("## 0. 产物清单\n")
    A("| 文件 | 内容 |")
    A("|---|---|")
    A(f"| `barley_env_matrix.tsv.gz` | 清洗后全部样点 ({len(full)} 行 x {full.shape[1]} 列) x 环境变量 + 元数据 + 标记列 |")
    A(f"| `qingke_plateau_subset.tsv.gz` | `is_qingke_candidate == True` 子集 ({len(plat)} 行) |")
    A(f"| `correlation_high.tsv` | \\|r\\| >= {RTHRESH} 的变量对 ({len(hi)} 对) |")
    A("| `correlation_matrix_full.tsv` / `correlation_matrix_plateau.tsv` | 全部变量两两相关(长表) |")
    A("| `correlation_matrix_full_square.tsv` | 相关矩阵(方阵, 便于绘图) |")
    A("| `variable_coverage.tsv` / `env_coverage.json` | 逐变量覆盖率与极值 |")
    A("| `cleaning_counts.json` | GBIF 清洗各步数量 |")
    A("| `recon_rasters.json` | 输入栅格 CRS/分辨率/范围盘点 |")
    A("| `rasters/srtm/srtm_mosaic.tif`, `srtm_slope.tif`, `srtm_aspect.tif` | SRTM 拼接 DEM (60-110E/20-45N) + 坡度/坡向 |")
    A("| `rasters/soil/sg_*.tif` | SoilGrids 8 属性 x 2 深度 区域拼接栅格 |\n")

    A("## 1. GBIF 分布点清洗\n")
    A("逐步数量 (每步在上一步结果上继续过滤):\n")
    st = pd.DataFrame(clean["steps"])[["step", "before", "after", "removed", "note"]]
    A(tsv_table(st, floatfmt="%s").replace("nan", ""))
    A("")
    A(f"- 原始记录 **{clean['raw_records']}** 条 -> 清洗后 **{clean['clean_records']}** 条 "
      f"(保留率 {100*clean['clean_records']/clean['raw_records']:.1f}%)")
    A(f"- 去重规则: 坐标四舍五入到 4 位小数后去重; 同一坐标保留 "
      f"`coordinateUncertaintyInMeters` 最小(最精确)、其次年份最新的记录")
    A(f"- 年份范围: {clean['year_range'][0]}-{clean['year_range'][1]}, 中位数 {clean['year_median']}")
    A(f"- basisOfRecord 保留集: {', '.join(clean['allowed_basis'])}")
    A(f"- 被剔除的记录类型: " +
      ", ".join(f"{k}={v}" for k, v in clean["basis_of_record_dropped"].items()))
    A(f"- occurrenceStatus 分布: {clean['occurrence_status']} (ABSENT 记录**未**剔除, 已保留该列供下游过滤)")
    A(f"- establishmentMeans 分布: {clean['establishment_means']} "
      f"(introduced 记录已保留, 建议栽培/逸生分析时过滤)")
    A(f"- 清理后点位空间范围: lon [{clean['lon_min']:.3f}, {clean['lon_max']:.3f}], "
      f"lat [{clean['lat_min']:.3f}, {clean['lat_max']:.3f}]")
    A(f"- 落在青藏高原矩形框 (70-105E, 25-40N) 内的点: {clean['n_on_plateau_box']}\n")
    A("清洗后记录数最多的国家:\n")
    tc = pd.DataFrame(list(clean["top_countries"].items()), columns=["country", "n"])
    A(tsv_table(tc, floatfmt="%s").replace("nan", ""))
    A("")
    A("鉴定名 (scientificName) 分布 (前 10):\n")
    sn = pd.DataFrame(list(clean["scientificName_top"].items())[:10],
                      columns=["scientificName", "n"])
    A(tsv_table(sn, floatfmt="%s").replace("nan", ""))
    A("")

    A("## 2. 栅格输入盘点 (CRS / 分辨率)\n")
    A("所有输入栅格经检查**均为 EPSG:4326** (WGS84 经纬度), 无需重投影。\n")
    A("| 数据源 | 文件数 | 分辨率 | 像素尺寸 | 覆盖范围 | 类型/NoData |")
    A("|---|---|---|---|---|---|")
    srtm_files = recon["srtm"]
    A(f"| SRTM 90m (CGIAR) | {len(srtm_files)} 瓦片 | 3 arcsec (0.0008333 deg, ~90 m) | "
      f"{srtm_files[0]['width']}x{srtm_files[0]['height']} /瓦片 | 60-110E, 20-45N | int16, -32768 |")
    wc = recon["worldclim_bio"]
    A(f"| WorldClim 2.1 bio | {len(wc)} | 30 arcsec (0.0083333 deg) | "
      f"{wc[0]['width']}x{wc[0]['height']} | 全球 -180..180, -90..90 | float32, -3.4e38 |")
    A(f"| WorldClim 2.1 elev | 1 | 30 arcsec | 43200x21600 | 全球 | int16, -32768 |")
    ch = recon["chelsa"]
    A(f"| CHELSA V2.1 bio | {len(ch)} | 30 arcsec | {ch[0]['width']}x{ch[0]['height']} | "
      f"全球 -180..180, -90..84 | uint16 (bio3 为 float32), 无 NoData 标记 |")
    sg = recon["soilgrids"]
    A(f"| SoilGrids 250m | {len(sg)} (8 属性 x 2 深度 x 2 经度块) | "
      f"{sg[0]['res'][0]:.6f} x {sg[0]['res'][1]:.6f} deg | {sg[0]['width']}x{sg[0]['height']} /块 | "
      f"60-110E, 20-45N | int16, 水体/无数据 = 0 |")
    A("")
    A("处理后的派生栅格:\n")
    A("| 文件 | 尺寸 | 分辨率 | 说明 |")
    A("|---|---|---|---|")
    A("| `rasters/srtm/srtm_mosaic.tif` | 60000x30000 | 3 arcsec | 49 个 SRTM 瓦片拼接 (60-110E, 20-45N), int16, nodata -32768, LZW |")
    A("| `rasters/srtm/srtm_slope.tif` | 60000x30000 | 3 arcsec | 坡度 (度), float32, Horn 3x3 窗口 |")
    A("| `rasters/srtm/srtm_aspect.tif` | 60000x30000 | 3 arcsec | 坡向 (度, 正北顺时针), float32, 坡度<0.1 度处为 NaN |")
    A("| `rasters/soil/sg_{prop}_{depth}cm_mean.tif` | 22118x10461 | ~250 m | 8 属性 x 2 深度, 两个经度块拼接, int16 |")
    A("")

    A("## 3. 变量提取与覆盖率\n")
    A(f"共 {len(ENV)} 个环境变量 (地形 3 + WorldClim 19 + CHELSA 19 + 土壤 16), "
      f"在 {len(full)} 个样点上采样。\n")
    A("### 3.1 地形\n")
    A(f"- `elev`: 先取 SRTM, 缺失时用 WorldClim elev 兜底。海拔来源统计: "
      f"{covj['elev_source_counts']}")
    A(f"- 最终 `elev` 非空 {covj['coverage']['elev']['n_valid']} / {covj['n_points']} "
      f"({covj['coverage']['elev']['pct_valid']}%)")
    A("- 坡度/坡向用 3x3 Horn 窗口计算; 经向像元间距按纬度换算为米 "
      "(`dx = 0.0008333 deg * 111194.9 * cos(lat)`), 纬向 `dy = 0.0008333 deg * 111194.9 m`")
    A("- 坡向在坡度 < 0.1 度处置为 NaN (平坦地面坡向无定义)\n")
    A("### 3.2 逐变量覆盖率\n")
    cov = pd.DataFrame(covj["coverage"]).T.reset_index().rename(columns={"index": "variable"})
    cov = cov[["variable", "unit", "n_valid", "pct_valid", "min", "median", "max", "source"]]
    A(tsv_table(cov, floatfmt="%.4g"))
    A("")
    A(f"- WorldClim 陆地像元 (wc_bio1 非空) 的点: {covj['n_wc_land']} / {covj['n_points']} "
      f"({100*covj['n_wc_land']/covj['n_points']:.1f}%); 其余为海洋/WorldClim 掩膜区域")
    A(f"- 高原矩形框内点: {covj['n_on_plateau']}; `is_qingke_candidate` "
      f"(框内 且 elev>=2000m): **{covj['n_qingke_candidates']}**")
    A("- 土壤变量仅覆盖 60-110E / 20-45N, 区域外为 NaN; 区域内水体像元 (SoilGrids 编码为 0) 亦置为 NaN")
    A("- CHELSA 转换公式 (已用 WorldClim 同点对照验证): "
      "bio1/5/6/8/9/10/11 = raw/10 - 273.15 (K->degC), bio2/bio7 = raw/10 (degC), "
      "bio3 = raw/10*100 (%), bio4/bio15 = raw/10, bio12-19 = raw/10 (mm)")
    A("- SoilGrids 换算系数: clay/sand/silt, soc, phh2o, cec = raw/10; nitrogen, bdod = raw/100\n")

    A(f"## 4. 相关性检查 (Pearson, |r| >= {RTHRESH})\n")
    A(f"在全部 {len(full)} 个清洗点上计算 {len(ENV)} 个变量的两两 Pearson 相关; "
      f"同时在高原候选子集 ({len(plat)} 点) 上重算。"
      f"共 **{len(hi)}** 对变量在至少一个集合中 |r| >= {RTHRESH}。\n")
    A("**注意: 变量未做任何删除, 仅在此标注。**\n")
    A("相关为成对完整观测 (pairwise complete) 计算, `n_full` / `n_plateau` 两列给出该变量对"
      "实际使用的样本量。土壤变量只在 60-110E/20-45N 区域内有值, 因此涉及土壤的变量对"
      "样本量远小于其他对, 其相关系数反映的是高原区域内的关系, 与其他对不可直接比较。\n")
    A("这 55 对绝大多数是**同一物理量的不同表达**, 可分三类: "
      "(a) 同一变量的两个土壤深度 (`sg_*_0_5` vs `sg_*_5_15`); "
      "(b) 同一变量的 WorldClim 与 CHELSA 两个来源 (`wc_bioN` vs `ch_bioN`); "
      "(c) 同一来源内季节含义重叠的生物气候变量 (如 bio13 最湿月 vs bio16 最湿季, bio5 最热月 vs bio10 最热季)。"
      "真正的跨类相关只有年降水与土壤 pH 两对 (`wc_bio12` vs `sg_phh2o_*`, r≈-0.91, "
      "即降水越多土壤越酸), 属于可解释的真实关系。"
      "**建模时应在每个共线簇内只保留一个变量** (如气候在 WorldClim/CHELSA 中二选一)。\n")
    show = hi[["var1", "var2", "r_full", "n_full", "r_plateau", "n_plateau", "high_in"]].copy()
    A(tsv_table(show.head(80), floatfmt="%.4f"))
    if len(hi) > 80:
        A(f"\n(仅显示前 80 对, 完整 {len(hi)} 对见 `correlation_high.tsv`)\n")
    A("")

    A("## 5. 已知局限\n")
    A("1. **无法区分青稞与其他大麦**: GBIF 记录仅有 `Hordeum vulgare` 种级鉴定, "
      "`scientificName` 中的变种名 (如 var. `distichon`/`nudum`) 不完整且不可靠。"
      "`is_qingke_candidate` 只是「高原框 + 海拔>=2000m」的地理代理, "
      "**不代表该记录一定是青稞**; 里面包含大量栽培大麦、野生/逸生个体与错误鉴定。")
    A("2. **栽培/逸生记录混杂**: 清洗后仍有 "
      f"{clean['establishment_means'].get('introduced', 0)} 条 `establishmentMeans=introduced` "
      "记录; 大麦全球广泛栽培, 欧洲/美洲/澳洲的大量 HUMAN_OBSERVATION 多为栽培或逸生个体, "
      "做生态位建模时建议结合 `establishmentMeans`、`occurrenceStatus`、`basisOfRecord` 进一步筛选。")
    A("3. **土壤变量仅高原区域**: SoilGrids 只下载了 60-110E / 20-45N, "
      "区域外样点的 16 个土壤变量全为 NaN, 区域内的水体像元也被置 NaN; "
      "因此涉及土壤变量的多元分析只能在区域内子集上做 (有效样本量见覆盖率表)。")
    A("4. **SRTM 拼接范围有限**: SRTM 只覆盖 60-110E / 20-45N (青藏高原及周边), "
      f"46102 个点中只有 {covj['coverage']['elev_srtm_m']['n_valid']} 个 "
      f"({covj['coverage']['elev_srtm_m']['pct_valid']}%) 落在该范围内并取到 SRTM 海拔, "
      f"其余 {covj['elev_source_counts'].get('WorldClim', 0)} 个点由 WorldClim 30s DEM 兜底 "
      "(`elev_source` 列标出); `slope`/`aspect` 因而也只有这 "
      f"{covj['coverage']['slope']['n_valid']} 个点有值 "
      f"(`aspect` 为 {covj['coverage']['aspect']['n_valid']} 个, 差额是坡度<0.1度处坡向无定义的点), "
      "且分辨率 90m 与气候变量 1km 不匹配 (提取时为点采样, 未做尺度聚合)。"
      "SRTM 拼接栅格本身约 9.8% 为 nodata (主要是区域南部的海洋)。")
    A("5. **气候背景年代不匹配**: WorldClim 2.1 bio 为 1970-2000 常年值, CHELSA V2.1 为 "
      "1981-2010 常年值, 而清洗后记录年份跨 1950-2026 (中位数 "
      f"{clean['year_median']}), 大部分记录晚于气候基准期。"
      "本数据集把 1950 年以后的点都按同一 1981-2010 气候背景赋值, "
      "未考虑气候变化趋势, 也**未**为每个点提取记录当年气候。")
    A("6. **CHELSA 海洋像元非掩膜值**: CHELSA V2.1 的 bio 栅格没有 NoData 标记, "
      "海洋/冰川上也有数值 (WorldClim 在同一位置为 NoData)。"
      "本数据集用 `wc_land` 列 (WorldClim bio1 是否非空) 标记陆地/海洋, "
      "需要严格陆地子集时请用 `wc_land == True` 过滤。")
    A("7. **坐标不确定性未过滤**: 清洗时未按 `coordinateUncertaintyInMeters` 设阈值 "
      f"(该列有 {int(full.coord_uncertainty_m.isna().sum())} 条缺失), "
      "只做了 4 位小数去重; 部分人类观测点坐标精度可达数公里。")
    A("8. **坡向是环形变量**: `aspect` (0-360 度) 在平坦处无定义, "
      "Pearson 相关系数对其不可解释, 相关性表中的 aspect 相关对请忽略或改用 sin/cos 变换。")
    A("9. **CHELSA/WorldClim 精度差异**: 两者在高海拔地区 (青藏高原) 差异明显 "
      "(如拉萨年降水 CHELSA 519mm vs WorldClim 418mm), "
      "相关性分析中二者同类变量高度共线 (见 `correlation_high.tsv`), 建模时应二选一。")
    A("10. **GBIF 采样偏差**: 记录空间分布高度不均 (欧洲/北美/澳洲城市与标本馆周边密集), "
      "且 1950 年后记录以人类观测为主, 存在明显的采样强度偏差, 未做空间稀疏化处理。\n")

    A("## 6. 复现\n")
    A("```bash")
    A("cd /ds1/workspace/ai/highland-barley-omics-gis")
    A("python3 gis_build/build_rasters.py     # SRTM 拼接 + 坡度坡向 + SoilGrids 拼接")
    A("python3 gis_build/clean_gbif.py        # GBIF 清洗 -> points_clean.tsv.gz")
    A("python3 gis_build/extract_env.py       # 环境变量提取 -> barley_env_matrix.tsv.gz")
    A("python3 gis_build/finalize_dataset.py  # 相关性 + 本报告")
    A("```")

    with open(os.path.join(PROC, "DATASET_BUILD_REPORT.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"report written, {len(hi)} high-corr pairs, {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
