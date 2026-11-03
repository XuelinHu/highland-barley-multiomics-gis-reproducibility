#!/usr/bin/env python3
"""Generate STRUCTURE_REPORT.md from the analysis artefacts (Chinese)."""
import os, json, csv
os.chdir(os.path.dirname(os.path.abspath(__file__)))

S = {r['gene_id']: r for r in json.load(open('structure_plddt_summary.json'))}
CF = json.load(open('colabfold_scores.json'))
DC = {r['gene_id']: r for r in json.load(open('domain_consistency.json'))}
AF = json.load(open('afdb/afdb_meta.json'))
MSA = json.load(open('msa_depth.json'))
MOB = json.load(open('disorder_annotation_mobidb.json'))
ANN = json.load(open('domain_annotation.json'))
LAB = json.load(open('figure_labels.json'))
targets = json.load(open('targets.json'))
DIR = {r['gene_id']: r for r in csv.DictReader(
    open(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/annotation/top100_annotation_full.tsv')),
    delimiter='\t')}

ORDER = list(LAB.keys())
DIS = ['HORVU.MOREX.r3.5HG0507410', 'HORVU.MOREX.r3.1HG0080020',
       'HORVU.MOREX.r3.7HG0637810', 'HORVU.MOREX.r3.2HG0197000']
ORD = [g for g in ORDER if g not in DIS]
NSHORT = {g: LAB[g].split(' (')[0] for g in ORDER}

def f(x, n=2):
    if x is None or x == '':
        return '—'
    return f'{x:.{n}f}' if isinstance(x, float) else str(x)

def srcname(g):
    s = S[g]['source']
    if 'exact' in s:
        return f"AFDB {s.split(':')[1].split()[0]}（序列完全一致）"
    if s.startswith('AFDB'):
        return f"AFDB 同源 {s.split(':')[1]}"
    return 'ColabFold 本地'

def modelname(g):
    """AFDB 官方条目没有 ColabFold 模型名；本地 ColabFold 才展示 rank_001 模型名。"""
    if 'exact' in S[g]['source']:
        return 'AFDB 官方模型'
    return CF.get(g, {}).get('model', '—')

def ags(x, n=2):
    """数值 + Å 单位；缺失值不留悬空单位。"""
    v = f(x, n)
    return f'{v} Å' if v != '—' else '—'

# "有结构域"的对照靶标。NAC 单列：它 N 端有序、C 端是家族固有的无序激活区，
# 整体低置信占比与无序组重叠，混在对照组里会污染对照（见 3.1 / 5.4）。
ORD_STRUCT = [g for g in ORD if g != 'HORVU.MOREX.r3.5HG0524540']
NAC = 'HORVU.MOREX.r3.5HG0524540'

L = []
w = L.append
w('# top100 干旱候选关键蛋白结构预测报告')
w('')
w('> TODO 阶段 4 的 AlphaFold 项。全部计算在 '
  '`/ds2/workspace/ai/highland-barley-omics-gis/processed/reanalysis/structures/` 内完成；')
w('> 未改动 `processed/reanalysis/` 下的任何既有文件，未执行任何 git 操作。')
w('')
w('---')
w('')
w('## 1. 方法：降级路径实际走到哪一步')
w('')
w('| 优先级 | 方法 | 实际使用情况 | 靶标数 |')
w('|---|---|---|---|')
w('| 1 | AlphaFold DB 直接下载（PDB + pLDDT） | 9 个同源条目全部下载；其中 2 个与靶标序列**逐字符相同**，直接作为该靶标结构采用 | 2 采用 / 9 下载 |')
w('| 2 | **ColabFold Local 1.6.2 本地 GPU 预测** | 成功。AF2 权重，5 模型 × 3 recycle，`--model-type auto`，MSA 走远程服务器 | **8** |')
w('| 3 | ESMFold（后备） | **未使用**——ColabFold 安装与 GPU 推理均一次成功，无需降级 | 0 |')
w('| 4 | 全部失败 | **未发生** | 0 |')
w('')
w('**运行环境**：RTX 3090 24 GB，驱动 580.126.09 / CUDA 13.0，无 sudo。')
w('用 `uv` 建独立 venv（Python 3.11）装 `colabfold[alphafold]`；')
w('默认解析到的 `jaxlib 0.10.2` 是 CPU 版，补装 `jax[cuda12]==0.10.2` 后 `jax.devices()` 返回 `CudaDevice(id=0)`，')
w('CUDA 12 运行时在 CUDA 13 驱动上通过向下兼容正常工作（安装约 16 分钟、约 5 GB nvidia wheel）。')
w('全流程 10 个靶标端到端耗时约 30 分钟（含 MSA 检索），峰值显存约 3.2 GB。')
w('')
w('### 1.1 `uniprot_ac` 列的一个关键说明（影响了靶标选择的前提）')
w('')
w('任务给的流程是"用 `top100_annotation_full.tsv` 的 `uniprot_ac` 列查 AlphaFold DB，查到就不再本地预测"。')
w('核对后确认：**该列的 AC 是 BLAST 比对到 Swiss-Prot 的最佳同源蛋白，不是大麦靶标自身的 UniProt 条目**。')
w('例如 `HORVU.MOREX.r3.5HG0524540` 的 AC 是 `Q75HE5` = 水稻 NAC domain-containing protein 2；')
w('`HORVU.MOREX.r3.2HG0099870` 是 `Q08277` = 玉米 HSP82。')
w('只有 DHN2(`P12951`) 与 HVA1(`P14928`) 本身就是大麦蛋白，序列与靶标完全一致。')
w('')
w('如果照字面执行，8 个靶标会拿到"水稻/玉米/拟南芥同源蛋白"的结构并当作大麦靶标结论使用——'
  '这是来源误述。因此实际处理为：')
w('')
w('- 序列完全一致的 2 个（DHN2、HVA1）→ 按原意直接用 AFDB，不做本地预测；')
w('- 其余 8 个 → 一律用 ColabFold 预测**大麦自身序列**；')
w('- 9 个同源 AFDB 结构仍全部下载保存在 `afdb/`，作为"该家族应该长什么样"的参照系，'
  '并在第 5 节用作结构一致性检查的比对对象。')
w('')
w('---')
w('')
w('## 2. 每个靶标的 pLDDT 汇总')
w('')
w('pLDDT 逐残基从 PDB 的 B-factor 列读取（AFDB 与 ColabFold 输出都用该列存 pLDDT）。')
w('读取口径已校验：对 `AF-P12951` 复算得均值 54.64，AFDB API 上报 `globalMetricValue` = 54.66。')
w('')
w('| 基因 ID | 注释 | 方向 | 长度 | 结构来源 | 最佳模型 | pLDDT 均值 | 中位数 | pTM | <50 | 50–70 | 70–90 | >90 |')
w('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for g in ORDER:
    r, c = S[g], CF.get(g, {})
    w(f"| {g} | {NSHORT[g]} | {DIR.get(g,{}).get('direction','')} | {r['seq_len']} | {srcname(g)} | "
      f"{modelname(g)} | **{f(r['mean'],1)}** | {f(r['median'],1)} | "
      f"{f(r.get('ptm'),3)} | {f(r['fraction_very_low_lt50'],3)} | {f(r['fraction_low_50_70'],3)} | "
      f"{f(r['fraction_confident_70_90'],3)} | {f(r['fraction_very_high_ge90'],3)} |")
w('')
w('MSA 深度（ColabFold 远程服务器返回）：'
  + '；'.join(f'{NSHORT[g]} {MSA[g]}' for g in ORDER) + ' 条。')
w('')
w('说明：DHN2 与 HVA1 采用的是 AFDB 官方模型；这两个序列本流程也各跑了一遍 ColabFold 作交叉验证，'
  '结果见 3.2（模型名与 pLDDT 见 `colabfold_scores.json`）。'
  'LEA_1 candidate B 的 MSA 只有 159 条，是全部靶标中最少的，其低置信有 MSA 偏浅的贡献——'
  '但它同时满足"<70 占比 100%""有序核心 0 残基""MobiDB-lite 全程判无序"三条，结论仍成立。')
w('')
w('图：`structure_plddt_summary.png`（逐残基 pLDDT 曲线，2×5）；'
  '`structure_gallery.png`（3D cartoon 按 pLDDT 着色，3Dmol.js + headless Chrome 离线渲染）。')
w('')
w('### 2.1 分档')
w('')
for lbl, lo, hi in [('高置信（均值 ≥ 85，结构可直接用于分析）', 85, 999),
                    ('中高（均值 70–85，折叠可信、部分 loop 需谨慎）', 70, 85),
                    ('低与无序（均值 < 70，多数残基无固定折叠）', -1, 70)]:
    gs = [g for g in ORDER if lo <= (S[g]['mean'] or 0) < hi]
    if gs:
        w(f"- **{lbl}**：" + '、'.join(f'{NSHORT[g]} {f(S[g]["mean"],1)}' for g in gs))
w('')
w('---')
w('')
w('## 3. 无序蛋白（DHN / LEA / HVA1）的无序性证据')
w('')
w('三条独立证据链指向同一结论：这几个蛋白在溶液中没有固定三维折叠。')
w('')
w('### 3.1 证据一：AlphaFold2 pLDDT 低置信区段占比')
w('')
w('| 基因 | 蛋白 | pLDDT 均值 | <50 占比 | 50–70 占比 | **<70 占比（无固定折叠）** | 最长连续 <50 区段 | pTM |')
w('|---|---|---|---|---|---|---|---|')
for g in DIS:
    r = S[g]; d = DC[g]
    w(f"| {g} | {NSHORT[g]} | {f(r['mean'],1)} | {f(r['fraction_very_low_lt50'],3)} | "
      f"{f(r['fraction_low_50_70'],3)} | **{f(r['fraction_very_low_lt50'] + r['fraction_low_50_70'],3)}** | "
      f"{d['longest_low_plddt_run']} aa | {f(r.get('ptm'),3)} |")
w('')
w('对照（`HORVU.MOREX.r3.4HG0386900` 之外的 5 个"有结构域"靶标）：')
w('<50 占比 ' + '、'.join(f'{NSHORT[g]} {f(S[g]["fraction_very_low_lt50"],3)}' for g in ORD_STRUCT) + '；')
w('**<70 占比 ' + '、'.join(f'{NSHORT[g]} {f(S[g]["fraction_very_low_lt50"]+S[g]["fraction_low_50_70"],3)}'
                            for g in ORD_STRUCT) + '**。')
w('')
w('即：这 5 个靶标 **82%–88% 的残基**都能达到 pLDDT≥70（折叠可信），')
w('而无序组 **77%–100% 的残基**落在 pLDDT<70（无固定折叠）。两组**完全不重叠**，分界清晰。')
w('')
w(f'（`NAC/NAM TF` 不算在内：它是唯一的中间情形，<70 占比 '
  f'{f(S[NAC]["fraction_very_low_lt50"]+S[NAC]["fraction_low_50_70"],3)}，与无序组重叠——'
  '因为它是"N 端有序结构域 + C 端长无序激活区"的嵌合架构（见 5.4），')
w('不属于纯无序蛋白，因此不计入无序组。）')
w('')
w('### 3.2 证据二：与 AlphaFold DB 官方预测的交叉验证')
w('')
w('DHN2 与 HVA1 既有 AFDB（官方 AF2 单体预测）又有本流程 ColabFold 独立跑的同一序列，')
w('两条完全独立的流程给出一致的结论：')
w('')
w('| 基因 | AFDB pLDDT 均值 | 本流程 ColabFold 均值 | 差值 | 两者 CA-RMSD（全链） |')
w('|---|---|---|---|---|')
for g in ['HORVU.MOREX.r3.5HG0507410', 'HORVU.MOREX.r3.1HG0080020']:
    d = DC[g]; cf = CF[g]['plddt_mean']
    w(f"| {g} | {f(S[g]['mean'],1)} | {f(cf,1)} | {f(abs(S[g]['mean']-cf),1)} | {f(d['ca_rmsd_A'],1)} Å |")
w('')
w('（DHN2 两条流程整体 RMSD 21.9 Å、HVA1 17.4 Å：折叠自由度大、几乎每个模型都不同，')
w('这正是无序蛋白的典型表现——**不是流程不稳定，而是本来就没有唯一构象**。）')
w('')
w('### 3.3 证据三：MobiDB-lite 独立无序预测')
w('')
w('上游 InterProScan 已带 MobiDB-lite 一致性无序预测，与 pLDDT 结果完全吻合：')
w('')
w('| 基因 | MobiDB-lite 判定的无序区（aa） | 覆盖 |')
w('|---|---|---|')
for g in ORDER:
    regs = MOB.get(g, [])
    if not regs:
        w(f"| {g} | 无 | — |")
    else:
        s = '; '.join(f"{x['start']}–{x['end']}" for x in sorted(regs, key=lambda y: y['start']))
        tot = len({i for x in regs for i in range(x['start'], x['end'] + 1)})
        w(f"| {g} | {s} | {tot} aa |")
w('')
w('（MobiDB-lite 只在前 60 aa 上跑过，见 6.1；因此"1–60 全段无序"实际是"已检测部分全段无序"。）')
w('')
w('### 3.4 结论')
w('')
for g in DIS:
    r = S[g]
    un = r['fraction_very_low_lt50'] + r['fraction_low_50_70']
    extra = f"，pTM 仅 {f(r['ptm'],3)}（pTM<0.5 说明没有可靠的链内相对取向）" if r.get('ptm') else ''
    w(f"- **{NSHORT[g]}**：pLDDT<70 的残基占 **{un*100:.1f}%**"
      f"（其中 <50 的 {f(r['fraction_very_low_lt50']*100,1)}%），pLDDT 均值仅 {f(r['mean'],1)}{extra}。")
w('')
w('- 两个 LEA_1 候选尤其极端：`2HG0197000` 在两个模型里**没有任何一个残基同时达到 pLDDT≥70**'
  '（有序核心 = 0 个残基），`7HG0637810` 也只有 18 个。'
  '`LEA_1`(Pfam PF03760) 本身是"偏亲水、无固定二级结构"的基序，不是可折叠的结构域——'
  '**"注释到 LEA_1 结构域"与"该区段无序"并不矛盾**。')
w('- DHN2 是唯一一个能看出一点结构的：K 片段（富含赖氨酸的重复）在部分模型里形成一段短螺旋，'
  '但只有约 10 个残基进入 pLDDT≥70，其余全为柔性链，与脱水素"高度亲水、高度无序"的经典特征一致。')
w('')
w('**结论**：对 DHN/LEA/HVA1 这类蛋白，结构预测的价值不是给出一个折叠，'
  '而是**定量地证明它没有折叠**；本次结果给出了可写进论文的无序性定量证据。')
w('')
w('---')
w('')
w('## 4. 有序 / 无序的对照意义')
w('')
ord_mean = [S[g]['mean'] for g in ORD_STRUCT]
dis_mean = [S[g]['mean'] for g in DIS]
ord_un = [S[g]['fraction_very_low_lt50'] + S[g]['fraction_low_50_70'] for g in ORD_STRUCT]
dis_un = [S[g]['fraction_very_low_lt50'] + S[g]['fraction_low_50_70'] for g in DIS]
w(f'有序组 5 个靶标的 pLDDT 均值 {f(min(ord_mean),1)}–{f(max(ord_mean),1)}、<70 残基占比 '
  f'{f(min(ord_un),2)}–{f(max(ord_un),2)}；'
  f'无序组 4 个为均值 {f(min(dis_mean),1)}–{f(max(dis_mean),1)}、<70 占比 '
  f'{f(min(dis_un),2)}–{f(max(dis_un),2)}——**两项指标上两组均不重叠**；'
  '说明这套流程具备区分能力，无序组的低分不是流程噪声或预测失败。')
w('')
w('---')
w('')
w('## 5. 与 CD-Search / InterPro 结构域注释的一致性检查')
w('')
w('### 5.1 方法')
w('')
w('把每个靶标的预测结构与同源参考结构做 **Smith-Waterman 局部比对**（BLOSUM62，gap −11/−1），')
w('在比对上的残基对上做 Kabsch 叠加，报两个指标：')
w('')
w('- `RMSD_core`：只取**两个模型都判为高置信**（pLDDT≥70）的残基做整体叠加；')
w('- `best_60aa`：核心区内最佳叠合的 60 残基窗口的 RMSD——对多结构域蛋白与螺线管蛋白，'
  '这直接回答"这段折叠本身对不对"，不受结构域间相对取向差异干扰。')
w('')
w('用局部比对而不是全局比对是必需的：大麦 LRR 蛋白只有 224 aa，而水稻 SERK2 有 628 aa，'
  '全局比对会把它强行拉到激酶域上，得出完全错误的结论（改用局部比对后该靶标 RMSD 从 11.66 Å 降到 1.74 Å）。')
w('')
w('### 5.2 结果')
w('')
w('> 下表的"注释到的结构域"列取自共享注释，**只覆盖前 60 aa**（见 6.1）。')
w('> "比对对象"列：8 个靶标比对的是 AFDB 上的**同源蛋白**（水稻/玉米/拟南芥，非大麦自身）；')
w('> DHN2 与 HVA1 因无独立参考，比对的是本流程对**同一序列**独立跑出的 ColabFold 模型，')
w('> 所以其序列一致性恒为 1.00——这两行的 RMSD 反映的是"两次独立预测差多少"，即构象自由度。')
w('')
w('| 靶标 | 注释到的结构域 | 比对对象 | 序列一致性 | 覆盖度 | RMSD_core | best_60aa | 判定 |')
w('|---|---|---|---|---|---|---|---|')
for g in ORDER:
    d = DC.get(g, {})
    doms = ANN.get(g, [])
    seen, ds = set(), []
    for x in doms:
        sig, desc = x['signature'], x['signature_desc']
        if sig in seen:
            continue
        seen.add(sig)
        ds.append(sig if desc in ('-', '') else f"{sig}（{desc[:30]}）")
    dstr = '; '.join(ds[:3]) if ds else '无正式结构域（仅无序预测）'
    if 'ca_rmsd_core_A' in d or 'ca_rmsd_A' in d:
        if d.get('cmp_kind') == 'AFDB homolog':
            obj = f"AFDB 同源 {d['afdb_ac']}（{d['afdb_organism'].split()[0]}）"
        else:
            obj = '本流程独立 ColabFold 模型（同一序列）'
        core, b60 = d.get('ca_rmsd_core_A'), d.get('best_60aa_window_rmsd_A')
        if core is None:
            verdict = '无序（无高置信核心）'
        elif g in DIS:
            verdict = '无序（构象不固定，RMSD 高符合预期）'
        elif core < 2:
            verdict = '一致'
        elif b60 is not None and b60 < 1:
            verdict = '折叠一致，域间取向不同'
        else:
            verdict = '需谨慎'
        w(f"| {NSHORT[g]} | {dstr} | {obj} | {f(d['seq_identity'],2)} | {f(d['coverage_of_target'],2)} | "
          f"{ags(core)} | {ags(b60)} | {verdict} |")
    else:
        w(f"| {NSHORT[g]} | {dstr} | 无同源条目（该基因无 uniprot_ac） | — | — | — | — | 仅凭折叠形态判断 |")
w('')
w('### 5.3 结构域区间内的 pLDDT')
w('')
w('| 靶标 | 注释的结构域（来源 / 区间） | 区间内平均 pLDDT | 区间内 ≥70 占比 |')
w('|---|---|---|---|')
for g in ORDER:
    for dom in DC.get(g, {}).get('domain_regions', []):
        w(f"| {NSHORT[g]} | {dom['source']} {dom['signature']} {dom['start']}–{dom['end']} "
          f"（{dom['signature_desc'][:34]}） | {f(dom['region_mean_plddt'],1)} | {f(dom['region_frac_ge70'],3)} |")
w('')
w('### 5.4 逐靶标结论与不一致点')
w('')
w('- **NAC/NAM（5H0524540）—— 一致性很好**：InterProScan 给出 Pfam PF02365 + ProSite PS51005 + '
  'SUPERFAMILY NAC + Gene3D NAC 四重支持；预测结构中 N 端折叠成紧凑的 α+β 球状域，'
  '与水稻 NAC2(AF-Q75HE5) 的高置信核心叠加 **RMSD 0.49 Å（153 残基）**，best_60aa 0.26 Å——'
  '**NAC 折叠成立**。同时 C 端 41.8% 残基 pLDDT<50、最长连续低置信 52 aa，'
  '对应 NAC 家族典型的无序转录激活区（TAR），也是该蛋白整体均值只有 67.1 的原因。'
  '这类"N 端有序 + C 端无序"正是 NAC 转录因子的真实架构，不是预测失败。')
w('- **水通道蛋白（1H0046190 TIP3-1 / 6H0596770 PIP2-2）—— 一致性很好**：'
  '核心 RMSD 分别 0.30 Å / 0.76 Å，best_60aa 0.12 / 0.16 Å，6 次跨膜螺旋束清晰。'
  'N 端有一段低置信柔性尾（MobiDB-lite 也判 1–20 / 1–22 无序），是水通道蛋白 N 端调控尾的正常特征。')
w('- **HSP90/HSP82（2H0099870）—— 折叠一致、域间取向不同**：716 aa 全长预测均值 83.7；'
  '与玉米 HSP82(AF-Q08277) 比对一致性 0.91、覆盖 0.99，但整体核心 RMSD 12.77 Å、best_60aa 只有 0.19 Å。'
  '分域拆开看：N 端 ATPase 域 7.8 Å、中间域 3.3 Å、C 端域 0.46 Å——'
  '**差异来自 N 端 ATPase 域与中间域的相对取向**（HSP90 本就是构象可变的分子开关），'
  '而不是二级/局部折叠不同。**全长预测在这里明显超出注释能力**：'
  '注释序列被截断在 60 aa，只能看到 N 端 ATPase 域的开头。')
w('- **LRR 受体激酶（3H0318380）—— 折叠一致，但只是"半截蛋白"**：'
  '预测结构 N 端帽 + 一串 LRR 重复，与水稻 SERK2(AF-Q7XV05) 的局部比对一致性 0.64、'
  '核心 RMSD 1.74 Å、best_60aa 0.22 Å，**LRR 螺线管折叠成立**。')
w('  但关键发现是：大麦这个蛋白只有 224 aa，**只对应 SERK2 的第 15–203 位**'
  '（即 LRRNT + 前几个 LRR 重复），既没有跨膜段之后的激酶域，也没有完整的 LRR 阵列。'
  '加上 SignalP/Phobius 给出信号肽与跨膜螺旋，**它更可能是基因模型不完整（或假基因），'
  '而不是一个完整的 SERK2 型受体激酶**。论文里不宜按"受体激酶"来写它的功能。')
w('- **⚠ 不一致点（重要）**：任务描述把 `HORVU.MOREX.r3.4HG0386900` 列为"NAM/NAC 结构域候选"，'
  '但 CD-Search 与 InterProScan **都没有任何 NAC 证据**：')
w('  CD-Search 只命中 `pfam09713`（E=7.7e-21），InterProScan 给出 NCBIfam TIGR01589 + '
  'PANTHER PTHR31871 + Pfam PF09713 三重支持，对应的是 '
  '**IPR006476 / CHP01589（植物保守未知功能蛋白）**，NAM 结构域（PF02365）未命中。')
w('  该蛋白仅 87 aa，预测折叠是一个紧凑的小 α/β 球状域（均值 84.7、中位数 93.8，87 个残基里只有 2 个 <50），'
  '与 NAC 折叠无关。**建议把"NAM/NAC 候选"改称"CHP01589/DUF 候选"，或从 NAC 相关结论中剔除。**')
w('')
w('---')
w('')
w('## 6. 局限与重要提醒')
w('')
w('### 6.1 共享注释表存在 60 aa 截断（影响面较大，建议上游修复）')
w('')
w('`processed/reanalysis/annotation/` 下的 InterProScan 与 CD-Search 结果，其输入 FASTA '
  '`processed/reanalysis/top100_proteins.fasta` 中 107 条序列**全部恰好 60 aa**；')
w('`iprscan_top100.tsv` 的 `sequence_length` 列 441 行全为 60，CD-Search 命中的 `To` 坐标最大值也恰为 60。')
w('即 **现有结构域注释只覆盖每个蛋白 N 端前 60 个残基**，')
w('而 `ANNOTATION_API_REPORT.md` 第 5 节"注意事项"未提及这一点。')
w('对 HSP90(716 aa)、NAC(330 aa)、水通道蛋白(263/286 aa) 这类蛋白注释是不完整的，')
w('下游引用（包括本报告 5.3 的结构域区间）都受此限制。')
w('本次全长结构预测与注释的差异大多源于此，**不是预测错误**。')
w('')
w('### 6.2 其他')
w('')
w('1. 结构预测用的是 AF2 **单体**模式（`--model-type auto` → `alphafold2_ptm`）。'
  '水通道蛋白是四聚体、HSP90 是二聚体，单体预测不能评价其组装界面。')
w('2. MSA 来自 ColabFold 公共服务器 `https://api.colabfold.com`，结果可复现但依赖该服务可用性。')
w('3. `afdb/` 下 9 个 AFDB 结构是**其他物种的同源蛋白**（其中 2 个为大麦自身），'
  '仅作折叠参照，不能当作大麦靶标的结构引用。')
w('4. 两个 LEA_1 候选与拟南芥参考的同源一致性只有 0.48 / 0.53，'
  '参考结构本身的代表性有限；其"无序"结论主要由 pLDDT 与 MobiDB-lite 两条内部证据支持。')
w('')
w('---')
w('')
w('## 7. 产出文件')
w('')
w('目录：`/ds2/workspace/ai/highland-barley-omics-gis/processed/reanalysis/structures/`')
w('')
w('| 文件 | 内容 |')
w('|---|---|')
w('| `pdb/<基因>_colabfold_rank001.pdb` | 本地 ColabFold 预测，每靶标 rank_001 最佳模型 |')
w('| `afdb/AF-<AC>-F1.pdb` + `afdb/afdb_meta.json` | AlphaFold DB 下载的 9 个条目及元数据 |')
w('| `PDB_LIST.tsv` | 全部 PDB 路径 / 来源 / 模型名对照 |')
w('| `structure_plddt_summary.png` | 10 靶标逐残基 pLDDT 曲线（2×5） |')
w('| `structure_gallery.png` | 10 靶标 3D cartoon 结构画廊（3Dmol.js + headless Chrome 离线渲染） |')
w('| `structure_gallery_catrace.png` | 备用画廊：matplotlib Cα 走向图（同色标） |')
w('| `structure_plddt_summary.tsv` | pLDDT 汇总表 |')
w('| `plddt_<基因>.tsv` | 逐残基 pLDDT 原始数据 |')
w('| `domain_consistency.tsv` / `.json` | 结构域一致性检查（含逐结构域区间 pLDDT） |')
w('| `colabfold_scores.json` / `msa_depth.json` | rank_001 的 pLDDT、pTM / MSA 深度 |')
w('| `cf_out/` | ColabFold 原始输出（5 模型 PDB、scores json、a3m、coverage 图、完整 log） |')
w('| `render/` | 3Dmol 渲染用的 HTML 与单靶标 PNG |')
w('| `*.py` | 全部分析脚本（见下） |')
w('')
w('**复现命令**')
w('')
w('```bash')
w('cd /ds2/workspace/ai/highland-barley-omics-gis/processed/reanalysis/structures')
w('# 1) ColabFold 本地 GPU 预测（.venv_cf 为已建好的环境）')
w('.venv_cf/bin/colabfold_batch input.fasta cf_out --num-recycle 3 \\')
w('    --model-type auto --host-url https://api.colabfold.com --rank plddt')
w('# 2) 汇总、画图、一致性检查')
w('python3 finalize_outputs.py && python3 collect_scores.py && python3 msa_depth.py')
w('python3 analyze_structures.py                      # pLDDT 统计（需 numpy）')
w('.venv_cf/bin/python build_domain_annotation.py     # 结构域注释汇总（需 biopython）')
w('python3 extract_mobidb.py')
w('.venv_cf/bin/python domain_check.py                # 结构比对一致性检查')
w('python3 make_figures.py && python3 render_gallery.py')
w('python3 build_report.py                            # 生成本报告')
w('```')
w('')
w('**在 ChimeraX / PyMOL 中查看**')
w('')
w('```bash')
w('chimerax pdb/*.pdb                 # 8 个本地预测')
w('chimerax afdb/AF-P12951-F1.pdb afdb/AF-P14928-F1.pdb   # DHN2 / HVA1')
w('# ChimeraX 中按 pLDDT 着色：')
w('#   color byattribute bfactor palette /alphafold   (或 palette ^RdYlBu_r)')
w('```')
w('')
w('（`structure_gallery.png` 已经是按 pLDDT 着色的 cartoon 渲染，可直接用作论文插图草稿。）')
w('')

open('STRUCTURE_REPORT.md', 'w').write('\n'.join(L))
print('wrote STRUCTURE_REPORT.md:', len('\n'.join(L)), 'chars,', len(L), 'lines')
