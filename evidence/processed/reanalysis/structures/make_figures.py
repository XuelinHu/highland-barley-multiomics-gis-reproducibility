#!/usr/bin/env python3
"""structure_plddt_summary.png (2x5 per-residue pLDDT) + structure_gallery.png (Ca trace, pLDDT colored)."""
import os, json, glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from plddt_lib import read_plddt

BASE = os.path.dirname(os.path.abspath(__file__)); os.chdir(BASE)

# --- design tokens (dataviz skill: surface / ink / recessive grid) ---
SURFACE   = '#fcfcfb'
INK       = '#1a1a1a'
INK_2     = '#555555'
INK_MUTED = '#8a8a8a'
GRID      = '#e3e3e0'
LINE      = '#0b4f9e'          # single-hue series color (sequential hue, dark step)

# AlphaFold pLDDT confidence bands - the field-standard convention, used ONLY as
# labelled background bands (text range labels carry identity, never colour alone).
BANDS = [(90, 100, '#0053D6', '>90 very high'),
         (70,  90, '#65CBF3', '70-90 confident'),
         (50,  70, '#FFDB13', '50-70 low'),
         ( 0,  50, '#FF7D45', '<50 very low')]

summary = json.load(open('structure_plddt_summary.json'))
LABEL = {
 'HORVU.MOREX.r3.4HG0386900': 'CHP01589 / DUF2017 (4H0386900)',
 'HORVU.MOREX.r3.5HG0524540': 'NAC/NAM TF (5H0524540)',
 'HORVU.MOREX.r3.3HG0318380': 'LRR-RK SERK2-like (3H0318380)',
 'HORVU.MOREX.r3.1HG0046190': 'Aquaporin TIP3-1 (1H0046190)',
 'HORVU.MOREX.r3.6HG0596770': 'Aquaporin PIP2-2 (6H0596770)',
 'HORVU.MOREX.r3.2HG0099870': 'HSP90 / HSP82 (2H0099870)',
 'HORVU.MOREX.r3.5HG0507410': 'Dehydrin DHN2 (5H0507410)',
 'HORVU.MOREX.r3.1HG0080020': 'LEA HVA1 (1H0080020)',
 'HORVU.MOREX.r3.7HG0637810': 'LEA_1 candidate A (7H0637810)',
 'HORVU.MOREX.r3.2HG0197000': 'LEA_1 candidate B (2H0197000)',
}
ORDER = list(LABEL.keys())

def load_curve(g):
    f = f'plddt_{g}.tsv'
    if not os.path.exists(f):
        return None
    x, y = [], []
    for ln in open(f).read().splitlines()[1:]:
        p = ln.split('\t')
        x.append(int(p[0])); y.append(float(p[2]))
    return x, y

# ---------------- Figure 1: per-residue pLDDT small multiples ----------------
fig, axes = plt.subplots(2, 5, figsize=(20, 7.6), facecolor=SURFACE)
axes = axes.ravel()
for ax, g in zip(axes, ORDER):
    ax.set_facecolor(SURFACE)
    rec = next((r for r in summary if r['gene_id'] == g), {})
    cur = load_curve(g)
    for lo, hi, c, _ in BANDS:
        ax.axhspan(lo, hi, color=c, alpha=0.13, lw=0, zorder=0)
    if cur:
        ax.plot(cur[0], cur[1], color=LINE, lw=2.0, solid_capstyle='round', zorder=3)
        ax.set_xlim(0.5, len(cur[0]) + 0.5)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 70, 90, 100])
    ax.set_yticklabels(['0', '50', '70', '90', '100'], fontsize=8.5, color=INK_2)
    ax.tick_params(axis='x', labelsize=8.5, colors=INK_2, length=3)
    ax.tick_params(axis='y', colors=INK_MUTED, length=3)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.grid(axis='y', color=GRID, lw=0.7, alpha=0.9, zorder=1)
    ax.set_axisbelow(True)
    src = rec.get('source', 'NA')
    if src.startswith('AFDB') and 'exact' not in src:
        sub = f"AFDB homolog {src.split(':')[1]} (reference only) · target {rec.get('seq_len','?')} aa"
    elif 'exact' in src:
        sub = f"AlphaFold DB {src.split(':')[1].split()[0]} (exact sequence) · {rec.get('seq_len','?')} aa"
    else:
        sub = f"ColabFold local (AF2) · {rec.get('seq_len','?')} aa"
    ax.set_title(f"{LABEL[g]}\n{sub}", fontsize=9.6, color=INK, pad=7, loc='left')
    if rec.get('mean') is not None:
        ax.text(0.985, 0.06,
                f"mean {rec['mean']:.1f}\nmedian {rec['median']:.1f}",
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=9.0, color=INK, zorder=5,
                bbox=dict(boxstyle='round,pad=0.32', fc=SURFACE, ec=GRID, lw=0.8))
    ax.set_xlabel('residue index', fontsize=8.5, color=INK_MUTED)
for ax in axes[::5]:
    ax.set_ylabel('pLDDT', fontsize=9, color=INK_2)

handles = [Patch(fc=c, alpha=0.30, ec='none', label=f'{lab}') for _, _, c, lab in BANDS]
handles.append(Line2D([], [], color=LINE, lw=2.0, label='per-residue pLDDT'))
fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False,
           fontsize=10, bbox_to_anchor=(0.5, -0.008), labelcolor=INK_2)
fig.suptitle('Per-residue pLDDT of 10 drought-candidate barley proteins  '
             '(AlphaFold2 / ColabFold, 5 models, 3 recycles, ranked by pLDDT; '
             'DHN2 & HVA1 from AlphaFold DB)',
             fontsize=12.5, color=INK, x=0.5, y=0.985)
fig.tight_layout(rect=[0, 0.045, 1, 0.955])
fig.savefig('structure_plddt_summary.png', dpi=200, facecolor=SURFACE, bbox_inches='tight')
plt.close(fig)
print('wrote structure_plddt_summary.png')

# ---------------- Figure 2: structure gallery (Ca trace, pLDDT coloured) ----------------
def ca_trace(pdb):
    pts = []
    for ln in open(pdb):
        if ln.startswith('ATOM') and ln[12:16].strip() == 'CA':
            pts.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54]), float(ln[60:66])))
    return pts

def band_color(v):
    for lo, hi, c, _ in BANDS:
        if lo <= v < hi or (hi == 100 and v >= 90):
            return c
    return BANDS[-1][2]

fig = plt.figure(figsize=(20, 8.4), facecolor=SURFACE)
for i, g in enumerate(ORDER):
    rec = next((r for r in summary if r['gene_id'] == g), {})
    ax = fig.add_subplot(2, 5, i + 1, projection='3d')
    ax.set_facecolor(SURFACE)
    ax.set_axis_off()
    if rec.get('pdb') and os.path.exists(rec['pdb']):
        pts = ca_trace(rec['pdb'])
        if len(pts) > 1:
            segs = [((pts[j][0], pts[j][1], pts[j][2]), (pts[j+1][0], pts[j+1][1], pts[j+1][2]))
                    for j in range(len(pts) - 1)]
            cols = [band_color((pts[j][3] + pts[j+1][3]) / 2) for j in range(len(pts) - 1)]
            ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=2.4))
            xs, ys, zs = [p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts]
            for f, arr in ((0, xs), (1, ys), (2, zs)):
                c = (min(arr) + max(arr)) / 2; r = max(max(arr) - min(arr), 1) * 0.62
                getattr(ax, f'set_{"xyz"[f]}lim')(c - r, c + r)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-62)
    m = rec.get('mean'); n = rec.get('n_residues', '?')
    ax.set_title(f"{LABEL[g]}\n{'pLDDT ' + format(m, '.1f') if m is not None else 'NA'} · {n} aa",
                 fontsize=9.6, color=INK, pad=-4, loc='center')
fig.legend(handles=[Patch(fc=c, alpha=0.95, ec='none', label=lab) for _, _, c, lab in BANDS],
           loc='lower center', ncol=4, frameon=False, fontsize=10.5,
           bbox_to_anchor=(0.5, 0.002), labelcolor=INK_2)
fig.suptitle('Structure gallery - Cα traces coloured by pLDDT (same colour scheme as AlphaFold DB)',
             fontsize=12.5, color=INK, y=0.975)
fig.tight_layout(rect=[0, 0.04, 1, 0.95])
fig.savefig('structure_gallery_catrace.png', dpi=190, facecolor=SURFACE, bbox_inches='tight')
plt.close(fig)
print('wrote structure_gallery_catrace.png')
