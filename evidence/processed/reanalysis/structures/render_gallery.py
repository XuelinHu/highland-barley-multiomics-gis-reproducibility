#!/usr/bin/env python3
"""Offline 3D render of every target with 3Dmol.js + headless Chrome, composed into
structure_gallery.png (2x5).  Falls back to a matplotlib Ca-trace gallery if the
headless browser render fails."""
import os, json, subprocess, shutil, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.patches import Patch

BASE = os.path.dirname(os.path.abspath(__file__)); os.chdir(BASE)
CHROME = shutil.which('google-chrome') or shutil.which('chromium')
THREEDMOL = 'assets/3Dmol-min.js'
SIZE = 640
os.makedirs('render', exist_ok=True)

LABEL = json.load(open('figure_labels.json'))
summary = {r['gene_id']: r for r in json.load(open('structure_plddt_summary.json'))}
ORDER = list(LABEL.keys())

HTML = """<!doctype html><html><head><meta charset="utf-8"><style>
html,body{margin:0;padding:0;background:__BG__;width:__S__px;height:__S__px;overflow:hidden}
#v{width:__S__px;height:__S__px;position:relative}
</style></head><body><div id="v"></div>
<script>__LIB__</script>
<script>
var PDB = __PDB__;
var cols = [[90,'#0053D6'],[70,'#65CBF3'],[50,'#FFDB13'],[0,'#FF7D45']];
function plddt(a){var b=(a.b===undefined)?0:a.b;for(var i=0;i<cols.length;i++){if(b>=cols[i][0])return cols[i][1];}return '#FF7D45';}
var v = $3Dmol.createViewer(document.getElementById('v'),{backgroundColor:'__BG__'});
v.addModel(PDB,'pdb');
v.setStyle({},{cartoon:{colorfunc:plddt, thickness:0.55, arrows:true}});
v.zoomTo(); v.zoom(0.74); v.render();
</script></body></html>"""

def render(gene, pdb_path):
    import hashlib
    tag = hashlib.md5(open(pdb_path, 'rb').read() + b'z074').hexdigest()[:8]
    out = f'render/{gene}_{tag}.png'
    if os.path.exists(out) and os.path.getsize(out) > 5000:
        return out
    lib = open(THREEDMOL).read()
    html = (HTML.replace('__LIB__', lib)
                .replace('__PDB__', json.dumps(open(pdb_path).read()))
                .replace('__BG__', '#fcfcfb').replace('__S__', str(SIZE)))
    hp = f'render/{gene}_{tag}.html'
    open(hp, 'w').write(html)
    cmd = [CHROME, '--headless=new', '--no-sandbox', '--disable-gpu',
           '--enable-unsafe-swiftshader', '--use-gl=angle', '--use-angle=swiftshader',
           '--hide-scrollbars', f'--window-size={SIZE},{SIZE}',
           '--virtual-time-budget=30000', f'--screenshot={out}', f'file://{os.path.abspath(hp)}']
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    return out if os.path.exists(out) and os.path.getsize(out) > 5000 else None

ok = {}
for g in ORDER:
    rec = summary.get(g, {})
    p = rec.get('pdb')
    if not p or not os.path.exists(p):
        continue
    try:
        r = render(g, p)
        if r:
            ok[g] = r
    except Exception as e:
        print('render failed', g, e)

if len(ok) < 4:
    print('3Dmol headless render unavailable -> falling back to make_figures.py Ca-trace gallery')
    sys.exit(2)

BANDS = [(90, '#0053D6', 'pLDDT >90 very high'), (70, '#65CBF3', '70-90 confident'),
         (50, '#FFDB13', '50-70 low'), (0, '#FF7D45', '<50 very low')]
SURFACE, INK, INK_2 = '#fcfcfb', '#1a1a1a', '#555555'

def source_line(rec):
    src = rec.get('source', 'NA')
    if src.startswith('AFDB') and 'exact' not in src:
        return f"AFDB homolog {src.split(':')[1]} (reference)"
    if 'exact' in src:
        return f"AlphaFold DB {src.split(':')[1].split()[0]} (exact sequence)"
    return 'ColabFold local (AlphaFold2)'

fig, axes = plt.subplots(2, 5, figsize=(20, 7.4), facecolor=SURFACE)
axes = axes.ravel()
for ax, g in zip(axes, ORDER):
    ax.set_facecolor(SURFACE); ax.set_axis_off()
    if g in ok:
        ax.imshow(imread(ok[g]))
    rec = summary.get(g, {})
    m = rec.get('mean')
    ax.set_title(f"{LABEL[g]}\n{source_line(rec)}\nmean pLDDT "
                 f"{format(m, '.1f') if m is not None else 'NA'} · {rec.get('seq_len','?')} aa",
                 fontsize=9.2, color=INK, pad=4, linespacing=1.45)
fig.legend(handles=[Patch(fc=c, ec='none', label=l) for _, c, l in BANDS],
           loc='lower center', ncol=4, frameon=False, fontsize=11,
           bbox_to_anchor=(0.5, 0.004), labelcolor=INK_2)
fig.suptitle('Structure gallery - cartoon coloured by pLDDT (AlphaFold colour scheme), '
             'rendered offline with 3Dmol.js', fontsize=12.5, color=INK, y=0.975)
fig.tight_layout(rect=[0, 0.035, 1, 0.94], h_pad=1.0)
fig.savefig('structure_gallery.png', dpi=190, facecolor=SURFACE, bbox_inches='tight')
print('wrote structure_gallery.png from', len(ok), 'renders')
