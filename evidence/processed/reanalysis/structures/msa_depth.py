#!/usr/bin/env python3
import os, json, glob
os.chdir(os.path.dirname(os.path.abspath(__file__)))
out = {}
for f in sorted(glob.glob('cf_out/*.a3m')):
    g = os.path.basename(f)[:-4]
    n = sum(1 for _ in open(f))
    out[g] = n
json.dump(out, open('msa_depth.json', 'w'), indent=1)
for k, v in out.items():
    print(f"{k:28s} MSA sequences = {v}")
