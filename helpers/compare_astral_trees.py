#!/usr/bin/env python3
"""
Compare ASTRAL species trees (e.g., from collapsed gene-tree sets).

Usage:
    python3 compare_astral_trees.py BS0=astral_M20_BS0.tre BS10=astral_M20_BS10.tre ... [REF=published_astral.tre]

The FIRST tree given is the baseline. Writes three tab-separated files:
  compare_RF.tsv          normalized RF distance of every tree to the baseline (and to REF if given)
  compare_key_nodes.tsv   LPP of the key clades in each tree ('absent' if the clade is not recovered)
  compare_all_splits.tsv  every bipartition found in any tree, with its LPP in each tree
Only the Python standard library is used. Node labels are read as LPP
(ASTRAL default output, e.g. ')0.98:0.12'); '[pp1=...]'-style labels are also handled.
"""
import sys, re

# ---------------- Newick parsing ----------------
def parse(s):
    s = s.strip().rstrip(';')
    pos = 0
    def node():
        nonlocal pos
        n = {'children': [], 'label': ''}
        if s[pos] == '(':
            pos += 1
            while True:
                n['children'].append(node())
                if s[pos] == ',': pos += 1
                elif s[pos] == ')': pos += 1; break
                else: raise ValueError(f'unexpected {s[pos]!r} at {pos}')
        if pos < len(s) and s[pos] == "'":          # quoted label
            end = s.index("'", pos + 1); n['label'] = s[pos+1:end]; pos = end + 1
        else:
            start = pos
            while pos < len(s) and s[pos] not in ':,()': pos += 1
            n['label'] = s[start:pos]
        if pos < len(s) and s[pos] == ':':
            pos += 1
            while pos < len(s) and s[pos] not in ',()': pos += 1
        return n
    return node()

def read_tree(path):
    with open(path) as f:
        txt = ''.join(l.strip() for l in f if l.strip())
    return parse(txt.split(';')[0])

def lpp(label):
    m = re.search(r'pp1=([0-9.eE-]+)', label)
    if m: return float(m.group(1))
    try: return float(label)
    except ValueError: return None

def splits(tree):
    """Return {frozenset(side without anchor taxon): LPP} for all non-trivial bipartitions."""
    leaves = []
    def collect(n):
        if not n['children']: leaves.append(n['label']); return
        for c in n['children']: collect(c)
    collect(tree)
    allt = frozenset(leaves)
    anchor = sorted(allt)[0]
    out = {}
    def rec(n):
        if not n['children']: return frozenset([n['label']])
        s = frozenset().union(*[rec(c) for c in n['children']])
        if 1 < len(s) < len(allt) - 1:
            side = s if anchor not in s else allt - s
            v = lpp(n['label'])
            if side not in out or out[side] is None: out[side] = v
        return s
    rec(tree)
    return allt, out

# ---------------- key clades (edit names here if needed) ----------------
B = 'Buddleja_'
CI   = {B+'asiatica_AB5', B+'subserrata_AB32'}
CII  = {B+x for x in ['colvilei_AB1','sessilifolia_AB2','forrestii_AB9','calostachya_AB29','fallowiana_AB14',
                       'davidii_AB7','albiflora_AB13','myriantha_AB12','candida_AB15','macrostachya_AB26']}
CORE3= {B+x for x in ['delavayi_AB11','delavayi_AB3','delavayi_AB30','microstachya_AB10','brachystachya_AB27',
                       'lindleyana_AB22','atropurpurea_AB28','curviflora_AB20','purdomii_AB31','yunnanensis_AB4',
                       'caryopteridifolia_AB19','officinalis_AB24','paniculata_AB25']}
CIII = CORE3 | {B+'crispa_AB8', B+'alternifolia_AB17', B+'tsetangensis_AB16'}
KEY = [
    ('Buddleja (Asian + OUT1 + OUT2)', CI|CII|CIII|{B+'spp_OUT1', B+'spp_OUT2'}),
    ('Asian Buddleja (Clades I-III)',  CI|CII|CIII),
    ('Clade I',  CI),
    ('Clade II + Clade III', CII|CIII),
    ('Clade II', CII),
    ('Clade III', CIII),
    ('B. forrestii + B. calostachya', {B+'forrestii_AB9', B+'calostachya_AB29'}),
    ('B. purdomii + caryopteridifolia/officinalis/paniculata/yunnanensis',
        {B+'purdomii_AB31', B+'caryopteridifolia_AB19', B+'officinalis_AB24', B+'paniculata_AB25', B+'yunnanensis_AB4'}),
    ('B. crispa + core Clade III', CORE3 | {B+'crispa_AB8'}),
    ('B. microstachya + B. brachystachya', {B+'microstachya_AB10', B+'brachystachya_AB27'}),
    # alternative (concatenation-tree) placements: present only if the taxon has moved
    ('ALT: B. forrestii + B. sessilifolia', {B+'forrestii_AB9', B+'sessilifolia_AB2'}),
    ('ALT: B. purdomii + lindleyana/atropurpurea/curviflora',
        {B+'purdomii_AB31', B+'lindleyana_AB22', B+'atropurpurea_AB28', B+'curviflora_AB20'}),
    ('ALT: B. crispa + B. alternifolia + B. tsetangensis', {B+'crispa_AB8', B+'alternifolia_AB17', B+'tsetangensis_AB16'}),
]

def short(side, allt):
    side = side if len(side) <= len(allt) / 2 else allt - side
    return ' + '.join(sorted(x.replace('Buddleja_', 'B.') for x in side))

def main():
    args = sys.argv[1:]
    if not args: sys.exit(__doc__)
    names, data, ref = [], {}, None
    for a in args:
        name, path = a.split('=', 1)
        allt, sp = splits(read_tree(path))
        if name.upper() == 'REF': ref = (allt, sp); continue
        names.append(name); data[name] = (allt, sp)
    base = names[0]
    allt0 = data[base][0]
    anchor = sorted(allt0)[0]
    for n in names:
        if data[n][0] != allt0: sys.exit(f'Taxon set of {n} differs from {base}')
    nint = len(allt0) - 3
    norm = lambda side: side if anchor not in side else allt0 - side

    with open('compare_RF.tsv', 'w') as f:
        hdr = ['tree', 'splits', f'RF_vs_{base}', f'nRF_vs_{base}'] + (['RF_vs_REF', 'nRF_vs_REF'] if ref else [])
        f.write('\t'.join(hdr) + '\n'); print('\t'.join(hdr))
        for n in names:
            s = set(data[n][1]); s0 = set(data[base][1])
            rf = len(s ^ s0); row = [n, str(len(s)), str(rf), f'{rf/(2*nint):.3f}']
            if ref:
                r = len(s ^ set(ref[1])); row += [str(r), f'{r/(2*nint):.3f}']
            f.write('\t'.join(row) + '\n'); print('\t'.join(row))

    print()
    with open('compare_key_nodes.tsv', 'w') as f:
        f.write('\t'.join(['clade'] + names) + '\n'); print('\t'.join(['clade'] + names))
        for lab, taxa in KEY:
            missing = taxa - allt0
            if missing: print(f'WARNING: {lab}: names not found: {missing}'); continue
            key = norm(frozenset(taxa))
            vals = []
            for n in names:
                sp = data[n][1]
                vals.append('absent' if key not in sp else (f'{sp[key]:.2f}' if sp[key] is not None else 'present'))
            f.write('\t'.join([lab] + vals) + '\n'); print('\t'.join([lab] + vals))

    every = []
    for n in names:
        for s in data[n][1]:
            if s not in every: every.append(s)
    with open('compare_all_splits.tsv', 'w') as f:
        f.write('\t'.join(['split (smaller side)'] + names) + '\n')
        for s in every:
            vals = [('-' if s not in data[n][1] else (f'{data[n][1][s]:.2f}' if data[n][1][s] is not None else '+')) for n in names]
            f.write('\t'.join([short(s, allt0)] + vals) + '\n')
    print('\nWrote compare_RF.tsv, compare_key_nodes.tsv, compare_all_splits.tsv')

if __name__ == '__main__':
    main()
