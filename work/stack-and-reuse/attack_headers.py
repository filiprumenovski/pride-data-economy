"""Attack on the header-sniff finding (model-assisted review, not a human audit).

Finding attacked: tables that the filename rules could not identify are often identifiable tool
outputs by their header row; 43% of 'ad hoc only' datasets with a read top-level table
(Header-sniff/Summary.json) and 72% of those with text tables inside tables-only zips
(Arsenal/H10-zip-headers/Summary.json) contain at least one identified table.

Inputs (offline, read only): Header-sniff/Headers.jsonl + Targets.csv.gz,
Arsenal/H10-zip-headers/Zip-headers.jsonl, Archive-contents/Dataset-before-after.csv.gz.
Matching is re-run with header_sniff.classify (identical to the frozen Header-matches files).

Sampling frame for precision = identified tables that count towards the two shares:
  top-level tables read in 'ad hoc only' datasets (the 5,988-dataset denominator) + all zip
  member tables (the 615-dataset denominator).
Strata = attributed signature (a table matching both 'DIA-NN matrix' and 'DIA-NN report' is in
the report stratum; tables matched only by a PRIDE header definition form one stratum).
  Stage 1 (requested sample): 200 tables, seed 20261001, all tables in strata with < 15,
    15 in the other small strata, the remainder split equally across the 9 large strata.
  Stage 2: the 9 strata that carry the dataset shares are topped up (simple random, same seed)
    to 36 tables, the smallest n at which 36/36 correct gives a Wilson 95% lower bound >= 0.90.
  Stage 3 (rule check, not used for signature-level precision): broad single rules topped up to
    20 tables each: PD ['master','accession'] alone; DIA-NN run-column rule; PRIDE
    'MaxQuant: peptides.txt' definition alone.
Miss sample: 150 read, unidentified tables (100 top-level from 'ad hoc only' datasets, 50 zip
members), seed 20261001; table-level miss rate weighted back to each source.

Labels were assigned by Claude (model-assisted review, not a human audit) from file name + full
header row; they are stored below in LABELS / MISS_LABELS keyed by a hash of url|member.

  python attack_headers.py           # write all outputs (needs the labels below)
  python attack_headers.py --dump    # print sampled headers for review
"""
import argparse, hashlib, json, re, sys
import numpy as np, pandas as pd
from io import StringIO
from scipy.stats import norm
from common import O, SEED
from header_sniff import classify, SIGNATURES, PRIDE, RUNCOL

OUT = O / 'Attacks/headers'
NOTE = 'model-assisted review, not a human audit'
SMALL_MIN, STAGE1_N, STAGE2_N, RULE_N = 15, 200, 36, 20
# missed tool outputs that are run summaries/QC statistics or spectral libraries rather than identification results
QC_OR_LIBRARY = {'fea8a6a105', '536c865abc', '2e213b0bf8', 'd7872338cf', '538b4f8d0b', '77a87d96f4', 'fc3aa363f6'}
BIG = ['MaxQuant table', 'Proteome Discoverer export', 'PLGS export', 'PEAKS export', 'DIA-NN matrix',
       'DIA-NN report (not PRIDE-complete)', 'Spectronaut report/pivot', 'FragPipe output', 'PRIDE definition only']


def key(url, member):
    return hashlib.md5(f'{url}|{member or ""}'.encode()).hexdigest()[:10]


def rules_fired(cols):
    s = {c.lower().strip() for c in cols}
    out = [f'{n} {alt}' for n, alts in SIGNATURES.items() for alt in alts if all(x in s for x in alt)]
    out += [f'PRIDE {n}' for n, req in PRIDE if all(r in s for r in req)]
    a, b = classify(cols)
    if 'DIA-NN matrix' in b and not any(o.startswith('DIA-NN matrix') for o in out):
        out.append('DIA-NN matrix [run-column rule]')
    if 'mzTab content in a text file' in b:
        out.append('mzTab [first cell]')
    return out


def stratum(sig):
    if sig == '':
        return 'PRIDE definition only'
    if 'DIA-NN report' in sig:
        return 'DIA-NN report (not PRIDE-complete)'
    if 'Perseus' in sig:
        return 'Perseus/limma statistics'
    return sig


def rule_cell(r):
    """Broad single rules examined separately (stage 3)."""
    if r.stratum == 'Proteome Discoverer export' and all('master' in x for x in r.rules if x.startswith('Proteome')):
        return "PD ['master','accession'] only"
    if 'DIA-NN matrix [run-column rule]' in r.rules:
        return 'DIA-NN run-column rule'
    if r.stratum == 'PRIDE definition only' and r.rules == ['PRIDE MaxQuant: peptides.txt']:
        return "PRIDE 'MaxQuant: peptides.txt' only"
    return ''


def load():
    H = pd.DataFrame([json.loads(l) for l in (O / 'Header-sniff/Headers.jsonl').open()]).drop_duplicates('url', keep='last')
    H = H[H.url.isin(pd.read_csv(O / 'Header-sniff/Targets.csv.gz', usecols=['url']).url)]
    H['source'], H['member'] = 'top-level', ''
    Z = pd.DataFrame([json.loads(l) for l in (O / 'Arsenal/H10-zip-headers/Zip-headers.jsonl').open()])
    Z['source'] = 'zip member'; Z['file_name'] = Z.member.str.split('/').str[-1]
    A = pd.concat([H, Z], ignore_index=True)
    A['read'] = A.status.isin(['ok', 'partial']) & (A['columns'].str.len() > 0)
    cl = A['columns'].apply(classify)
    A['pride_def'] = cl.str[0].str.join('; '); A['signature'] = cl.str[1].str.join('; ')
    A['identified'] = (A.pride_def != '') | (A.signature != '')
    A['key'] = [key(u, m) for u, m in zip(A.url, A.member)]
    D = pd.read_csv(O / 'Archive-contents/Dataset-before-after.csv.gz', usecols=['accession', 'adhoc_only_after'])
    adhoc = set(D.loc[D.adhoc_only_after, 'accession'])
    A['frame'] = (A.source == 'zip member') | A.accession.isin(adhoc)
    return A.sort_values(['source', 'url', 'member']).reset_index(drop=True)


def draw(pool, n, rng, exclude=()):
    pool = pool[~pool.key.isin(set(exclude))]
    return pool.iloc[rng.choice(len(pool), size=min(n, len(pool)), replace=False)] if n > 0 else pool.iloc[:0]


def sample_identified(A):
    I = A[A.frame & A.read & A.identified].copy()
    I['stratum'] = I.signature.map(stratum)
    I['rules'] = I['columns'].apply(rules_fired)
    I['rule_cell'] = I.apply(rule_cell, axis=1)
    rng = np.random.default_rng(SEED)
    pop = I.stratum.value_counts()
    small = [s for s in sorted(pop.index) if s not in BIG]
    alloc = {s: min(SMALL_MIN, pop[s]) for s in small}
    rest = STAGE1_N - sum(alloc.values())
    for i, s in enumerate(sorted(BIG)):
        alloc[s] = rest // len(BIG) + (i < rest % len(BIG))
    parts = []
    for s in sorted(pop.index):
        parts.append(draw(I[I.stratum == s], alloc[s], rng).assign(stage=1))
    S = pd.concat(parts)
    for s in sorted(BIG):
        parts.append(draw(I[I.stratum == s], STAGE2_N - alloc[s], rng, S.key).assign(stage=2))
    S = pd.concat(parts)
    for c in sorted(I.rule_cell.unique()):
        if c:
            have = (S.rule_cell == c).sum()
            parts.append(draw(I[I.rule_cell == c], RULE_N - have, rng, S.key).assign(stage=3))
    S = pd.concat(parts)
    return I, S


def sample_unidentified(A):
    U = A[A.frame & A.read & ~A.identified]
    rng = np.random.default_rng(SEED + 1)
    S = pd.concat([draw(U[U.source == 'top-level'], 100, rng), draw(U[U.source == 'zip member'], 50, rng)])
    return U, S


def wilson(k, n, z=norm.ppf(0.975)):
    if n == 0:
        return np.nan, np.nan
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def hdr(cols, n=60):
    cols = [str(c) for c in cols]
    return ' | '.join(cols[:n]) + (f' | ... (+{len(cols) - n} more)' if len(cols) > n else '')


def dump(S, U):
    for r in S.itertuples():
        print(f'\n## {r.key} st{r.stage} [{r.stratum}] {r.rule_cell} {r.source} {r.accession}\n   {r.member or r.file_name}\n   rules: {r.rules}\n   {hdr(r.columns)}')
    print('\n\n######## UNIDENTIFIED')
    for r in U.itertuples():
        print(f'\n## {r.key} {r.source} {r.accession}\n   {r.member or r.file_name}\n   {hdr(r.columns)}')


def parse_labels(text, cols):
    return pd.read_csv(StringIO(text.strip()), names=cols, dtype=str, keep_default_na=False, skipinitialspace=True).set_index('key')


PD_COMPANION = re.compile(r'(# ?psms|number of psms|protein fdr confidence|sum pep score|score sequest|score mascot|score ms amanda|# peptides|'
                          r'number of peptides|coverage \[%\]|coverage in percent|exp\. q-value|abundance|found in sample|# aas|number of aas|calc\. ?pi)')
PEAKS_COMPANION = re.compile(r'^(#peptides|#unique|avg\. mass|protein group|protein id|#spec|#feature|ascore|found by|source file|ptm|'
                             r'ppm|coverage \(%\)|denovo score|alc \(%\)|tag length|quality score|#vector)$')


def companion_screen(I):
    """Full-population check of the broad rules: how many matched tables lack any other column of the attributed tool."""
    out = {}
    x = I[I.rule_cell == "PD ['master','accession'] only"]
    k = x['columns'].map(lambda c: sum(bool(PD_COMPANION.search(str(s).lower())) for s in c))
    out["PD ['master','accession'] only"] = {'tables': int(len(x)), 'with_no_other_PD_column': int((k == 0).sum()),
                                             'with_fewer_than_3_other_PD_columns': int((k < 3).sum())}
    y = I[I.stratum == 'PEAKS export']
    k = y['columns'].map(lambda c: sum(bool(PEAKS_COMPANION.search(str(s).lower().strip())) for s in c))
    out["PEAKS ['-10lgp']"] = {'tables': int(len(y)), 'with_fewer_than_3_other_PEAKS_columns': int((k < 3).sum()),
                               'note': 'low-companion tables are PEAKS layouts with short headers, mostly DB search psm.csv '
                                       '(Peptide/-10lgP/Mass/ppm/m/z/RT/Scan/Accession/PTM) and PEAKS Q exports'}
    z = I[I.stratum == 'PLGS export']
    out['PLGS'] = {'tables': int(len(z)), 'with_protein.key': int(z['columns'].map(lambda c: 'protein.key' in {s.lower() for s in c}).sum())}
    r = I[I.rule_cell == 'DIA-NN run-column rule']
    nm = r.file_name.str.contains('matrix|diann|dia-nn|dia_nn', case=False)
    odd = r[~nm]
    sn = odd['columns'].map(lambda c: 'proteinaccessions' in {s.lower() for s in c})
    out['DIA-NN run-column rule'] = {'tables': int(len(r)), 'file_name_says_matrix_or_diann': int(nm.sum()),
                                     'other_names': int((~nm).sum()),
                                     'other_names_with_Spectronaut_style_ProteinAccessions_column': int(sn.sum()),
                                     'note': 'ProteinAccessions/Genes/ProteinDescriptions + run columns look like Spectronaut pivots with the PG. '
                                             'prefix removed, i.e. probable false DIA-NN attributions; the rest are DIA-NN matrices'}
    m = I[I.stratum == 'mzTab content in a text file']
    out['mzTab first-cell rule'] = {'tables_in_frame': int(len(m)),
                                    'first_row_MTD_submitter_name': int(m['columns'].map(lambda c: len(c) > 1 and c[1].lower() == 'submitter_name').sum())}
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dump', action='store_true'); a = ap.parse_args()
    A = load(); I, S = sample_identified(A); U, M = sample_unidentified(A)
    if a.dump:
        dump(S, M); return
    OUT.mkdir(parents=True, exist_ok=True)
    L = parse_labels(LABELS, ['key', 'label', 'form', 'reason'])
    ML = parse_labels(MISS_LABELS, ['key', 'label', 'tool', 'fixable', 'reason'])
    missing = set(S.key) - set(L.index) | set(M.key) - set(ML.index)
    if missing:
        sys.exit(f'{len(missing)} sampled tables have no label: {sorted(missing)[:5]}')
    assert L.label.isin(['correct', 'wrong', 'undecidable']).all() and L.form.isin(['native', 'derived', 'na']).all()
    assert ML.label.isin(['missed tool output', 'generic processed table', 'non-result', 'undecidable']).all()

    # --- Sample-identified.csv
    S = S.join(L, on='key')
    S['header_row'] = S['columns'].map(lambda c: ' | '.join(map(str, c)))
    S['n_columns'] = S['columns'].str.len()
    S['rules_fired'] = S.rules.str.join(' ; ')
    S.insert(0, 'review', NOTE)
    cols = ['review', 'key', 'stage', 'stratum', 'rule_cell', 'source', 'accession', 'file_name', 'member', 'signature', 'pride_def',
            'rules_fired', 'label', 'form', 'reason', 'n_columns', 'header_row', 'url']
    S[cols].to_csv(OUT / 'Sample-identified.csv', index=False)

    # --- Precision by signature (stages 1+2 are a simple random sample within each stratum)
    pop = I.stratum.value_counts()
    rows = []
    def prec(name, kind, g, population, note=''):
        k, w, u = int((g.label == 'correct').sum()), int((g.label == 'wrong').sum()), int((g.label == 'undecidable').sum())
        n = len(g); lo, hi = wilson(k, n); lo2, hi2 = wilson(k, k + w)
        rows.append({'review': NOTE, 'unit': kind, 'name': name, 'population_tables_in_frame': int(population), 'n_judged': n,
                     'correct': k, 'wrong': w, 'undecidable': u,
                     'precision_strict': round(k / n, 3) if n else np.nan, 'wilson_low': round(lo, 3), 'wilson_high': round(hi, 3),
                     'precision_decidable_only': round(k / (k + w), 3) if k + w else np.nan,
                     'wilson_low_decidable': round(lo2, 3), 'wilson_high_decidable': round(hi2, 3),
                     'correct_but_derived_pct': round(100 * ((g.label == 'correct') & (g.form == 'derived')).sum() / max(k, 1), 1),
                     'note': note})
    for s in pop.index:
        g1 = S[(S.stratum == s) & (S.stage == 1)]; g12 = S[(S.stratum == s) & (S.stage <= 2)]
        prec(s, 'signature, stage 1 (200-table sample)', g1, pop[s])
        if s in BIG:
            prec(s, 'signature, stages 1+2 (n up to 36)', g12, pop[s])
    for c in sorted(x for x in I.rule_cell.unique() if x):
        prec(c, 'broad rule, all stages', S[S.rule_cell == c], (I.rule_cell == c).sum(),
             'includes stage 3 top-up; a random sample within this rule only')
    P = pd.DataFrame(rows)
    P.to_csv(OUT / 'Precision-by-signature.csv', index=False)
    # the estimate used per signature: stages 1+2 for the 9 large strata, stage 1 (= all or 15) for the small ones
    sig = P[(P.unit == 'signature, stages 1+2 (n up to 36)') | ((P.unit == 'signature, stage 1 (200-table sample)') & ~P.name.isin(BIG))].set_index('name')
    w = pop / pop.sum()
    overall = float((sig.precision_strict * w).sum())
    se = float(np.sqrt((w ** 2 * sig.precision_strict * (1 - sig.precision_strict) / sig.n_judged).sum()))
    lb_overall = float((sig.wilson_low * w).sum())   # conservative: population-weighted mean of stratum lower bounds
    pd_other_file = S[(S.stratum == 'PRIDE definition only') & (S.stage <= 2)]
    pride_mismatch = {'n': int(len(pd_other_file)), 'other_MaxQuant_file_than_definition': int(pd_other_file.reason.str.contains('file type differs').sum())}
    s1 = S[S.stage == 1]; s1w = s1.stratum.map(w / s1.stratum.value_counts())
    overall_stage1 = float((s1w * (s1.label == 'correct')).sum() / s1w.sum())

    # --- Misses
    M = M.join(ML, on='key')
    idds = {(src, acc) for src, acc in A.loc[A.frame & A.read & A.identified, ['source', 'accession']].itertuples(index=False)}
    M['dataset_has_identified_table'] = [(s_, a_) in idds for s_, a_ in zip(M.source, M.accession)]
    M['header_row'] = M['columns'].map(lambda c: ' | '.join(map(str, c)))
    M.insert(0, 'review', NOTE)
    M[['review', 'key', 'source', 'accession', 'file_name', 'member', 'dataset_has_identified_table', 'label', 'tool', 'fixable', 'reason',
       'header_row', 'url']].to_csv(OUT / 'Sample-unidentified.csv', index=False)
    miss = {}
    upop = U.source.value_counts()
    for src, g in M.groupby('source'):
        k = int((g.label == 'missed tool output').sum()); lo, hi = wilson(k, len(g))
        kf = int(((g.label == 'missed tool output') & (g.fixable == 'yes')).sum())
        h = g[~g.dataset_has_identified_table]; kh = int((h.label == 'missed tool output').sum())
        miss[src] = {'population_unidentified_tables': int(upop[src]), 'n_judged': int(len(g)),
                     'labels': g.label.value_counts().to_dict(), 'missed_tool_output_pct': round(100 * k / len(g), 1),
                     'wilson_95': [round(100 * lo, 1), round(100 * hi, 1)],
                     'missed_and_fixable_by_header_signature_pct': round(100 * kf / len(g), 1),
                     'n_from_datasets_with_no_identified_table': int(len(h)),
                     'missed_tool_output_pct_in_those': round(100 * kh / max(len(h), 1), 1),
                     'missed_tools': g.loc[g.label == 'missed tool output', 'tool'].value_counts().to_dict()}
    # dataset-weighted: each sampled table weighted by 1 / (read tables of its dataset in that source), restricted to
    # datasets with no identified table -> mean per-dataset share of tables that are missed tool outputs
    nread = A[A.frame & A.read].groupby(['source', 'accession']).size()
    for src, g in M.groupby('source'):
        h = g[~g.dataset_has_identified_table]
        wd = 1 / np.array([nread[(src, a_)] for a_ in h.accession])
        miss[src]['missed_results_excluding_QC_or_library_pct'] = round(
            100 * ((g.label == 'missed tool output') & ~g.key.isin(QC_OR_LIBRARY)).mean(), 1)
        y = (h.label == 'missed tool output').to_numpy()
        est = float((wd * y).sum() / wd.sum())
        rng = np.random.default_rng(SEED)
        bs = [float((wd[i] * y[i]).sum() / wd[i].sum()) for i in (rng.integers(0, len(h), len(h)) for _ in range(2000))]
        miss[src]['dataset_weighted_missed_share_in_unidentified_datasets_pct'] = round(100 * est, 1)
        miss[src]['dataset_weighted_bootstrap_95'] = [round(100 * float(np.percentile(bs, 2.5)), 1), round(100 * float(np.percentile(bs, 97.5)), 1)]
        base = 43.2 if src == 'top-level' else 72.0
        miss[src]['implied_share_with_a_tool_output_rough'] = {
            'pct': round(base + (100 - base) * est, 1),
            'range_from_bootstrap': [round(base + (100 - base) * float(np.percentile(bs, q)) / 1, 1) for q in (2.5, 97.5)],
            'caveat': 'mean per-dataset share of missed tool-output tables is a lower bound on P(dataset has >= 1 such table); '
                      'model-assisted labels, table-level sample of up to 20 tables per dataset; indicative only'}
    wt = upop / upop.sum()
    miss['pooled_table_weighted_missed_pct'] = round(float(sum(wt[s_] * miss[s_]['missed_tool_output_pct'] for s_ in upop.index)), 1)

    # --- Dataset-level shares under signature filters
    keep = {'all signatures (as published)': set(pop.index),
            'Wilson lower bound >= 0.90 (strict: undecidable counted as not correct)': set(sig.index[sig.wilson_low >= 0.90]),
            'point precision >= 0.90': set(sig.index[sig.precision_strict >= 0.90]),
            'drop signatures with any wrong judgement': set(sig.index[sig.wrong == 0]),
            'drop mzTab rule and DIA-NN run-column rule': set(pop.index) - {'mzTab content in a text file'}}
    A2 = A[A.frame & A.read].copy()
    A2['stratum'] = np.where(A2.identified, A2.signature.map(stratum), '')
    A2['rule_cell'] = ''
    ix = A2.identified
    A2.loc[ix, 'rule_cell'] = A2.loc[ix].assign(rules=A2.loc[ix, 'columns'].apply(rules_fired)).apply(rule_cell, axis=1)
    D = pd.read_csv(O / 'Archive-contents/Dataset-before-after.csv.gz', usecols=['accession', 'adhoc_only_after'])
    den_top = set(A.loc[(A.source == 'top-level') & A.read, 'accession']) & set(D.loc[D.adhoc_only_after, 'accession'])
    den_zip = set(A.loc[(A.source == 'zip member') & (A.status == 'ok'), 'accession'])
    shares = {}
    for name, ks in keep.items():
        t = A2[A2.stratum.isin(ks)]
        if name.startswith('drop mzTab rule and DIA-NN'):
            t = t[t.rule_cell != 'DIA-NN run-column rule']
        nt = t[t.source == 'top-level'].accession.nunique(); nz = t[t.source == 'zip member'].accession.nunique()
        shares[name] = {'signatures_kept': sorted(ks),
                        'top_level_datasets_identified': nt, 'top_level_denominator': len(den_top),
                        'top_level_pct': round(100 * nt / len(den_top), 1),
                        'zip_datasets_identified': nz, 'zip_denominator': len(den_zip), 'zip_pct': round(100 * nz / len(den_zip), 1)}
    pr = sig.precision_strict.to_dict()
    def expected(src, den):
        t = A2[(A2.source == src) & (A2.stratum != '')]
        q = t.assign(f=1 - t.stratum.map(pr)).groupby('accession').f.prod()
        return round(100 * float((1 - q).sum()) / len(den), 1)
    shares['expected share, precision-weighted (assumes independent errors)'] = {
        'top_level_pct': expected('top-level', den_top), 'zip_pct': expected('zip member', den_zip)}
    # only native (not re-saved/derived) tool exports: re-weight by stratum share of derived among correct
    nat = (S[S.stage <= 2].groupby('stratum').apply(lambda g: ((g.label == 'correct') & (g.form == 'native')).mean())).to_dict()
    shares['note'] = ('Dataset shares use the same denominators as the frozen summaries: top-level = ad hoc only datasets with a read '
                      'top-level table; zip = datasets with at least one zip member table read with status ok.')

    summ = {'review': NOTE, 'seed': SEED,
            'frame': {'identified_tables': int(len(I)), 'by_source': I.source.value_counts().to_dict(), 'by_stratum': pop.to_dict(),
                      'definition': 'identified (read) top-level tables in ad hoc only datasets + identified zip-member tables'},
            'sample_sizes': {'stage1': int((S.stage == 1).sum()), 'stage2': int((S.stage == 2).sum()), 'stage3': int((S.stage == 3).sum())},
            'labels_all_stages': S.label.value_counts().to_dict(),
            'forms_among_correct_all_stages': S.loc[S.label == 'correct', 'form'].value_counts().to_dict(),
            'native_share_by_stratum_stages_1_2': {k_: round(v, 3) for k_, v in nat.items()},
            'overall_table_weighted_precision_strict': round(overall, 4), 'overall_se': round(se, 4),
            'overall_table_weighted_precision_stage1_only': round(overall_stage1, 4),
            'overall_conservative_weighted_wilson_low': round(lb_overall, 3),
            'pride_definition_only_tables_that_are_a_different_MaxQuant_file': pride_mismatch,
            'broad_rule_full_population_screen': companion_screen(I),
            'misses': miss, 'dataset_shares': shares}
    (OUT / 'Summary.json').write_text(json.dumps(summ, indent=2, default=int))
    write_verdict(P, sig, summ)
    print(P.drop(columns=['review', 'note']).to_string()); print(json.dumps(summ, indent=1, default=int))


def write_verdict(P, sig, summ):
    """Verdict.md: numbers are filled from this run; the judgement text was written after reading the labels."""
    sh = summ['dataset_shares']; ms = summ['misses']; scr = summ['broad_rule_full_population_screen']
    full, lb = sh['all signatures (as published)'], sh['Wilson lower bound >= 0.90 (strict: undecidable counted as not correct)']
    nar = sh['drop mzTab rule and DIA-NN run-column rule']
    big = sig.loc[[b for b in BIG if b in sig.index]]
    small = sig.loc[[x for x in sig.index if x not in BIG]]
    rule = P[P.unit == 'broad rule, all stages'].set_index('name')
    t, z = ms['top-level'], ms['zip member']
    pdm, pk = scr["PD ['master','accession'] only"], scr["PEAKS ['-10lgp']"]
    def row(name, r):
        return (f"| {name} | {r.population_tables_in_frame:,} | {r.n_judged} | {r.correct} | {r.wrong} | "
                f"{r.precision_strict:.2f} ({r.wilson_low:.2f}-{r.wilson_high:.2f}) | {r.correct_but_derived_pct:.0f}% |")
    lines = [
        '# Attack on the header-sniff finding: verdict',
        '',
        f'*{NOTE.capitalize()}. Claude judged each table from its file name and full header row. Seed {SEED}. '
        'Labels and reasons are in Sample-identified.csv and Sample-unidentified.csv, so a person can check them.*',
        '',
        '**Finding attacked.** Many tables that the filename rules could not identify turn out to be tool outputs when their header row is read. '
        f"The published figures were {full['top_level_pct']}% of 'ad hoc only' datasets with a read top-level table "
        f"({full['top_level_datasets_identified']:,}/{full['top_level_denominator']:,}) and {full['zip_pct']}% of datasets with "
        f"text tables inside tables-only zips ({full['zip_datasets_identified']}/{full['zip_denominator']}). Each of these datasets contains at least one identified table.",
        '',
        '## Verdict: SURVIVES (as a lower bound), with two rules to narrow',
        '',
        '### 1. Precision (are the identified tables really that tool?)',
        f"Sampling frame: {summ['frame']['identified_tables']:,} identified tables. These are top-level tables in ad hoc only datasets plus tables inside zips. "
        f"The requested stratified sample was 200 tables, at least 15 per signature where available. A top-up raised the 9 signatures that carry the dataset shares to 36 tables each. "
        f"36 is the smallest n at which 36/36 correct gives a Wilson lower bound of at least 0.90. "
        f"In total {summ['sample_sizes']['stage1'] + summ['sample_sizes']['stage2'] + summ['sample_sizes']['stage3']} tables were judged: "
        f"{summ['labels_all_stages'].get('correct', 0)} correct, {summ['labels_all_stages'].get('wrong', 0)} wrong and "
        f"{summ['labels_all_stages'].get('undecidable', 0)} undecidable.",
        '',
        '| Signature | Tables in frame | Judged | Correct | Wrong | Precision (Wilson 95%) | Correct but re-saved/derived |',
        '|---|---|---|---|---|---|---|']
    lines += [row(n, r) for n, r in big.iterrows()] + [row(n, r) for n, r in small.iterrows()]
    lines += ['', 'Broad single rules (random within the rule; includes stage-3 top-up):', '',
              '| Rule | Tables | Judged | Correct | Wrong | Precision (Wilson 95%) | Derived |', '|---|---|---|---|---|---|---|']
    lines += [row(n, r) for n, r in rule.iterrows()]
    lines += ['',
        f"**Overall precision, table-weighted:** {summ['overall_table_weighted_precision_strict']:.4f}. "
        f"A conservative version averages the stratum Wilson lower bounds, weighted by population: {summ['overall_conservative_weighted_wilson_low']:.3f}.",
        '',
        '**Full-population screens of the broad rules.** These are mechanical checks over every matched table, not judgements:',
        f"- PD `['master','accession']` alone: all {pdm['tables']:,} tables also carry at least 3 other PD columns "
        "(Protein FDR Confidence, # PSMs, Score Sequest HT, Abundance). The rule is broad on paper but clean in practice. "
        "`['identifying node']` never fires alone; it always co-occurs with `['confidence','annotated sequence']`.",
        f"- PEAKS `['-10lgp']`: {pk['tables']:,} tables. The {pk['with_fewer_than_3_other_PEAKS_columns']} with short headers are still PEAKS layouts, "
        "mostly `DB search psm.csv`.",
        f"- PLGS: every one of the {scr['PLGS']['tables']:,} tables has `protein.key`, which is a PLGS-specific name.",
        f"- DIA-NN run-column rule: {scr['DIA-NN run-column rule']['file_name_says_matrix_or_diann']} of {scr['DIA-NN run-column rule']['tables']} tables are named "
        f"`*matrix*`/DIA-NN. {scr['DIA-NN run-column rule']['other_names_with_Spectronaut_style_ProteinAccessions_column']} have "
        "`ProteinAccessions | Genes | ProteinDescriptions` plus run columns. These look like Spectronaut pivots with the `PG.` prefix removed, so they are probably false DIA-NN attributions (about 1%).",
        f"- mzTab first-cell rule: all {scr['mzTab first-cell rule']['tables_in_frame']} tables in the frame start `MTD | submitter_name | <person>`. "
        "That is the ProteomeXchange submission-summary (`submission.px`) layout. It is neither mzTab nor a results table. "
        "Outside this frame, the rule mostly catches real mzTab sections.",
        '',
        '**Correct tool, but not always a native export.** '
        f"{summ['forms_among_correct_all_stages'].get('derived', 0)} of {summ['labels_all_stages'].get('correct', 0)} correct tables are re-saved or derived. "
        'These keep the tool\'s column vocabulary but have renamed, dropped or added columns, statistics or a supplementary-table layout. '
        'The share is about 20% for DIA-NN matrices and Spectronaut, about 11% for MaxQuant, and 100% for the limma rule. '
        'So "identified" means "carries a tool\'s column vocabulary". It does not mean "an unmodified tool export".',
        f"Also, {summ['pride_definition_only_tables_that_are_a_different_MaxQuant_file']['other_MaxQuant_file_than_definition']} of "
        f"{summ['pride_definition_only_tables_that_are_a_different_MaxQuant_file']['n']} tables matched only by a PRIDE header definition are MaxQuant files of another type. "
        '`modificationSpecificPeptides.txt` matches the `peptides.txt` definition, and `allPeptides.txt` matches the `evidence.txt` definition. '
        'The tool is right but the PRIDE file type is not, which matters for the "meets a PRIDE header definition" share.',
        '',
        '### 2. Misses (unidentified tables in ad hoc only datasets)',
        f"Top-level (n={t['n_judged']} of {t['population_unidentified_tables']:,}): {t['missed_tool_output_pct']}% (Wilson {t['wilson_95'][0]}-{t['wilson_95'][1]}%) "
        f"are recognisable tool outputs that we missed. Excluding run summaries, QC statistics and spectral libraries, the figure is {t['missed_results_excluding_QC_or_library_pct']}%. "
        f"Labels: {t['labels']}.",
        f"Inside zips (n={z['n_judged']} of {z['population_unidentified_tables']:,}): {z['missed_tool_output_pct']}% (Wilson {z['wilson_95'][0]}-{z['wilson_95'][1]}%). "
        f"Excluding QC and libraries, {z['missed_results_excluding_QC_or_library_pct']}%. Labels: {z['labels']}. "
        'Most zip non-results are PLGS logs and MaxQuant parameters.txt files.',
        'The most frequent misses are:',
        '- PD protein exports without the `Master` column (PD 1.x, R-friendly headers, `ProteinGroups`);',
        '- Mascot CSV exports, which open with a `Header` banner and a dashed line;',
        '- Comet and DTASelect files with a version banner in row 1;'
        '- ProteinPilot summaries;',
        '- MaxQuant `msms.txt` and `summary.txt`;',
        '- Byonic, pGlyco, pLink, Sipros and Progenesis QI outputs (Progenesis has a two-row header).',
        'Roughly half of all sampled unidentified tables could be caught by an additional header-row signature.',
        f"Rough dataset-level implication: among datasets with no identified table, the mean per-dataset share of missed tool-output tables is "
        f"{t['dataset_weighted_missed_share_in_unidentified_datasets_pct']}% (bootstrap {t['dataset_weighted_bootstrap_95'][0]}-{t['dataset_weighted_bootstrap_95'][1]}%) "
        f"at top level and {z['dataset_weighted_missed_share_in_unidentified_datasets_pct']}% in zips. "
        f"This suggests the true shares are at least about {t['implied_share_with_a_tool_output_rough']['pct']}% "
        f"(range {t['implied_share_with_a_tool_output_rough']['range_from_bootstrap'][0]}-{t['implied_share_with_a_tool_output_rough']['range_from_bootstrap'][1]}%) "
        f"and {z['implied_share_with_a_tool_output_rough']['pct']}% "
        f"(range {z['implied_share_with_a_tool_output_rough']['range_from_bootstrap'][0]}-{z['implied_share_with_a_tool_output_rough']['range_from_bootstrap'][1]}%). "
        'That is indicative only and should not be quoted without a human audit.',
        '',
        '### 3. Dataset-level shares using only high-precision signatures',
        '| Filter | Top-level share | Zip share |', '|---|---|---|',
        f"| All signatures (published) | {full['top_level_pct']}% ({full['top_level_datasets_identified']:,}/{full['top_level_denominator']:,}) | {full['zip_pct']}% ({full['zip_datasets_identified']}/{full['zip_denominator']}) |",
        f"| Wilson lower bound >= 0.90 (9 signatures; drops Mascot, Crux, Scaffold, limma and mzTab, mostly for small n) | {lb['top_level_pct']}% ({lb['top_level_datasets_identified']:,}) | {lb['zip_pct']}% ({lb['zip_datasets_identified']}) |",
        f"| Drop mzTab rule and DIA-NN run-column rule only | {nar['top_level_pct']}% ({nar['top_level_datasets_identified']:,}) | {nar['zip_pct']}% ({nar['zip_datasets_identified']}) |",
        '',
        'The shares move by at most 0.5 percentage points. Precision errors do not drive the finding. If anything, the missed tool outputs mean the figures understate it.',
        '',
        '## Corrected wording',
        f"> Reading header rows identifies at least one table with a recognised tool's column vocabulary in {lb['top_level_pct']}-{full['top_level_pct']}% of 'ad hoc only' datasets "
        f"with a read top-level table (n={full['top_level_denominator']:,}). The same holds for {lb['zip_pct']}-{full['zip_pct']}% of datasets with text tables inside tables-only zips "
        f"(n={full['zip_denominator']}). A model-assisted check of {summ['labels_all_stages'].get('correct', 0) + summ['labels_all_stages'].get('wrong', 0)} identified tables found the tool attribution correct in all of them except a 7-table rule. "
        "About 10% of the correct tables are re-saved or derived rather than native exports. "
        "These are lower bounds: about half of the remaining unidentified tables are also recognisable tool outputs "
        "(Mascot CSV, PD 1.x, ProteinPilot, Byonic, pGlyco, pLink, MaxQuant msms/summary), so the rules under-detect.",
        '',
        '## Signatures to drop or narrow',
        "- **Drop or narrow the `mzTab content in a text file` rule.** Require `MTD | mzTab-version` / `mzTab-mode`, or a `PRH`/`PEH`/`PSH`/`SMH` section header. "
        "Exclude `MTD | submitter_name`, which is a ProteomeXchange submission summary. In this frame it was 0/7 correct.",
        "- **Narrow the DIA-NN run-column rule.** Exclude tables with `ProteinAccessions`/`ProteinDescriptions` columns, which are Spectronaut-style, "
        "or require a DIA-NN file name (`*_matrix.tsv`). It was 20/20 in the sample, but about 1% of the population looks like Spectronaut.",
        "- **Rename or exclude `Perseus/limma statistics` as a \"tool output\".** It flags limma `topTable` statistics, which are downstream statistics, not a search or quantification tool (4 tables).",
        "- **Keep PD `['master','accession']`, PEAKS `['-10lgp']` and PLGS.** They are broad on paper but clean in this frame (36/36 each, plus the full-population screens). `['identifying node']` never fires alone.",
        "- **Report PRIDE-definition matches as \"MaxQuant table\".** Do not report them as the specific PRIDE file type unless a file-type-specific column is present "
        "(for example, `N-term cleavage window` for peptides.txt). 10/36 were modificationSpecificPeptides or allPeptides files.",
        '',
        '## Limits',
        '- These are model-assisted labels, not a human audit.',
        '- Each judgement uses only the file name and header row; file contents beyond the header were not checked.',
        "- Small strata (Mascot, Crux, Scaffold, limma) have n <= 15, so their lower bounds cannot reach 0.90 even at 100% observed precision.",
        '- Up to 20 tables were read per dataset. .xls/parquet and nested xlsx were not read.',
        '- The miss sample is table-level, so its dataset-level extrapolation is rough.',
        '- Claims here are descriptive only.']
    (OUT / 'Verdict.md').write_text('\n'.join(lines) + '\n')

# Labels below: model-assisted review (Claude), not a human audit. Judged from file name + full header row.
# key, label (correct / wrong / undecidable), form (native = a tool export as written; derived = a tool export
# re-saved with renamed/added/dropped columns or statistics; na), reason
LABELS = """
209f87705e,correct,derived,Crux percolator PSM columns (file_idx; distinct matches/spectrum; flanking aa) plus appended feature columns
024907c380,correct,native,Crux percolator.target.psms columns; file name says percolator.target.psms
ca3c93b5ff,correct,derived,Crux percolator PSM columns plus appended feature columns
a01d8b6acc,correct,derived,Crux percolator PSM columns plus appended feature columns
f39a8f6730,correct,native,Crux percolator.target.psms columns
cec0c96eb8,correct,derived,Crux percolator PSM columns plus appended feature columns
ee79b47b5f,correct,derived,Crux percolator PSM columns plus appended feature columns
d8778973d8,correct,derived,Crux percolator PSM columns plus appended feature columns
9fbef0a66c,correct,derived,Crux percolator PSM columns plus appended feature columns
893667e1ed,correct,native,Crux percolator.target.psms columns
8c15a2da54,correct,native,Crux percolator.target.psms columns
c4a164485a,correct,derived,Crux percolator PSM columns plus appended feature columns
d47bbcaab7,correct,derived,Crux percolator PSM columns plus appended feature columns
05b2bf0c98,correct,derived,Crux percolator PSM columns plus appended feature columns
acb1f0016c,correct,native,Crux percolator.target.psms columns
404573412a,correct,native,DIA-NN pg_matrix columns (Protein.Group/Protein.Ids/Protein.Names/Genes/First.Protein.Description + .raw runs) saved as xlsx
7127e47740,correct,native,DIA-NN pg_matrix with .d run columns; file name says DIANN pg_matrix
2259084157,correct,native,DIA-NN gene matrix (Genes + .raw run paths); file name gg_matrix
e66d0f1573,correct,native,DIA-NN pr_matrix (FragPipe-DIA runs DIA-NN for quantification)
9073ab01e5,correct,native,DIA-NN pg_matrix
75ddc5c0ae,correct,native,DIA-NN pr_matrix
a3ec67e0ef,correct,native,DIA-NN pg_matrix (.raw.dia run columns)
7089b14360,correct,native,DIA-NN unique_genes_matrix (Genes + .d runs)
980be5fb5b,correct,native,DIA-NN gg_matrix (Genes + .raw runs)
40f33e5848,correct,native,DIA-NN pg_matrix
455412923b,correct,native,DIA-NN gg_matrix (.wiff.dia run columns)
e6c2d96242,correct,derived,DIA-NN pg_matrix columns with run columns renamed to sample names
14725ddec4,correct,native,DIA-NN unique_genes_matrix
b3e86ec753,correct,native,DIA-NN pr_matrix
5d8b142252,correct,native,DIA-NN pg_matrix
83721f8621,correct,native,DIA-NN gg_matrix
9cff9c06fe,correct,native,DIA-NN unique_genes_matrix
7f26fde1a2,correct,native,DIA-NN pr_matrix (.wiff runs)
3f78f5f0f9,correct,native,DIA-NN pg_matrix
749a6ff4dc,correct,native,DIA-NN gg_matrix
df622efeed,correct,native,DIA-NN unique_genes_matrix
625b876333,correct,native,DIA-NN pg_matrix
6a737f1d1d,correct,derived,DIA-NN pr_matrix columns with run columns shortened
d37fb2204f,correct,derived,DIA-NN pr_matrix columns with run columns renamed to samples
d7b8abf2e9,correct,native,DIA-NN pr_matrix
07bda4ac3d,correct,derived,statistics table built on DIA-NN protein-group columns (PValue/log2FoldChange added)
4c5c8dd602,correct,native,DIA-NN pr_matrix
e1327ad5ec,correct,native,DIA-NN pg_matrix (.wiff.dia runs) saved as xlsx
d9964faf11,correct,derived,DIA-NN protein-group columns incl. Lib.PG.Q.Value/Global.PG.Q.Value plus renamed and z-scored samples
4264243f78,correct,derived,DIA-NN pg_matrix columns with runs renamed CT1..
5ebe70b7bf,correct,native,DIA-NN pr_matrix
2a66ebf873,correct,native,DIA-NN gg_matrix
acab737700,correct,native,DIA-NN pr_matrix
ec008ba9f6,correct,native,DIA-NN pg_matrix
d0c1067bea,correct,derived,DIA-NN pg_matrix columns (Protein.Ids dropped; peptide count column added)
7df260b314,correct,native,DIA-NN unique_genes_matrix
879f906368,correct,native,DIA-NN gg_matrix (run-column rule; Genes + .d runs)
d5b2bd70c5,correct,native,DIA-NN gg_matrix (run-column rule)
8c51959834,correct,native,DIA-NN unique_genes_matrix (run-column rule)
8684b4f044,correct,native,DIA-NN unique_genes_matrix (run-column rule)
7d0ae6e3ee,correct,native,DIA-NN unique_genes_matrix (run-column rule)
89f081d806,correct,native,DIA-NN unique_genes_matrix (run-column rule; .mzML runs)
09bf38ba93,correct,native,DIA-NN unique_genes_matrix (run-column rule)
68d9229065,correct,native,DIA-NN gg_matrix (run-column rule)
2573dcd008,correct,native,DIA-NN unique_genes_matrix (run-column rule)
08865b600b,correct,native,DIA-NN main report (early DIA-NN 1.6/1.7 column set: File.Name/Precursor.Id/Protein.Q.Value/Gene.Group.Quantity)
a3bcd270e6,correct,native,DIA-NN main report (early DIA-NN 1.6/1.7 column set: File.Name/Precursor.Id/Protein.Q.Value/Gene.Group.Quantity)
976f89273c,correct,native,DIA-NN main report (early DIA-NN 1.6/1.7 column set: File.Name/Precursor.Id/Protein.Q.Value/Gene.Group.Quantity)
148150ff90,correct,native,DIA-NN main report (early DIA-NN 1.6/1.7 column set: File.Name/Precursor.Id/Protein.Q.Value/Gene.Group.Quantity)
72dbf4885b,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
dde3e16f83,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
fc25c80267,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
469bb78f80,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
9266c8d8e5,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
2d45d562a0,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
108b0a21aa,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
c8acbb2c28,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
e781f57990,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
d94f7d21e9,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
a66a473ac5,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
48accf6cd0,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
8b0afe438c,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
7c65475c8a,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
643541bed4,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
ac024a6673,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
113fe702ec,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
845d25ef78,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
8aa11fa45e,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
1645ab1a2b,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
fb22c0abf2,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
5e96f6c2d8,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
1e1da25850,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
bf9bc8a715,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
01a1de2883,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
f65d9918a3,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
5fcabaa868,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
bb76513a5b,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
3917cb2e0a,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
c15c473ec1,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
a8885ef264,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
a017130a4f,correct,native,DIA-NN main report (File.Name/Run/Precursor.Id/PG.MaxLFQ/Lib.Q.Value)
2556ff37c3,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
2a9835ca97,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
c35298f35b,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
bb6519a965,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
2e13ca2245,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
269f13f924,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
59885b9aad,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
83390cf423,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
0191920ba3,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
7b47d0c51c,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
ae25642a42,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
ad3374a49e,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
1b3a86181c,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
98ef1753b9,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
a591cb24a7,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
2cbbdf6eb2,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
0a7983fae3,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
c4f6bea87c,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
245fac77dc,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
9018b7a290,correct,native,FragPipe/MSFragger psm.tsv columns (Spectrum File/Hyperscore/PeptideProphet Probability/Mapped Proteins)
65d28ea7f6,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
7f5153be82,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
0a3b8bd8f0,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
834c9d96af,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
f5b013eed4,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
fca23650f2,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
134350599c,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
44e5edaa09,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
e158ff1e49,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
94ec7d2fa9,correct,native,FragPipe/Philosopher protein.tsv columns (Protein ID/Protein Probability/Top Peptide Probability/Razor Spectral Count)
9c96b07785,correct,native,FragPipe combined_protein.tsv columns (Combined Spectral Count/MaxLFQ Intensity/Indistinguishable Proteins)
be71ea56e6,correct,native,FragPipe combined_protein.tsv columns (Combined Spectral Count/MaxLFQ Intensity/Indistinguishable Proteins)
de8dc15b7d,correct,native,FragPipe combined_protein.tsv columns (Combined Spectral Count/MaxLFQ Intensity/Indistinguishable Proteins)
7476157482,correct,native,FragPipe combined_protein.tsv columns (Combined Spectral Count/MaxLFQ Intensity/Indistinguishable Proteins)
2c9356503c,correct,native,FragPipe combined_protein.tsv columns (Combined Spectral Count/MaxLFQ Intensity/Indistinguishable Proteins)
edbb9a4906,correct,derived,FragPipe combined_protein columns re-saved as xlsx with a subset of quantity columns
a68f2755e5,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
0368dded1b,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
c0205e1dae,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
44f14586c4,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
4f3d2e295a,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
1b0446b45d,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
38cf053c28,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
00518510c1,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
84ebdb8d13,correct,native,Mascot CSV export (query_number/pep_seq/pep_expect/prot_hit_num; F###### Mascot result number in file name)
ea3fb7391a,correct,native,Mascot CSV protein-peptide export (prot_hit_num/prot_acc/pep_seq/pep_scan_title)
cef1cb5c6d,correct,native,Mascot CSV protein-peptide export (prot_hit_num/prot_acc/pep_seq/pep_scan_title)
f9d050036f,correct,native,Mascot CSV protein-peptide export (prot_hit_num/prot_acc/pep_seq/pep_scan_title)
a605c80eae,correct,native,Mascot CSV protein-peptide export (prot_hit_num/prot_acc/pep_seq/pep_scan_title)
b498c42f59,correct,derived,Mascot export columns (prot_hit_num/pep_seq/pep_expect) re-saved as supplementary xlsx
8dc2094292,correct,derived,Mascot export columns (prot_hit_num/pep_seq/pep_expect) re-saved as supplementary xlsx
03d6efa108,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
acd66ba327,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
d22781983a,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
b2dc5cde9c,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
3d7942e9ff,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
3a7f18ae0a,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
74882bc172,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
cf56a55338,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
236ed604e7,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
7ab0276865,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
290c042e91,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
ffb84ef5f0,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
95b60db4a3,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
58bb573744,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
7e0d1bd00c,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
4f697fcedd,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
ce52c43849,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
06df1f2c0d,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
cbd9ea409c,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
5387f727f8,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
216ac8dfb0,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
b187b2703c,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
78b16436f8,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
b9e7de889d,correct,native,MaxQuant modification Sites table (Positions within proteins/Localization prob/Score diff)
307c52d8b3,correct,native,MaxQuant modification Sites table (Positions within proteins/Localization prob/Score diff)
ab02f6dab2,correct,native,MaxQuant modification Sites table (Positions within proteins/Localization prob/Score diff)
fcdd3da7fb,correct,derived,MaxQuant proteinGroups columns (subset) re-saved as xlsx
3f18f80e3d,correct,derived,MaxQuant proteinGroups processed in Perseus (Significance B; GO annotation columns)
0a41dacdb7,correct,derived,MaxQuant proteinGroups with an added TriTrypDB annotation column saved as csv
09f261a770,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
551927b752,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
688dd02f05,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
54c204db8d,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
1d2dba78b4,correct,native,MaxQuant proteinGroups.txt columns (Majority protein IDs/Razor + unique peptides/Only identified by site)
aa981fc196,correct,derived,MaxQuant proteinGroups columns with summary columns removed (PRIDE-submission copy)
9d5d1fa639,correct,native,MaxQuant evidence.txt columns (Leading razor protein/Raw file/Match time difference)
ab502d4f9a,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
f51ece7bad,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
11100a086f,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
b3075080f9,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
35fc51d3df,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
94d422add4,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
6c38e964e1,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
93f093f882,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
9835a4e7c6,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
1311171d48,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
3f2f8f2b32,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
64cacbc0bc,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
a552b83280,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
a4a188839d,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
6cd512085d,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
2b7394391a,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
17dc0903da,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
62a1f3dbb1,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
f5f13b461e,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
e68fd3727d,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
232557ddbe,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
4d64c3c70a,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
0fe4e5a7cc,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
b54bd774e2,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
4cf6475ce1,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
e0cdce3386,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
937f8c3a55,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
d48866f7b8,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
8ef5112814,correct,native,PEAKS Studio peptide/psm export (Peptide/-10lgP/Source File/#Spec or #Feature/AScore)
91f04e8a88,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
c0c69bffe9,correct,native,PEAKS Studio protein-peptides export (Protein Group/-10lgP/#Spec or #Feature/AScore)
6da039c07d,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
99f26472b7,correct,native,PEAKS Studio proteins export (Protein Group/Protein ID/-10lgP/#Peptides/#Unique/Avg. Mass)
7b6610974d,correct,native,PEAKS Studio peptide export (Peptide/-10lgP/Source File/#Feature/AScore/Found By)
77c8ab2b35,correct,native,PEAKS Studio peptide export (Peptide/-10lgP/Source File/#Feature/AScore/Found By)
8bde783696,correct,native,PEAKS Studio peptide export (Peptide/-10lgP/Source File/#Feature/AScore/Found By)
60fdf5860e,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
041d485986,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
e852156efa,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
0e08f0541b,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
e48dfbdd66,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
7fbfbb4874,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
96072e3f4e,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
ea78401ed5,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
cfb249f6c8,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
0777fee309,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
ef0bf8355c,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
26f52d466b,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
021d42019c,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
f96494e95a,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
1a9f8f88ec,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
78bd811764,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
5dab38b62d,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
9ef34ec9f7,correct,native,PLGS (ProteinLynx Global Server) IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
f8e25f696d,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
1916ec8846,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
e72ffae877,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
d51f758c06,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
9ec11fcb7c,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
9d0391e43e,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
634f28d96d,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
bdc2e7f527,correct,native,PLGS IA_final_fragment export (protein.key/peptide.seq + product-ion columns)
ac214c86d8,correct,native,PLGS IA_final_fragment export (protein.key/peptide.seq + product-ion columns)
cdf634e373,correct,native,PLGS IA_final_fragment export (protein.key/peptide.seq + product-ion columns)
32ba52c646,correct,native,PLGS IA_final_fragment export (protein.key/peptide.seq + product-ion columns)
a1be9e566b,correct,native,PLGS IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
bbf7c7da96,correct,native,PLGS IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
fcdf24075b,correct,native,PLGS IA_final_peptide export (protein.key/protein.fmolOnColumn/peptide.seq/peptide.ModelRT)
e5d034860c,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
39ebd8b6ed,correct,native,PLGS IA_final_protein export (protein.key/protein.AutoCurate/protein.Key_ForHomologs)
1f7c03f458,correct,native,PLGS IA_final_protein columns in a file renamed by participant/condition
e649cbe6b7,correct,native,PLGS IA_final_fragment export
dbd74cf887,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
bd8c132b5d,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
6ddec3d96a,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
5cffc843c9,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
3079bb217f,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
f0a7b4756f,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
c817c1fcad,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
7070cc2c9e,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
cddad200c2,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
796f691608,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
94dc262687,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
3dcf53cdf9,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
4dcfae1495,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
257a514232,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
14691741d2,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
ee73102ab2,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
72ff274dd0,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
b3ee485993,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
4a443d1e94,correct,native,MaxQuant peptides.txt (older MaxQuant 1.x column set: Slice/Leading Razor Protein/Ratio H/L)
71a8f39e20,correct,native,MaxQuant peptides.txt (older column set with Ratio M/L)
d2f51debb5,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
765276a77e,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
62a057b265,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
8d6c94e07f,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
9cbfbdcc51,correct,native,MaxQuant allPeptides.txt (Uncalibrated m/z/Isotope pattern/MSMS isotope indices); PRIDE evidence/peptides definitions matched but file type differs
8319ab79d8,correct,native,MaxQuant allPeptides.txt (Uncalibrated m/z/Isotope pattern/MSMS isotope indices); PRIDE evidence/peptides definitions matched but file type differs
05cedc0180,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
aa7bfe6703,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
cd6ffd6e1c,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
33abf693b7,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
c2b6dbf81d,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
7877509b4b,correct,native,MaxQuant peptides.txt columns (cleavage windows/amino-acid counts/Leading razor protein/Unique (Groups))
ba31d9fe09,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
0d8cbc6206,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
5b1a1f8f00,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
9d1a85a0d0,correct,native,MaxQuant modificationSpecificPeptides.txt (Mass Fractional Part/Protein Groups/Delta score); PRIDE peptides.txt definition matched but file type differs
0c9185d25a,correct,derived,limma topTable columns (logFC/AveExpr/t/P.Value/adj.P.Val/B); statistics package not a search engine
3c8c1403f7,correct,derived,limma topTable columns (logFC/AveExpr/t/P.Value/adj.P.Val/B); statistics package not a search engine
0d6bae457a,correct,derived,limma topTable columns (logFC/AveExpr/t/P.Value/adj.P.Val/B); statistics package not a search engine
bb4b7ac7f7,correct,derived,limma statistics (logFC/adj.P.Val) on DIA-NN protein-group columns; both attributions (DIA-NN matrix and limma) apply
24730a73c5,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
f2e3e3de2a,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
c02d7b816a,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
72e0c32c13,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
6535d8b8ba,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
2a5a0929c3,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
e1368c039a,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
082acfd75a,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
8c031d4f1a,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or Mascot/Abundance)
1bcbc3dc66,correct,native,Proteome Discoverer Proteins export with R-friendly headers (Master/Accession/Number of PSMs/Score Sequest HT)
e0d9149aa9,correct,native,Proteome Discoverer Proteins export with R-friendly headers and reordered columns (Master/Protein FDR Confidence Combined/Score Sequest HT)
63b0b41308,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence/Spectrum File)
2655cd2e20,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence/Spectrum File)
6b0e623d5f,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence/Spectrum File)
7357daf8ab,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence/Spectrum File)
6d0ff0cf07,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence/Spectrum File)
642328ffa0,correct,native,Proteome Discoverer PeptideGroups export (Annotated Sequence/Qvality PEP/Master Protein Accessions/by Search Engine Sequest HT)
b42aaf98f2,correct,native,Proteome Discoverer PeptideGroups export (Annotated Sequence/Qvality PEP/Master Protein Accessions/by Search Engine Sequest HT)
fc9d7c52f9,correct,native,Proteome Discoverer PeptideGroups export (Annotated Sequence/Qvality PEP/Master Protein Accessions/by Search Engine Sequest HT)
b5f05cd500,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
1bd1043029,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
172f174939,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
690c97c93a,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
3cd231d5fa,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
0fcb021fae,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
d420fb892f,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
968edef9ef,correct,native,Proteome Discoverer Proteins export (Master/Accession plus Protein FDR Confidence/# PSMs/Score Sequest HT or MS Amanda/Abundance)
a944d5c258,correct,native,Proteome Discoverer Proteins export with R-friendly headers (Master/Protein FDR Confidence Combined/Score Sequest HT)
49c7a3897d,correct,derived,Proteome Discoverer Proteins columns re-saved as csv with an added index column
5ba21f27ab,correct,native,Proteome Discoverer PeptideGroups export (Annotated Sequence/Master Protein Accessions/Sequest HT columns)
e0edc23baa,correct,native,Proteome Discoverer PeptideGroups export (Annotated Sequence/Master Protein Accessions/Sequest HT columns)
e0fe19a7fc,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence)
d748895a18,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence)
4d8b9e7676,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence)
dfa9023923,correct,native,Proteome Discoverer PSMs export (Confidence/Identifying Node/PSM Ambiguity/Annotated Sequence)
2d71880305,correct,native,Proteome Discoverer PSMs export with R-friendly headers (Identifying Node/MS Amanda/CharmeRT)
2250788266,correct,native,Scaffold Peptide Report (Biological sample name/Protein identification probability/Exclusive unique peptide count)
b372be22a6,correct,native,Scaffold Peptide Report (Biological sample name/Protein identification probability/Exclusive unique peptide count)
3aa8a586ec,correct,native,Scaffold Peptide Report (Biological sample name/Protein identification probability/Exclusive unique peptide count)
a3f3a6469e,correct,native,Scaffold Protein Report (Protein identification probability/Exclusive unique peptide count)
453f60d036,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
8fdf55d04a,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
6a8461343c,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
7f0d6e3a8e,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
6a3ad094fe,correct,native,Scaffold Peptide/Samples report (Biological sample category/Protein identification probability/Best Mascot Ion score)
aa35dd7f55,correct,native,Scaffold Samples Report text export (Protein Identification Probability/Exclusive Unique Peptide Count blocks per sample)
7a819a5b2b,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
b322dabb8c,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
fe5645ab93,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
777ef7706c,correct,native,Scaffold Spectrum Report (Protein identification probability/Peptide identification probability/SEQUEST XCorr)
e5a1b6d9b1,correct,native,Scaffold Peptide Report (Biological sample name/Protein identification probability/Exclusive unique peptide count)
aa5eff1e47,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
588ad7a516,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
16d84afe0c,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
ecd84fe245,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
72fbfff4b5,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
8dd5cc3cc8,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
a78b478f2e,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
28e1710fe7,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
010c18cb51,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
565f97648f,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
e1e98268c7,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
a902c3c985,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
1f9c30fa9d,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
e36790c22c,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
816a493134,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
c1b2119d70,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
6cc906321f,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
dbe5584ab9,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
ccc855e308,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
2605a5db3d,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
6be9c704a4,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
181f23792f,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
14a4eee2b4,correct,native,Spectronaut pivot report (PG./PEP./EG. prefixed columns with [n] run.raw/.d/.htrms.PG.Quantity or EG.TotalQuantity)
2416189101,correct,native,Spectronaut long-format report (R.FileName/R.Condition/PG./PEP./EG. columns; BGS Factory Report)
f51e106a5b,correct,native,Spectronaut long-format report (R.FileName/R.Condition/PG./PEP./EG. columns; BGS Factory Report)
bce1742089,correct,native,Spectronaut long-format report (R.FileName/R.Condition/PG./PEP./EG. columns; BGS Factory Report)
74e06b1a68,correct,native,Spectronaut long-format report (R.FileName/R.Condition/PG./PEP./EG. columns; BGS Factory Report)
3397133f66,correct,native,Spectronaut long-format report (R.FileName/R.Condition/PG./PEP./EG. columns; BGS Factory Report)
49d0bebba0,correct,derived,Spectronaut protein-group columns (PG.ProteinGroups/PG.Genes) with log2-processed renamed sample columns
98f434d206,correct,derived,differential statistics (log2fc/padjust) built on Spectronaut PG. columns
d47118bb71,correct,derived,Spectronaut PG.ProteinAccessions/PG.Genes columns with unnamed sample columns (supplementary xlsx)
9fdd1f1376,correct,derived,Spectronaut PG./EG.PrecursorId columns with run columns renamed to samples
7ca8fcfe89,correct,derived,PG.ProteinGroups (Spectronaut naming) plus renamed sample columns and protein/gene names; weakly Spectronaut-derived
0150c69576,correct,derived,Spectronaut PG. columns with renamed samples; file name says SN16 (Spectronaut 16)
bdba814ccc,correct,derived,Spectronaut PG./PEP. columns with renamed samples (Report_... file)
62c9ddaa07,correct,derived,MSstats-style summary (BioReplicate/nFeatures) of Spectronaut output; file name says Spectronaut
209a50a31b,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
fe90ec9c74,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
cfb6f64d76,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
c5837f9d1f,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
8c0bbd73fd,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
7ee46d3e36,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
bbc089669e,wrong,na,first row is 'MTD submitter_name <person>': ProteomeXchange submission-summary (submission.px) format; not mzTab and not a results table
"""

# key, label (missed tool output / generic processed table / non-result / undecidable), tool, fixable by a
# header-row signature (yes/no), reason
MISS_LABELS = """
629d6e6ab5,undecidable,,no,first row is a 'Search Name:' key-value line; no column header to judge
94d09e34f1,missed tool output,Proteome Discoverer (by file name),no,header row empty in first sheet; file name 'MultiConsensus_proteingroups' is PD vocabulary
51a34c6a7f,non-result,,no,free-text readme line ('includes data of the following paper')
ae87f9f283,missed tool output,Comet,yes,first line 'CometVersion 2019.01' is the Comet txt-output banner
24790e3db9,missed tool output,X!Tandem (GPM),yes,peptide/modifications/charge/expect/scan/intensity/proteins from a .t.xml converted to text
b4bf07486e,missed tool output,DTASelect,yes,first line 'DTASelect v2.0.39' banner of a DTASelect-filter file
1ee2de9885,missed tool output,Mascot,yes,'Header' then dashed line is the opening of a Mascot CSV export
31b951cf42,generic processed table,,no,PSM-level table (Seq/Mod/z/score/FDR/Raw scan) from an unidentified in-house tool
cf4106a6f3,generic processed table,,no,protein quant table (Collapsed/Parsimony/scaled/sum/mean/median) from an unidentified pipeline
480609411d,non-result,,no,raw-file to sample mapping
f1dbfaa1f0,missed tool output,Spectronaut,yes,run columns '[1] X.raw.PG.Quantity' / PG.IsSingleHit are Spectronaut pivot names; PG. prefix stripped from leading columns
69a732f631,non-result,MaxQuant (probable),no,Parameter/Value parameter file (likely MaxQuant parameters.txt); not a results table
9bb23375e7,non-result,,no,TMT sample reference / file nomenclature sheet
73d1f81c9d,missed tool output,Peptide-PRISM on PEAKS de novo,yes,PEAKS de novo columns (Denovo score/ALC (%)/Tag length) plus Peptide-PRISM annotation (Genome/ORF/netMHC)
dcea422352,non-result,,no,README text
3c853f7180,generic processed table,,no,supplementary results table whose first row is a caption
5c22709978,undecidable,,no,first row is ' | 5'; header not readable
fea8a6a105,missed tool output,MaxQuant,yes,MaxQuant summary.txt (Raw file/Enzyme/MS/MS identified/Peaks sequenced); run summary rather than results
a31282bd0d,non-result,MaxQuant (input),no,MaxQuant experimental design template (Name/Fraction/Experiment/PTM)
633a0d2e98,missed tool output,Proteome Discoverer,yes,PD ProteinGroups export (Checked/Protein Group ID/# Proteins/# PSMs/Group Description)
8ecabcaa8b,non-result,,no,FASTA sequence
49e3126959,generic processed table,,no,PSM-level table (Seq/Mod/z/score/FDR/Raw.scan) from an unidentified in-house tool
249d83b9ab,missed tool output,Byonic (probable),yes,first sheet starts 'Summary:' as Byonic xlsx reports do; no column header in row 1
d5fea34e16,generic processed table,,no,SEQUEST-based TMT PSM table with link columns (ScanF/XCorr/Redun Link) from an in-house pipeline; no public-tool signature
6198bd6060,missed tool output,Byonic,yes,file name says Byonic; first sheet starts 'Summary:'
db8726f18c,missed tool output,pLink,yes,Peptide_Order/Peptide/Peptide_Mass/Modifications/Proteins is the pLink 2 cross-link report; file name Looplink
64b6d62898,generic processed table,,no,PSM table (Scan/Reversed/Expectation/Theoretical m/z) whose tool could not be identified
b908b0852f,generic processed table,,no,supplementary cross-link table whose first row is a caption
5f49f1cbdb,missed tool output,FragPipe,yes,FragPipe combined_site_STY report (Index/Best Localization Probability/MaxLFQ Intensity)
ff0316f757,missed tool output,ProteinPilot,yes,ProteinPilot peptide summary (N/Unused/Total/%Cov(95)/Conf/Sc/iTRAQ ratio columns)
14aad13dad,missed tool output,OpenMS MetaProSIP,yes,RIA/LR of RIA/INT/Cor. and Global Peptide LR columns are MetaProSIP output
9092459388,missed tool output,Byonic (probable),yes,'<raw>.raw.xlsx' search output whose first sheet starts 'Summary:' as Byonic reports do
24e21a5cf7,missed tool output,Sipros,yes,.pro2pep.txt with BalancedSpectrumCounts columns is Sipros Ensemble output
6f4e62b3e8,generic processed table,,no,matrix of sample columns only
43ce393a98,non-result,,no,file list with sample descriptions
fe99bd8a0c,missed tool output,MSstats input (upstream not named),yes,MSstats input format (ProteinName/PeptideSequence/FragmentIon/IsotopeLabelType/BioReplicate/Run)
1a51125151,missed tool output,pFind,yes,pFind spectra columns (File_Name/Scan_No/Exp.MH+/Raw_Score/Final_Score/Avg.Frag.Mass.Shift)
a4ac10bec0,generic processed table,,no,supplementary cross-link table whose first row is a caption
f6e966a555,missed tool output,Mascot,yes,'Header' then dashed line is the opening of a Mascot CSV export
80af9d3a32,generic processed table,,no,protein matrix with sample columns only
44df07e97e,missed tool output,Mascot,yes,'Header' then dashed line is the opening of a Mascot CSV export (F###### file name)
49650fd872,undecidable,,no,header row empty in first sheet (file name 'Protein_Groups_Report')
25a6970b5d,non-result,,no,sample annotation sheet
fc9e614179,missed tool output,X!Tandem (rTANDEM),yes,sample-prefixed rTANDEM columns (.pep.id/.expect.value/.Iraw)
5d5e5ce2a9,missed tool output,Proteome Discoverer,yes,PD Proteins export with the Master column removed (Protein FDR Confidence: Combined/Abundance Ratio/Score Sequest HT)
d8a8ca49ff,non-result,,no,raw file to sample/search-result mapping
b691a92c7c,missed tool output,Progenesis QI,no,Progenesis QI protein measurements (Normalized abundance/Raw abundance/Spectral counts banner row); real header is row 2
289376bc49,non-result,,no,contact/metadata form (E-mail)
43a6f9b508,missed tool output,Spectronaut (derived),yes,PEP.IsProteotypic and 'SN' in file name: Spectronaut peptide table with renamed columns
536c865abc,missed tool output,FragPipe EasyPQP library,yes,spectral library TSV (PrecursorMz/ProductMz/LibraryIntensity/NormalizedRetentionTime/PrecursorIonMobility); a library not results
ac7477d3ee,missed tool output,pLink,yes,pLink 2 cross-linked sites report (Protein_Order/Unique_Peptide_Number/Spectrum_Number); file name says pLink
09eeee6464,generic processed table,,no,supplementary dataset whose first row is a caption
5a26db1992,generic processed table,,no,no header row: first row is data (peptide sequence with TAIR annotation)
69e4b678f0,missed tool output,Mascot,yes,banner '(Generated by Mascot version ...)'
01aeda7574,missed tool output,Byonic,yes,file name says Byonic; first sheet starts 'Summary:'
b1b356490a,non-result,MaxQuant (probable),no,Parameter/Value parameter file (likely MaxQuant parameters.txt)
e025c5a0d6,generic processed table,,no,log10 abundance matrix (Protein/Accession/sample columns)
5690024951,missed tool output,StrucGP,yes,structure-specific glycopeptide columns (LowEnergy_MS2Scan/Structure_coding/Corefucose/Bisection/GP_FDR) of StrucGP
64bd025f5a,missed tool output,pGlyco,yes,pGlyco GP-FDR-Pro columns (GlySpec/PepSpec/GlyIonRatio/GlycanFDR)
de0b78467b,non-result,MaxQuant (probable),no,Parameter/Value parameter file (likely MaxQuant parameters.txt)
7820b73630,missed tool output,pGlyco,yes,pGlyco3 GP-FDR-Pro-Quant-Site columns (GlySpec/PepSpec/GlycanPEP/LocalizedSiteGroups)
2ed5f83a48,missed tool output,BioPharma Finder,yes,peptide-mapping export (Level/No./Identification/Confidence Score/Best Average Structural Resolution/ID Type)
1bf59a79a6,undecidable,,no,first readable row is a dashed separator; tool not determinable from it
f9bb12b0fa,non-result,FragPipe (input),no,FragPipe experiment_annotation.tsv (file/sample/condition/replicate)
267a796e9b,missed tool output,Sipros,yes,.pro.txt with BalancedSpectrumCounts columns is Sipros Ensemble output
17be5ca1e9,missed tool output,pGlyco,yes,pGlyco3 glycopeptide columns plus reporter-ion columns
4b1d0adafc,missed tool output,SCIEX PeakView,yes,'X.wiff (sample 1))' protein quantity columns are the PeakView SWATH export
fcb7dcaa34,generic processed table,,no,first row holds only sample names (multi-row header) of a quantification report
9abc2410d7,generic processed table,,no,in-house TMT site table (site_id/best_search_id/default~126_sn_sum); no public-tool signature
0bf5ea0939,missed tool output,DTASelect,yes,DTASelect peptide columns (Conf%/SpR/IonProportion/Redundancy/DeltCN)
78768522b8,missed tool output,Peptide-PRISM on PEAKS de novo,yes,PEAKS de novo columns (ALC (%)/Tag length) plus Peptide-PRISM annotation
71eb723fa8,missed tool output,Byonic,yes,Byonic spectrum table (PID/Prot. Rank/|Log Prob|/PEP 2D/Off-By-X/Glycans)
000a42bb13,missed tool output,ProteinPilot,yes,ProteinPilot ProteinSummary (N/Unused/Total/%Cov(95)/PVal/EF/LowerCI)
2ed1609d54,missed tool output,MaxQuant,yes,MaxQuant msms.txt (Scan number/Matches/Intensity coverage/All modified sequences)
e76efb4d2e,non-result,,no,run to sample annotation
2e213b0bf8,missed tool output,FragPipe EasyPQP library,yes,spectral library TSV (PrecursorMz/ProductMz/LibraryIntensity/PrecursorIonMobility); a library not results
6bf2a45390,non-result,,no,bare protein sequence
b16026a365,generic processed table,,no,supplementary table whose first row is a caption
c3e74f0dc9,generic processed table,,no,cross-link quantification table (Linker Modification/Best Coupling Site/XIC Area Sum) from a tool I could not identify
10769672ee,missed tool output,Proteome Discoverer,yes,PD 1.x protein export (Accession/Score/Coverage/# Proteins/# PSMs/Area/calc. pI) without Master column
aa01e7bbc1,generic processed table,,no,glycan list whose first row is a caption
2e0c9730ca,generic processed table,,no,no header: list of peptide sequences
88defb0cd2,undecidable,SCIEX PeakView (by file name),no,first row ' | 5'; file name '..._SWATH_library__FDR.xlsx' suggests a PeakView FDR report but header not readable
d86819ff1c,generic processed table,,no,SEQUEST-based PSM table with link columns (ScanF/XCorr/Redun Link) from an in-house pipeline
c44a742e86,generic processed table,,no,protein table (Amount Proteins/Amount quanted PSMs/MS1 precursor area/tmt10plex) from a tool I could not identify
07c1522457,missed tool output,Comet,yes,first line 'CometVersion 2018.01' is the Comet txt-output banner
b622ac7d79,missed tool output,DTASelect,yes,first line 'DTASelect v1.9' banner; file name DTASelect-filter
c9f5971946,generic processed table,,no,supplementary table whose first row is a caption
b9791bab7a,missed tool output,Proteome Discoverer,yes,PD Proteins export with Master removed and abundances renamed (Exp. q-value: Combined/Sum PEP Score/Score Sequest HT)
0ec5c78080,undecidable,,no,no header: first row is a glycopeptide data row
12cc958064,non-result,,no,ORF/transcript annotation used as a search database
a2e5ca368b,generic processed table,,no,protein table (Group_ID/Probability_Score/Unique_PSMs/Coverage_PCT) from a tool I could not identify
6a71ef17b5,missed tool output,ProteinPilot,yes,ProteinPilot ProteinSummary (N/Unused/Total/%Cov(95)/PVal/EF)
2391a382f1,non-result,Thermo Xtract,no,deconvoluted mass list ('Program Version:' Xtract banner); spectra-level not results
76eec90b16,undecidable,,no,header row empty
6aa47c2035,missed tool output,Proteome Discoverer,yes,PD 1.x protein export with Light/Heavy ratio columns (Accession/Score/# PSMs/Area/calc. pI)
daf7624567,missed tool output,ECL (cross-link search),yes,*.inter.target.csv columns (C13_correction/delta_C/protein_annotation_1/e_value) of ECL
33627b09ce,missed tool output,Scaffold,yes,banner 'Peptide report created on ...' of a Scaffold text export
59ad33fd46,missed tool output,PEAKS,yes,PEAKS de novo peptides (Tag length/ALC (%)/local confidence (%)); no -10lgP column
9e74283117,missed tool output,MaxQuant,yes,MaxQuant msms.txt (Scan number/Matches/Intensity coverage/All modified sequences)
eb31f3f922,missed tool output,GlycoDash (by path),no,GlycoDash summary file whose first row is a 'Parameter Settings' banner
d7872338cf,missed tool output,Proteome Discoverer,yes,PD ResultStatistics export (Name/Arithmetic Mean/Median/SD/IQR/CV) in a 'PD output' folder; QC statistics not identifications
c1b5e02390,non-result,Spectronaut (settings),no,Spectronaut report setup file (version banner); settings not results
d4a1fad8e6,undecidable,,no,no header: first row is numbers
a3506b1ceb,non-result,MaxQuant,no,MaxQuant parameters.txt (Parameter/Value)
d8a986fe85,missed tool output,Progenesis QI,no,Progenesis QI protein measurements (Normalized abundance/Raw abundance banner row); real header is row 2
38cfb26148,non-result,PLGS (log),no,PLGS IA_Log file (binary-looking content)
1f2193a6b2,missed tool output,DirectMS1 (ms1searchpy),yes,DirectMS1 PFMs table (seqs/md/mods/iorig/Isotopes/rt_diff/qpreds)
538b4f8d0b,missed tool output,MaxQuant,yes,MaxQuant summary.txt (Raw file/Enzyme/MS/MS Identified/Peaks Sequenced); run summary rather than results
10302da323,missed tool output,pLink,yes,pLink 2 report (Order/Title/Peptide_Type/Refined_Score/SVM_Score/isComplexSatisfied)
ffa480a3f2,missed tool output,pLink,yes,pLink 2 filtered proteins report (Protein_Order/Unique_Peptide_Number/Spectrum_Number)
9ca7a1a27c,missed tool output,GlycoDash (by path),no,GlycoDash summary file whose first row is a 'Parameter Settings' banner
c07b0b3b06,missed tool output,Mascot,yes,'Header' then dashed line is the opening of a Mascot CSV export
dfaf64d9bc,missed tool output,mokapot (MS2Rescore),yes,mokapot protein table (mokapot protein group/mokapot q-value/mokapot PEP)
4cad9c028f,undecidable,DIA-NN (by file name),no,read returned non-text bytes; file name pg_matrix.tsv suggests DIA-NN but header unreadable
c54f70ae03,missed tool output,PLGS,yes,PLGS IA_final_protein columns (protein.Key_ForHomologs/AutoCurate) mis-split by the delimiter guess
999936fc12,non-result,MaxQuant,no,MaxQuant parameters.txt (Parameter/Value)
012f90d093,missed tool output,MODplus/MODi,yes,modp.ranking.id output (SpectrumFile/ObservedMW/Probability/Modification) in a 'Modi' folder
377a5dfbe1,missed tool output,Mascot,yes,'Header' then dashed line is the opening of a Mascot CSV export
27dd5b1332,non-result,PLGS (log),no,PLGS IA_Log text
3a224b85c2,non-result,PLGS (log),no,PLGS Apex3D log
77a87d96f4,missed tool output,DIA-NN,yes,DIA-NN stats.tsv (Precursors.Identified/MS2.Mass.Instability); run QC not results
61b77ea7ed,missed tool output,Pytheas,yes,RNA modification search report (5'-end/3'-end/Score (Sp)/dSp/isotopologue match) of Pytheas
9e2781e548,non-result,PLGS (log),no,PLGS IA_Log text
14d941c060,non-result,MaxQuant,no,MaxQuant parameters.txt (Parameter/Value)
3af37844b8,non-result,PLGS (log),no,PLGS IA_Log text
205eec1450,undecidable,,no,first row '[Header]' only
35da21119f,missed tool output,pLink,yes,pLink 2 mono-linked peptides report (Peptide_Order/Peptide_Mass)
09045cdd7d,non-result,PLGS (log),no,PLGS IA_Log text
622294cd6c,missed tool output,PEAKS,yes,PEAKS de novo candidates (Tag length/ALC (%)/local confidence (%)); no -10lgP column
3b259f25aa,missed tool output,Proteome Discoverer,yes,PD ProteinGroups export with R-friendly headers (Protein Group ID/Number of PSMs/Group Description)
f5228c1568,missed tool output,Progenesis QI,no,Progenesis QI peptide ion measurements (Normalized abundance/Raw abundance banner row); real header is row 2
59933c1cd2,missed tool output,DIA-NN (derived),yes,DIA-NN pr_matrix columns with dots replaced by spaces (Protein Group/Precursor Id/Proteotypic)
4d6a97ab53,missed tool output,MSPathFinder (Informed Proteomics),yes,IcTda.tsv top-down PrSM table (Ms1Features/#MatchedFragments/SpecEValue)
819a4d51e8,missed tool output,NeoDisc pipeline,yes,immunopeptidomics identifications (ms_pipeline/HLAType/Max.XCorr/Min.Comet.lFDR/Max.Hyperscore) of the NeoDisc pipeline
56569d8fcd,non-result,,no,peak list with MGF-style COM= comment line
1e32f641ee,non-result,,no,inclusion list (RT/z/m/z/t start/t stop)
1625d2251c,non-result,MaxQuant,no,MaxQuant parameters.txt (Parameter/Value)
830b4469f7,missed tool output,MaxQuant (derived),yes,MaxQuant phospho-site table columns (Position within protein/Number of p(STY)/p(STY) probabilities) with renamed ratio columns
9d111badda,non-result,PLGS (log),no,ProteinLynx memory log
32257899c1,non-result,PLGS (log),no,ProteinLynx automation log
1b52e78012,non-result,PLGS (log),no,PLGS IA_Log text
9909942c8b,missed tool output,pLink,yes,pLink 2 loop-linked sites report (Protein_Order/Unique_Peptide_Number)
2ab0db3122,missed tool output,Spectrum Mill,yes,Spectrum Mill tagSummary (backbone_cleavage_score/deltaForwardReverseScore/searchCycle)
89cf9401d6,non-result,,no,list of modification mass shifts (search parameter input)
61a6158e01,non-result,DIA-NN (log),no,DIA-NN log banner
42c9eca069,non-result,MaxQuant,no,MaxQuant parameters.txt (Parameter/Value)
4072d96515,generic processed table,,no,protein-protein interaction network (prot1/prot2/score) from a downstream analysis
c1a6daf7eb,non-result,,no,KEGG annotation README (garbled encoding)
fc3aa363f6,missed tool output,MaxQuant,yes,MaxQuant summary.txt (Raw file/Enzyme/MS/MS Identified/Peaks Sequenced); run summary rather than results
"""

if __name__ == '__main__':
    main()
