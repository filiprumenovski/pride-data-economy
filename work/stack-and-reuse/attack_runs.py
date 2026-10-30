"""Adversarial test of ledger finding H1 (experiment scale; arsenal_scale_instruments.py).

Claim under attack (Ledger H1.2/H1.3): median raw runs per dataset fell from 48 (2013) to 16 (2025);
the share with >=100 runs fell 28.9% -> 8.7%; groups seen in 2015-18 and 2022-25 did not shrink;
datasets from new groups are smaller.

Alternatives tested (each written to its own CSV; verdicts in Verdict.md):
  A0  anchor: tiny 2013 cohort; pooled 2013-16 vs 2023-25; whole distribution; release vs submission year
  A1  splitting: one study deposited as several datasets (publication components from DOI/PMID/OpenAlex,
      same-group submission episodes, group-years; datasets per paper over time)
  A2  run-counting artefacts: archives, multi-file runs, loose folder contents, fractionation naming;
      hand audit of 30 + 30 random datasets (file lists printed to Audit-file-lists.csv)
  A3  composition: species/human, subfield, DIA, TMT, submission type, vendor, unpublished, PRD/PAD,
      very small (1-3 run) deposits; reweighting to the 2013-16 mix
  A4  geography: per-country trends, excluding China / newly active countries, within-country shift-share
  A5  embargo censoring: large studies may still be private at the snapshot (fixed release window)
  A6  turnover: within-group slope (group fixed effects), group-entry cohorts, newcomer gap like-for-like

Descriptive only. Run rule identical to arsenal_scale_instruments.py (re-implemented and checked equal).
CIs: 95% percentile intervals from resampling hashed pi_group clusters (accession when missing),
seed 20261001; single rates via common.cluster_boot, medians and differences via the same group resampling.
Writes only to outputs/Stack-and-reuse/Attacks/runs/.
"""
import json, re, sys
import numpy as np, pandas as pd
from common import R, O, SEED, cluster_boot
from layer_rules import RULES
sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Attacks/runs'; OUT.mkdir(parents=True, exist_ok=True)
c = connect(); c.execute("set threads=1; set memory_limit='2GB'")
S = {'notes': ['Descriptive only; no causal or predictive claim.',
               'Run counts: lower bound (runs_lb) unless stated; exact = datasets without unopened bundles.',
               'CIs: group bootstrap over hashed pi_group clusters (accession when missing), seed 20261001.']}
YEARS = list(range(2013, 2026))
EARLY, LATE = (2013, 2016), (2023, 2025)


def w(df, name):
    df.to_csv(OUT / name, index=False)


def wmedian(x, wt):
    x = np.asarray(x, float); wt = np.asarray(wt, float); k = wt > 0
    x, wt = x[k], wt[k]
    if len(x) == 0: return np.nan
    o = np.argsort(x, kind='stable'); x, wt = x[o], wt[o]; cw = np.cumsum(wt); h = cw[-1] / 2
    i = np.searchsorted(cw, h)
    # midpoint when the half-weight falls exactly on a boundary (matches pandas median for unit weights)
    if np.isclose(cw[i], h) and i + 1 < len(x): return (x[i] + x[i + 1]) / 2
    return x[i]


def gboot(df, fn, reps=1000, seed=SEED):
    """Point estimate and 95% interval of fn(df, weights), resampling clusters df.g with replacement."""
    g = pd.factorize(df.g)[0]; ng = g.max() + 1; rng = np.random.default_rng(seed)
    est = np.asarray(fn(df, np.ones(len(df))), float)
    out = np.array([fn(df, np.bincount(rng.integers(0, ng, ng), minlength=ng)[g].astype(float)) for _ in range(reps)], float)
    return est, np.nanquantile(out, .025, axis=0), np.nanquantile(out, .975, axis=0)


def rate_ci(mask, groups):
    m = np.asarray(mask, float)
    if len(m) == 0: return (np.nan,) * 3
    lo, hi = cluster_boot(m, np.asarray(groups), reps=2000, seed=SEED)
    return m.mean(), lo, hi


def period_contrast(df, col='runs_lb', reps=1000, label=None, early=EARLY, late=LATE, ycol='release_year'):
    """Early (2013-16) vs late (2023-25): medians, >=100 shares and their differences, group bootstrap."""
    z = df[df[ycol].between(*early) | df[ycol].between(*late)].copy()
    z['late'] = z[ycol] >= late[0]
    e = z.late.to_numpy(); v = z[col].to_numpy(float)
    if (~e).sum() < 10 or e.sum() < 10:
        return {'label': label, 'n_early': int((~e).sum()), 'n_late': int(e.sum())}

    def fn(_, wt):
        me, ml = wmedian(v[~e], wt[~e]), wmedian(v[e], wt[e])
        se = (wt[~e] * (v[~e] >= 100)).sum() / wt[~e].sum(); sl = (wt[e] * (v[e] >= 100)).sum() / wt[e].sum()
        return [me, ml, ml / me, se, sl, sl - se]
    est, lo, hi = gboot(z, fn, reps=reps)
    k = ['median_early', 'median_late', 'median_ratio', 'ge100_early', 'ge100_late', 'ge100_diff']
    r = {'label': label, 'early': f'{early[0]}-{early[1]}', 'late': f'{late[0]}-{late[1]}', 'n_early': int((~e).sum()), 'n_late': int(e.sum()),
         'groups_early': int(z[~e].g.nunique()), 'groups_late': int(z[e].g.nunique())}
    for i, kk in enumerate(k):
        r[kk] = est[i]; r[kk + '_lo'] = lo[i]; r[kk + '_hi'] = hi[i]
    return r


# ============================================================ data: features + run counts (re-implemented rule)
d = pd.read_parquet(O / 'Deposit-features.parquet')
d['g'] = d.pi_group.fillna('acc:' + d.accession)
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
RUNEXT = "('raw','d','wiff','wiff2','lcd','uimf')"
c.execute(rf"""create temp table nf as
select project_accession a, file_name fn, lower(file_name) n,
       lower(regexp_replace(file_name, '\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
       coalesce(file_category_value, '') cat, coalesce(file_size_bytes, 0) sz
from pride_files
where project_accession in (select accession from study)
  and coalesce(primary_public_location_value, '') not like '%/generated/%'""")
c.execute(f"create temp table nfl as select *, {case} layer from nf")
c.execute(rf"""create temp table vr as
with h as (select a, fn, n, cat, sz, regexp_replace(n, '(\.(zip|gz|bz2|7z|tar|rar|xz|tgz))+$', '') s1 from nfl where layer = 'vendor_raw'),
k as (select *, n <> s1 comp, regexp_extract(s1, '\.([a-z0-9_]+)$', 1) ext from h)
select *,
  case when regexp_matches(s1, '\.wiff2?\.(scan|mtd)$') then 'companion'
       when ext in {RUNEXT} or regexp_matches(s1, '\.d\.\d+\.dar$') then 'run'
       when ext in ('baf','tdf_bin','yep','fid') then 'run_loose'
       when ext in ('t2d','t2d_') then 'maldi_spot'
       when comp then 'bundle' else 'other' end rrole,
  case when regexp_matches(s1, '\.d\.\d+\.dar$') then regexp_replace(s1, '\.d\.\d+\.dar$', '')
       when ext in {RUNEXT} then regexp_replace(s1, '\.[a-z0-9_]+$', '') else s1 end stem
from k""")
runs = c.execute("""select a accession,
  count(distinct stem) filter (where rrole = 'run') runs_named,
  count(*) filter (where rrole = 'run_loose') runs_loose,
  count(*) filter (where rrole = 'bundle') bundles,
  count(*) filter (where rrole = 'run' and ext in ('wiff','wiff2')) f_wiff,
  count(*) filter (where rrole = 'run' and ext = 'd') f_d,
  count(*) filter (where rrole = 'run' and ext = 'raw') f_raw,
  count(*) filter (where rrole = 'run' and comp) f_run_compressed,
  sum(sz) bytes_raw_layer
  from vr group by 1""").df()
d = d.merge(runs, on='accession', how='left', validate='one_to_one')
for col in runs.columns[1:]: d[col] = d[col].fillna(0)
d['runs_exact'] = d.runs_named + d.runs_loose
d['runs_lb'] = d.runs_exact + d.bundles
d['exact'] = (d.bundles == 0) & (d.runs_exact > 0)
# check against the frozen per-dataset table of the claim under attack
ref = pd.read_csv(O / 'Arsenal/H1-H3-scale-instruments/H1-dataset-scale.csv.gz', usecols=['accession', 'runs_lb', 'runs_exact', 'status'])
chk = d[['accession', 'runs_lb', 'runs_exact']].merge(ref, on='accession', suffixes=('', '_ref'))
S['rule_reimplementation_check'] = {'datasets': len(chk), 'runs_lb_mismatch': int((chk.runs_lb != chk.runs_lb_ref).sum()),
                                    'runs_exact_mismatch': int((chk.runs_exact != chk.runs_exact_ref).sum())}
assert S['rule_reimplementation_check']['runs_lb_mismatch'] == 0
d = d.merge(ref[['accession', 'status']], on='accession')
R_ = d[(d.bytes_raw_layer > 0) & (d.runs_lb > 0)].copy()     # the claim's population (35,141)
S['population'] = {'cohort': len(d), 'with_runs': len(R_), 'excluded_no_runs': len(d) - len(R_), 'groups': int(R_.g.nunique())}
geo = json.loads((R / 'work/stack-and-reuse/ne_50m_admin_0_countries.geojson').read_text())
cont = {}
for f in geo['features']:
    p = f['properties']
    for k in ('ADMIN', 'NAME', 'NAME_LONG'): cont[p[k]] = p['CONTINENT']
R_['continent'] = R_.country.map(cont)
R_['published'] = R_.work_id.notna()
print(S['population'], S['rule_reimplementation_check'])

# ============================================================ A2a  random audit sample: print raw file lists
rng = np.random.default_rng(SEED)
smp = []
for lab_, (a, b) in {'2013-15': (2013, 2015), '2023-25': (2023, 2025)}.items():
    pool = R_[R_.release_year.between(a, b)].accession.sort_values().to_numpy()
    smp += [(lab_, x) for x in rng.choice(pool, 30, replace=False)]
smp = pd.DataFrame(smp, columns=['window', 'accession'])
c.register('smp', smp)
fl = c.execute("""select s.window, f.a accession, f.fn file_name, f.layer, coalesce(v.rrole, '-') rrole, coalesce(v.stem, '') stem,
  round(f.sz / 1e6, 1) mb from nfl f join smp s on s.accession = f.a
  left join vr v on v.a = f.a and v.fn = f.fn order by s.window, f.a, f.layer, f.fn""").df()
w(fl, 'Audit-file-lists.csv')

# hand audit (one reviewer, file lists in Audit-file-lists.csv). count: ok = every listed vendor file maps to one
# injection and companions are not double counted; bundle = count is a lower bound of unknown tightness.
# frac: explicit = names show gel bands / SCX / fraction indices; unclear = serial numbering that may be fractions.
AUDIT = {
 'PXD000131': ('ok', 'explicit', 'SCX F33.. fractions x proteases x 2 instruments (.raw and .wiff); wiff.scan/.mtd not double counted'),
 'PXD000197': ('ok', 'unclear', 'serial AP-MS run numbers'),
 'PXD000249': ('bundle', 'unclear', 'single T123.tar.gz; runs inside unknown'),
 'PXD000326': ('ok', 'no', '8 samples x 3 replicates (.wiff + .wiff.scan)'),
 'PXD000362': ('ok', 'no', 'conditions x replicates'),
 'PXD000375': ('ok', 'no', '3 strains x 3 replicates'),
 'PXD000389': ('ok', 'no', 'replicates on two dates'),
 'PXD000487': ('ok', 'explicit', '6 lanes x 20 gel bands = 120 runs'),
 'PXD000654': ('ok', 'unclear', 'serial 1..66'),
 'PXD000700': ('bundle', 'unclear', 'single AG_RAW.zip; runs inside unknown'),
 'PXD000853': ('ok', 'no', 'patient samples'),
 'PXD000916': ('ok', 'unclear', 'a-d suffixes per sample'),
 'PXD000918': ('ok', 'explicit', 'whole-proteome b01..b09 fractions x 2 replicates + 2 diGly runs'),
 'PXD000954': ('ok', 'no', '6 SWATH .wiff with companions'),
 'PXD001001': ('ok', 'no', 'conditions D/U'),
 'PXD001032': ('ok', 'no', '1 .RAW + 3 zips each named after one .wiff'),
 'PXD001162': ('ok', 'explicit', '2 samples x 44 2D fractions (not caught by name rule)'),
 'PXD001246': ('ok', 'no', 'membrane/soluble samples'),
 'PXD001264': ('ok', 'no', '.wiff + .wiff.scan'),
 'PXD001422': ('ok', 'unclear', 'cyt1a..cyt4b per donor'),
 'PXD001490': ('ok', 'no', 'samples'),
 'PXD001495': ('ok', 'unclear', 'instrument sequence numbers'),
 'PXD001641': ('ok', 'no', 'fibre samples'),
 'PXD001726': ('ok', 'no', 'samples'),
 'PXD002051': ('ok', 'no', 'samples'),
 'PXD002075': ('bundle', 'no', '2 MALDI digest zips; count plausible'),
 'PXD002170': ('ok', 'explicit', 'mixed method runs incl. iTRAQ_fractionation_1..3'),
 'PXD002244': ('ok', 'no', 'samples'),
 'PXD002308': ('ok', 'no', 'samples'),
 'PXD002320': ('ok', 'unclear', '330 plate wells (P1A01..)'),
 'PXD022018': ('bundle', 'unclear', 'single Bruker raw-data zip; runs inside unknown'),
 'PXD029215': ('ok', 'no', 'co-IP x 3'),
 'PXD031077': ('ok', 'no', 'cells/EP samples; some re-acquired with timestamp suffix (counted twice, correctly as injections)'),
 'PXD032959': ('ok', 'no', '2 x 3'),
 'PXD035055': ('ok', 'no', 'sequence numbers'),
 'PXD035863': ('ok', 'no', 'patient samples'),
 'PXD036869': ('ok', 'unclear', 'TMT sets x 3 (_001.._003)'),
 'PXD036899': ('ok', 'no', '2 samples'),
 'PXD038231': ('ok', 'no', '2 samples'),
 'PXD039817': ('ok', 'explicit', 'several studies; TMT-phospho F01..F12 and fraction01.. sets'),
 'PXD040610': ('ok', 'no', '6 samples'),
 'PXD040947': ('bundle', 'unclear', '5 .rar packs of DDA/DIA pools; runs inside unknown (likely many)'),
 'PXD043202': ('ok', 'no', '3 x 3'),
 'PXD044191': ('ok', 'no', 'trypsin/GluC per sample'),
 'PXD044637': ('ok', 'no', 'Bruker .d.7z, one per run'),
 'PXD046332': ('ok', 'no', '5 samples'),
 'PXD047557': ('ok', 'no', '8 samples'),
 'PXD051929': ('ok', 'no', 'Waters .raw.zip folders, one per run'),
 'PXD052299': ('ok', 'no', '2 x 3'),
 'PXD054309': ('ok', 'no', '12 samples'),
 'PXD057432': ('ok', 'no', '6 samples'),
 'PXD058654': ('ok', 'no', '6 samples'),
 'PXD059562': ('ok', 'explicit', 'label-free S1..S10 plus TMT sets p1..p9 / s'),
 'PXD059981': ('bundle', 'no', '6 zips named per sample (N1-3, S1-3); count plausible'),
 'PXD061302': ('ok', 'explicit', 'TMT Fraction_1..15 plus 15 later runs'),
 'PXD063359': ('ok', 'unclear', 'serial 1..14'),
 'PXD063593': ('bundle', 'unclear', '23 zips named *_raw_files; several runs each likely'),
 'PXD064007': ('ok', 'unclear', 'serial 001..024'),
 'PXD066893': ('ok', 'explicit', '2 TMT18 plexes x 6 fractions (36 samples in 12 runs)'),
 'PXD068555': ('ok', 'no', 'Bruker .d.zip, one per run'),
}
aud = smp.merge(R_[['accession', 'release_year', 'runs_lb', 'runs_exact', 'bundles', 'vendor_single']], on='accession')
aud[['count_check', 'fractionation', 'note']] = [AUDIT[a] for a in aud.accession]
w(aud, 'A2-hand-audit-60.csv')
S['A2_hand_audit'] = {wd: {'datasets': len(x), 'count_ok': int((x.count_check == 'ok').sum()), 'bundle_lower_bound': int((x.count_check == 'bundle').sum()),
                           'fractionation_explicit': int((x.fractionation == 'explicit').sum()), 'fractionation_unclear': int((x.fractionation == 'unclear').sum()),
                           'median_runs_lb': float(x.runs_lb.median())} for wd, x in aud.groupby('window')}

# ============================================================ A0  anchor, distribution, time axis
base = {}
qs = []
for y, x in R_.groupby('release_year'):
    m, lo, hi = rate_ci(x.runs_lb >= 100, x.g); s3, l3, h3 = rate_ci(x.runs_lb <= 3, x.g)
    est, mlo, mhi = gboot(x, lambda z, wt: [wmedian(z.runs_lb, wt)], reps=500)
    qs.append({'release_year': y, 'datasets': len(x), 'groups': x.g.nunique(), **{f'p{int(q*100)}': x.runs_lb.quantile(q) for q in (.1, .25, .5, .75, .9)},
               'median_ci_lo': mlo[0], 'median_ci_hi': mhi[0], 'share_ge100': m, 'ge100_lo': lo, 'ge100_hi': hi,
               'share_le3': s3, 'le3_lo': l3, 'le3_hi': h3, 'exact_share': x.exact.mean(),
               'median_exact_only': x[x.exact].runs_exact.median(), 'share_ge100_exact_only': (x[x.exact].runs_exact >= 100).mean()})
qs = pd.DataFrame(qs); w(qs, 'A0-distribution-by-release-year.csv')
CONTRAST = [period_contrast(R_, label='all datasets with runs (claim population)'),
            period_contrast(R_[R_.exact], col='runs_exact', label='exact counts only')]
S['A0_year_2013_vs_2014_16'] = {'2013_n': int((R_.release_year == 2013).sum()), '2014_16_median': float(R_[R_.release_year.between(2014, 2016)].runs_lb.median()),
                                '2013_16_median': float(R_[R_.release_year.between(2013, 2016)].runs_lb.median())}
# submission year (closer to acquisition); release-year bins hide submission years that are not yet public
R_['sub_year'] = R_.submission_date.dt.year
R_['lag_days'] = (R_.release_date - R_.submission_date).dt.days
sy = []
for (lagcap, y), x in pd.concat([R_.assign(lagcap='any lag'), R_[R_.lag_days <= 365].assign(lagcap='released <=365 d after submission')]).groupby(['lagcap', 'sub_year']):
    if y < 2012 or (lagcap == 'released <=365 d after submission' and y > 2024): continue
    m, lo, hi = rate_ci(x.runs_lb >= 100, x.g)
    est, mlo, mhi = gboot(x, lambda z, wt: [wmedian(z.runs_lb, wt)], reps=300)
    sy.append({'window': lagcap, 'submission_year': y, 'datasets': len(x), 'median_runs': x.runs_lb.median(), 'median_lo': mlo[0], 'median_hi': mhi[0],
               'share_ge100': m, 'ge100_lo': lo, 'ge100_hi': hi})
sy = pd.DataFrame(sy); w(sy, 'A5-by-submission-year-fixed-release-window.csv')
# runs vs release lag (are long-embargo datasets larger? if so recent release years miss large datasets)
R_['lag_bin'] = pd.cut(R_.lag_days, [-1e9, 90, 180, 365, 730, 1e9], labels=['<=90 d', '91-180 d', '181-365 d', '1-2 y', '>2 y'])
lagt = []
for (per, lb), x in R_[R_.sub_year.between(2016, 2022)].assign(per='submitted 2016-22').groupby(['per', 'lag_bin'], observed=True):
    m, lo, hi = rate_ci(x.runs_lb >= 100, x.g)
    lagt.append({'population': per, 'lag_bin': lb, 'datasets': len(x), 'median_runs': x.runs_lb.median(), 'share_ge100': m, 'ge100_lo': lo, 'ge100_hi': hi})
lagt = pd.DataFrame(lagt); w(lagt, 'A5-runs-by-release-lag.csv')
print(qs[['release_year', 'datasets', 'p25', 'p50', 'p75', 'p90', 'median_ci_lo', 'median_ci_hi', 'share_ge100', 'share_le3']].round(3).to_string(index=False))
print(sy.round(3).to_string(index=False)); print(lagt.round(3).to_string(index=False))

# ============================================================ A1  splitting: publication components, episodes, group-years
pp = c.execute(r"""
  select project_accession accession, 'doi:' || lower(regexp_replace(trim(publication_doi), '^(https?://(dx\.)?doi\.org/|doi:)', '')) k
    from pride_publications where nullif(trim(publication_doi), '') is not null
  union select project_accession, 'pmid:' || publication_pubmed_id from pride_publications where coalesce(publication_pubmed_id, 0) > 0
  union select accession, 'oa:' || work_id from lit.project_publications where work_id is not null""").df()
pp = pd.concat([pp, d[d.work_id.notna()].assign(k='oa:' + d.work_id)[['accession', 'k']]]).drop_duplicates()
pp = pp[pp.accession.isin(set(d.accession))]
parent = {}
def find(x):
    while parent.setdefault(x, x) != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for k, accs in pp.groupby('k').accession:
    a0 = accs.iloc[0]; find(a0)          # registers single-dataset publications too
    for a in accs.iloc[1:]:
        ra, rb = find(a0), find(a)
        if ra != rb: parent[rb] = ra
R_['pub_comp'] = [find(a) if a in parent else None for a in R_.accession]
R_['has_pub'] = R_.pub_comp.notna()
csize = R_[R_.has_pub].groupby('pub_comp').size()
S['A1_publication_components'] = {'datasets_with_runs_linked': int(R_.has_pub.sum()), 'of': len(R_), 'components': len(csize),
    'components_multi_dataset': int((csize > 1).sum()), 'datasets_in_multi': int(csize[csize > 1].sum()), 'largest_component': int(csize.max()),
    'components_gt20': int((csize > 20).sum()), 'link_keys': 'PRIDE DOI + PubMed id + OpenAlex work ids (lit.project_publications and features work_id), connected components'}
comp = R_[R_.has_pub].groupby('pub_comp').agg(datasets=('accession', 'size'), runs_lb=('runs_lb', 'sum'), release_year=('release_year', 'min'),
                                              g=('g', 'first')).reset_index()
# group episodes: same hashed group, consecutive submissions <= 30 days apart; group-years
ep = R_[R_.pi_group.notna()].sort_values(['pi_group', 'submission_date']).copy()
ep['gap'] = ep.groupby('pi_group').submission_date.diff().dt.days
ep['ep'] = (ep.gap.isna() | (ep.gap > 30)).groupby(ep.pi_group).cumsum()
epi = ep.groupby(['pi_group', 'ep']).agg(datasets=('accession', 'size'), runs_lb=('runs_lb', 'sum'), release_year=('release_year', 'min'), g=('g', 'first')).reset_index()
gy = R_[R_.pi_group.notna()].groupby(['pi_group', 'release_year']).agg(datasets=('accession', 'size'), runs_lb=('runs_lb', 'sum'), g=('g', 'first')).reset_index()
sp = []
for grain, df in {'dataset (all)': R_, 'dataset (publication-linked only)': R_[R_.has_pub], 'publication component': comp,
                  'group submission episode (<=30 d)': epi, 'group-year': gy}.items():
    for y, x in df.groupby('release_year'):
        m, lo, hi = rate_ci(x.runs_lb >= 100, x.g)
        r = {'grain': grain, 'release_year': y, 'units': len(x), 'median_runs': x.runs_lb.median(), 'share_ge100': m, 'ge100_lo': lo, 'ge100_hi': hi}
        r['mean_datasets_per_unit'] = x.datasets.mean() if 'datasets' in x else 1.0
        r['share_units_multi_dataset'] = (x.datasets > 1).mean() if 'datasets' in x else 0.0
        sp.append(r)
    CONTRAST.append(period_contrast(df, label=f'grain: {grain}'))
sp = pd.DataFrame(sp); w(sp, 'A1-splitting-grains-by-year.csv')
print(sp.pivot(index='release_year', columns='grain', values='median_runs').to_string())
print(sp.pivot(index='release_year', columns='grain', values='mean_datasets_per_unit').round(2).to_string())
# share of 2023-25 datasets whose publication component also holds a dataset released in another year
dsplit = []
for y, x in R_.groupby('release_year'):
    xp = x[x.has_pub]; cs = xp.pub_comp.map(xp.groupby('pub_comp').size())
    dsplit.append({'release_year': y, 'datasets': len(x), 'publication_linked': len(xp), 'share_linked': len(xp) / len(x),
                   'share_linked_in_multi_dataset_publication': (cs > 1).mean(),
                   'share_in_multi_dataset_group_episode': float((ep[ep.release_year == y].groupby(['pi_group', 'ep']).accession.transform('size') > 1).mean())})
w(pd.DataFrame(dsplit), 'A1-datasets-per-publication-by-year.csv')

# ============================================================ A2  systematic run-counting checks
st = c.execute("select distinct a accession, stem from vr where rrole = 'run'").df()
FRAC = re.compile(r'(?:^|[_\-.\s])(?:f|fr|frac|fraction|fxn|fx|band|gb|slice|scx|sax|hph|hprp|bprp|hilic|sf)[_\-\s]?0*\d{1,2}(?=$|[_\-.\s])')
st['base'] = [FRAC.sub('', x) for x in st.stem]
st['tok'] = st.base != st.stem
fam = st[st.tok & (st.base.str.strip('_-. ') != '')].groupby(['accession', 'base']).size()
st['fam_n'] = [fam.get((a, b), 0) if t else 0 for a, b, t in zip(st.accession, st.base, st.tok)]
st['in_frac_family'] = st.fam_n >= 3                          # >= 3 runs differing only by a fraction token
st['key'] = np.where(st.in_frac_family, 'F|' + st.base, 'S|' + st.stem)
sm = st.groupby('accession').agg(samples_named=('key', 'nunique'), frac_runs=('in_frac_family', 'sum'), named=('stem', 'size')).reset_index()
R_ = R_.merge(sm, on='accession', how='left')
for col_ in ('samples_named', 'frac_runs', 'named'): R_[col_] = R_[col_].fillna(0)
R_['samples_proxy'] = R_.samples_named + R_.runs_loose + R_.bundles        # fraction families collapsed
R_['has_frac'] = R_.frac_runs > 0
cnt = []
for y, x in R_.groupby('release_year'):
    m, lo, hi = rate_ci(x.samples_proxy >= 100, x.g)
    cnt.append({'release_year': y, 'datasets': len(x), 'share_exact': x.exact.mean(), 'share_with_bundles': (x.bundles > 0).mean(),
                'bundle_runs_share_of_runs_lb': x.bundles.sum() / x.runs_lb.sum(),
                'share_with_wiff_runs': (x.f_wiff > 0).mean(), 'share_with_bruker_d_runs': (x.f_d > 0).mean(), 'share_with_loose_bruker': (x.runs_loose > 0).mean(),
                'share_datasets_with_fraction_family': x.has_frac.mean(), 'share_named_runs_in_fraction_family': x.frac_runs.sum() / max(x.named.sum(), 1),
                'median_runs_lb': x.runs_lb.median(), 'median_samples_proxy': x.samples_proxy.median(),
                'share_ge100_samples_proxy': m, 'ge100sp_lo': lo, 'ge100sp_hi': hi,
                'median_runs_no_fraction_family_datasets': x[~x.has_frac].runs_lb.median(),
                'share_ge100_no_fraction_family_datasets': (x[~x.has_frac].runs_lb >= 100).mean(),
                'median_runs_thermo_raw_exact': x[x.exact & (x.vendor_single == 'Thermo Fisher')].runs_exact.median()})
cnt = pd.DataFrame(cnt); w(cnt, 'A2-counting-checks-by-year.csv'); print(cnt.round(3).to_string(index=False))
CONTRAST += [period_contrast(R_, col='samples_proxy', label='fraction families collapsed (sample proxy)'),
             period_contrast(R_[~R_.has_frac], label='datasets without fraction-token families'),
             period_contrast(R_[R_.exact & (R_.vendor_single == 'Thermo Fisher')], col='runs_exact', label='Thermo, exact counts only'),
             period_contrast(R_[~R_.accession.str.startswith(('PRD', 'PAD'))], label='PXD accessions only')]

# ============================================================ A3  composition
R_['vendor4'] = R_.vendor_single.where(R_.vendor_single.isin(['Thermo Fisher', 'Bruker', 'SCIEX', 'Waters']), 'other/multiple/unknown')
R_['human_s'] = R_.human.astype('object').map({True: 'human', False: 'non-human'}).fillna('species unknown')
R_['dia_s'] = np.where(R_.dia_expanded, 'DIA detected', 'DIA not detected')
R_['tmt_s'] = np.where(R_.tmt_expanded, 'TMT detected', 'TMT not detected')
R_['pub_s'] = np.where(R_.has_pub, 'publication linked', 'no publication linked')
R_['cont_s'] = R_.continent.fillna('country unresolved')
R_['small_s'] = np.where(R_.runs_lb <= 3, '1-3 runs', '>=4 runs')
for fac in ['human_s', 'dia_s', 'tmt_s', 'submission_type', 'vendor4', 'pub_s', 'cont_s', 'small_s']:
    for lev, x in R_.groupby(fac):
        CONTRAST.append({**period_contrast(x, reps=500, label=f'{fac}={lev}'), 'factor': fac})
topsf = R_[R_.has_pub].primary_subfield.value_counts().head(8).index
for lev in topsf:
    CONTRAST.append({**period_contrast(R_[R_.primary_subfield == lev], reps=500, label=f'subfield={lev}'), 'factor': 'primary_subfield'})
LIKE = {'Thermo, DIA & TMT not detected, human, publication linked':
            (R_.vendor4 == 'Thermo Fisher') & ~R_.dia_expanded & ~R_.tmt_expanded & (R_.human_s == 'human') & R_.has_pub,
        'Thermo, DIA & TMT not detected, non-human': (R_.vendor4 == 'Thermo Fisher') & ~R_.dia_expanded & ~R_.tmt_expanded & (R_.human_s == 'non-human'),
        'Thermo, exact, >=4 runs, no fraction-token family': (R_.vendor4 == 'Thermo Fisher') & R_.exact & (R_.runs_lb >= 4) & ~R_.has_frac}
for k, m_ in LIKE.items():
    CONTRAST.append({**period_contrast(R_[m_], reps=500, label=f'like-for-like: {k}'), 'factor': 'like_for_like'})
# multi-way reweighting of each year to the 2013-16 mix
SVARS = ['vendor4', 'dia_s', 'tmt_s', 'human_s', 'submission_type', 'pub_s']
R_['stratum'] = R_[SVARS].astype(str).agg('|'.join, axis=1)


def reweighted(z, wt, target_mask, ref_mask):
    """Weighted median and >=100 share of target rows reweighted to the stratum mix of ref rows."""
    ref = pd.Series(wt[ref_mask]).groupby(z.stratum.to_numpy()[ref_mask]).sum()
    tw = pd.Series(wt[target_mask]).groupby(z.stratum.to_numpy()[target_mask]).sum()
    common = ref.index.intersection(tw.index[tw > 0])
    f = (ref[common] / ref[common].sum()) / (tw[common] / tw[common].sum())
    ww = wt[target_mask] * pd.Series(z.stratum.to_numpy()[target_mask]).map(f).fillna(0).to_numpy()
    v = z.runs_lb.to_numpy()[target_mask]
    return wmedian(v, ww), (ww * (v >= 100)).sum() / ww.sum(), ref[common].sum() / ref.sum()


refm = R_.release_year.between(*EARLY).to_numpy()
rw = []
for y in YEARS:
    tm = (R_.release_year == y).to_numpy()
    med, sh, cov = reweighted(R_, np.ones(len(R_)), tm, refm)
    rw.append({'release_year': y, 'crude_median': R_[tm].runs_lb.median(), 'reweighted_median_2013_16_mix': med,
               'crude_ge100': (R_[tm].runs_lb >= 100).mean(), 'reweighted_ge100_2013_16_mix': sh, 'ref_mix_covered': cov})
rw = pd.DataFrame(rw); w(rw, 'A3-reweighted-to-2013-16-mix-by-year.csv'); print(rw.round(3).to_string(index=False))
zz = R_[R_.release_year.between(*EARLY) | R_.release_year.between(*LATE)].reset_index(drop=True)
em, lm = zz.release_year.between(*EARLY).to_numpy(), zz.release_year.between(*LATE).to_numpy()
def rw_fn(z, wt):
    me = wmedian(z.runs_lb.to_numpy()[em], wt[em]); se = (wt[em] * (z.runs_lb.to_numpy()[em] >= 100)).sum() / wt[em].sum()
    ml, sl, cov = reweighted(z, wt, lm, em)
    return [me, ml, ml / me, se, sl, sl - se, cov]
est, lo, hi = gboot(zz, rw_fn, reps=300)
S['A3_reweighted_late_to_early_mix'] = {'strata_vars': SVARS, 'strata_early': int(zz[em].stratum.nunique()), 'strata_late': int(zz[lm].stratum.nunique()),
    **{k: [float(est[i]), float(lo[i]), float(hi[i])] for i, k in enumerate(['median_early', 'median_late_reweighted', 'median_ratio',
                                                                           'ge100_early', 'ge100_late_reweighted', 'ge100_diff', 'early_mix_covered'])}}
print(S['A3_reweighted_late_to_early_mix'])
mixt = pd.concat([pd.crosstab(R_.release_year, R_[f], normalize='index').add_prefix(f + '=') for f in SVARS + ['cont_s', 'small_s']], axis=1)
w(mixt.reset_index(), 'A3-composition-by-year.csv')

# ============================================================ A4  geography
ctry = R_.country.fillna('Unresolved')
R_['ctry'] = ctry
ce = R_[R_.release_year.between(*EARLY)].ctry.value_counts(); cl = R_[R_.release_year.between(*LATE)].ctry.value_counts()
for k in [k for k in ce.index if ce[k] >= 30 and cl.get(k, 0) >= 30]:
    CONTRAST.append({**period_contrast(R_[R_.ctry == k], reps=500, label=f'country={k}'), 'factor': 'country'})
newly = set(cl.index) - set(ce[ce >= 20].index)
CONTRAST += [{**period_contrast(R_[R_.ctry != 'China'], reps=500, label='excluding China'), 'factor': 'geography'},
             {**period_contrast(R_[~R_.ctry.isin(newly - {'Unresolved'})], reps=500, label='excluding countries with <20 datasets in 2013-16'), 'factor': 'geography'},
             {**period_contrast(R_[R_.ctry != 'Unresolved'], reps=500, label='country resolved only'), 'factor': 'geography'}]
S['A4_newly_active'] = {'countries': len(newly - {'Unresolved'}), 'late_datasets': int(R_[R_.release_year.between(*LATE) & R_.ctry.isin(newly - {'Unresolved'})].shape[0]),
                        'late_total': int(lm.sum()), 'china_share_early': float((R_[R_.release_year.between(*EARLY)].ctry == 'China').mean()),
                        'china_share_late': float((R_[R_.release_year.between(*LATE)].ctry == 'China').mean())}
keep = set(ce[ce >= 20].index) & set(cl[cl >= 20].index)
zz['cs'] = zz.country.fillna('Unresolved').where(zz.country.fillna('Unresolved').isin(keep), 'Other countries')
zz['l2'] = np.log2(zz.runs_lb); zz['ge'] = (zz.runs_lb >= 100).astype(float)
def shift(z, wt, col):
    t = pd.DataFrame({'s': z.cs, 'p': np.where(lm, 'l', 'e'), 'w': wt, 'wy': wt * z[col].to_numpy()}).groupby(['p', 's'])[['w', 'wy']].sum()
    a, b = t.loc['e'], t.loc['l']; idx = a.index.intersection(b.index)
    ma, mb = a.wy / a.w, b.wy / b.w; sa, sb = a.w / a.w.sum(), b.w / b.w.sum()
    tot = (sb * mb).sum() - (sa * ma).sum(); within = (((sa[idx] + sb[idx]) / 2) * (mb[idx] - ma[idx])).sum()
    return [tot, within, tot - within]
geo_ss = {}
for col_ in ('l2', 'ge'):
    est, lo, hi = gboot(zz, lambda z, wt: shift(z, wt, col_), reps=300)
    geo_ss[col_] = {k: [float(est[i]), float(lo[i]), float(hi[i])] for i, k in enumerate(['total', 'within_country', 'between_country_composition'])}
S['A4_shift_share_country'] = {'mean_log2_runs': geo_ss['l2'], 'share_ge100': geo_ss['ge'], 'strata': len(keep) + 1,
                               'note': 'countries with >=20 datasets in both periods; rest pooled; Unresolved kept as its own stratum'}
print(S['A4_shift_share_country'])

# ============================================================ A6  turnover: within-group change, entry cohorts, newcomer gap
gg = R_[R_.pi_group.notna()].copy(); gg['l2'] = np.log2(gg.runs_lb)
ny = gg.groupby('pi_group').release_year.nunique(); gg = gg[gg.pi_group.map(ny) >= 2].reset_index(drop=True)
def slopes(z, wt, excl_first=False):
    m = np.ones(len(z), bool) if not excl_first else (z.release_year > z.group_first_year).to_numpy()
    z = z[m]; wt = wt[m]; gidx = pd.factorize(z.pi_group)[0]
    sw = np.bincount(gidx, wt); yb = np.bincount(gidx, wt * z.l2) / np.where(sw > 0, sw, 1); tb = np.bincount(gidx, wt * z.release_year) / np.where(sw > 0, sw, 1)
    dy, dt = z.l2.to_numpy() - yb[gidx], z.release_year.to_numpy() - tb[gidx]
    within = (wt * dy * dt).sum() / (wt * dt * dt).sum()
    py, pt = np.average(z.l2, weights=wt), np.average(z.release_year, weights=wt)
    pooled = (wt * (z.l2 - py) * (z.release_year - pt)).sum() / (wt * (z.release_year - pt) ** 2).sum()
    return [within, pooled]
gg['g'] = gg.pi_group
tv = {}
for k, ex in (('all datasets of multi-year groups', False), ('excluding each group\'s first-year datasets', True)):
    est, lo, hi = gboot(gg, lambda z, wt: slopes(z, wt, ex), reps=300)
    tv[k] = {'groups': int(gg.pi_group.nunique()), 'datasets': int(len(gg) if not ex else (gg.release_year > gg.group_first_year).sum()),
             'within_group_slope_log2_per_year': [float(est[0]), float(lo[0]), float(hi[0])],
             'pooled_slope_same_datasets': [float(est[1]), float(lo[1]), float(hi[1])],
             'within_fold_per_decade': float(2 ** (10 * est[0])), 'pooled_fold_per_decade': float(2 ** (10 * est[1]))}
al = R_.assign(l2=np.log2(R_.runs_lb))
tv['pooled slope, all datasets with runs'] = float(np.polyfit(al.release_year, al.l2, 1)[0])
S['A6_within_group_slopes'] = tv; print(json.dumps(tv, indent=1))
R_['entry'] = pd.cut(R_.group_first_year, [2012, 2015, 2018, 2021, 2025], labels=['entered 2013-15', 'entered 2016-18', 'entered 2019-21', 'entered 2022-25'])
coh = R_[R_.entry.notna()].groupby(['entry', 'release_year'], observed=True).agg(datasets=('accession', 'size'), groups=('g', 'nunique'),
      median_runs=('runs_lb', 'median'), share_ge100=('runs_lb', lambda s: (s >= 100).mean())).reset_index()
w(coh, 'A6-entry-cohorts-by-year.csv'); print(coh.pivot(index='release_year', columns='entry', values='median_runs').to_string())
# H1.3 re-test with the same windows, dataset level, group bootstrap
both = set(R_[R_.release_year.between(2015, 2018)].pi_group.dropna()) & set(R_[R_.release_year.between(2022, 2025)].pi_group.dropna())
CONTRAST.append({**period_contrast(R_[R_.pi_group.isin(both)], early=(2015, 2018), late=(2022, 2025), label='groups present in 2015-18 and 2022-25 (H1.3 set)'), 'factor': 'turnover'})
CONTRAST.append({**period_contrast(R_, early=(2015, 2018), late=(2022, 2025), label='all datasets, H1.3 windows'), 'factor': 'turnover'})
CONTRAST.append({**period_contrast(R_[R_.release_year > R_.group_first_year], label='returning-group datasets only (not group first year)'), 'factor': 'turnover'})
CONTRAST.append({**period_contrast(R_[R_.release_year == R_.group_first_year], label='group first-year datasets only'), 'factor': 'turnover'})
# newcomer vs returning, like-for-like strata (year x vendor x DIA x human x continent), 2018-2025
nc = R_[R_.pi_group.notna() & (R_.release_year >= 2018)].copy().reset_index(drop=True)
nc['new'] = (nc.release_year == nc.group_first_year).to_numpy()
nc['st2'] = nc.release_year.astype(str) + '|' + nc.vendor4 + '|' + nc.dia_s + '|' + nc.human_s + '|' + nc.cont_s
nc['l2'] = np.log2(nc.runs_lb); nc['ge'] = (nc.runs_lb >= 100).astype(float)
def gap(z, wt):
    t = pd.DataFrame({'s': z.st2, 'n': z.new, 'w': wt, 'a': wt * z.l2.to_numpy(), 'b': wt * z['ge'].to_numpy()}).groupby(['s', 'n'])[['w', 'a', 'b']].sum().unstack('n').dropna()
    t = t[(t[('w', True)] > 0) & (t[('w', False)] > 0)]; wn = t[('w', True)]
    dl = (t[('a', True)] / wn - t[('a', False)] / t[('w', False)]); dg = (t[('b', True)] / wn - t[('b', False)] / t[('w', False)])
    return [(dl * wn).sum() / wn.sum(), (dg * wn).sum() / wn.sum(), wn.sum() / wt[z.new.to_numpy()].sum()]
est, lo, hi = gboot(nc, gap, reps=300)
S['A6_newcomer_gap_like_for_like'] = {'datasets': len(nc), 'newcomer_datasets': int(nc.new.sum()),
    'mean_log2_runs_diff_new_minus_returning': [float(est[0]), float(lo[0]), float(hi[0])], 'fold': float(2 ** est[0]),
    'ge100_diff': [float(est[1]), float(lo[1]), float(hi[1])], 'newcomer_share_in_matched_strata': float(est[2]),
    'crude_mean_log2_diff': float(nc[nc.new].l2.mean() - nc[~nc.new].l2.mean())}
print(S['A6_newcomer_gap_like_for_like'])
nshare = (R_.release_year == R_.group_first_year)[R_.pi_group.notna()].groupby(R_.release_year).mean()
S['A6_share_datasets_from_group_first_year'] = nshare.round(4).to_dict()

# ============================================================ outputs
CT = pd.DataFrame(CONTRAST); w(CT, 'Contrasts-2013-16-vs-2023-25.csv')
print(CT[['label', 'n_early', 'n_late', 'median_early', 'median_late', 'median_ratio', 'median_ratio_lo', 'median_ratio_hi',
          'ge100_early', 'ge100_late', 'ge100_diff', 'ge100_diff_lo', 'ge100_diff_hi']].round(3).to_string(index=False))
S['contrasts'] = CT.set_index('label')[['n_early', 'n_late', 'median_early', 'median_late', 'median_ratio', 'median_ratio_lo', 'median_ratio_hi',
                                         'ge100_early', 'ge100_late', 'ge100_diff', 'ge100_diff_lo', 'ge100_diff_hi']].round(4).to_dict(orient='index')
S['A0_by_year'] = qs.set_index('release_year')[['datasets', 'p25', 'p50', 'p75', 'p90', 'median_ci_lo', 'median_ci_hi', 'share_ge100', 'ge100_lo', 'ge100_hi', 'share_le3']].round(4).to_dict(orient='index')
S['A5_fixed_window'] = sy[sy.window != 'any lag'].set_index('submission_year')[['datasets', 'median_runs', 'median_lo', 'median_hi', 'share_ge100', 'ge100_lo', 'ge100_hi']].round(4).to_dict(orient='index')
S['A5_lag'] = lagt.set_index('lag_bin')[['datasets', 'median_runs', 'share_ge100', 'ge100_lo', 'ge100_hi']].round(4).to_dict(orient='index')
S['A1_grains'] = sp.pivot(index='release_year', columns='grain', values='median_runs').to_dict(orient='index')
S['A2_counting'] = cnt.set_index('release_year')[['share_exact', 'bundle_runs_share_of_runs_lb', 'share_datasets_with_fraction_family',
    'share_named_runs_in_fraction_family', 'median_samples_proxy', 'share_ge100_samples_proxy']].round(4).to_dict(orient='index')
S['A1_datasets_per_publication'] = pd.DataFrame(dsplit).set_index('release_year')[['share_linked', 'share_linked_in_multi_dataset_publication',
                                                                                    'share_in_multi_dataset_group_episode']].round(4).to_dict(orient='index')
(OUT / 'Summary.json').write_text(json.dumps(S, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))

# ============================================================ Verdict.md
C_ = CT.set_index('label')
def f(label, med=True):
    r = C_.loc[label]
    a = (f"median {r.median_early:g} → {r.median_late:g} (ratio {r.median_ratio:.2f}, {r.median_ratio_lo:.2f}–{r.median_ratio_hi:.2f}); " if med else '')
    return (a + f"≥100 runs {100*r.ge100_early:.1f}% → {100*r.ge100_late:.1f}% ({100*r.ge100_diff:+.1f} pp, {100*r.ge100_diff_lo:+.1f} to {100*r.ge100_diff_hi:+.1f}); "
            f"n = {r.n_early:,} / {r.n_late:,}")
def ci3(v, k=1, pct=False):
    m = 100 if pct else 1
    return f"{m*v[0]:.{k}f} ({m*v[1]:.{k}f} to {m*v[2]:.{k}f})"
q = qs.set_index('release_year'); y5 = sy[sy.window != 'any lag'].set_index('submission_year'); cn = cnt.set_index('release_year')
RW = S['A3_reweighted_late_to_early_mix']; TV = S['A6_within_group_slopes']; NG = S['A6_newcomer_gap_like_for_like']; GS = S['A4_shift_share_country']
wi = TV['all datasets of multi-year groups']; frac_within = wi['within_group_slope_log2_per_year'][0] / wi['pooled_slope_same_datasets'][0]
md = f"""# Attack on H1: experiment scale (raw runs per dataset)

Script: `work/stack-and-reuse/attack_runs.py` (about 2 min). Outputs: this folder. Descriptive claims only.
Population: the claim's own population, 35,141 of 36,249 datasets (released 2013–2025) with ≥1 identified run or raw bundle;
13,425 contributor clusters. The run rule was re-implemented and matches the frozen per-dataset table exactly (0 mismatches).
Intervals: 95% group bootstrap over hashed `pi_group` clusters (accession when missing), seed 20261001; single rates use
`common.cluster_boot`, medians and differences use the same cluster resampling. Periods: early = 2013–16, late = 2023–25 (release year).

## Claim tested

Ledger H1.2/H1.3: "median raw runs per dataset fell from 48 (2013) to 16 (2025); share with ≥100 runs 28.9% → 8.7%;
groups seen in both 2015–18 and 2022–25 did not shrink; new groups' datasets are smaller."

Baseline, pooled periods: {f('all datasets with runs (claim population)')}.

## Alternatives, evidence, verdicts

| # | Alternative | Evidence | Verdict |
|---|---|---|---|
| A0 | The 2013 anchor is a small, unusual year | 2013 has {int(q.loc[2013,'datasets'])} datasets: median {q.loc[2013,'p50']:g} (CI {q.loc[2013,'median_ci_lo']:g}–{q.loc[2013,'median_ci_hi']:.0f}). The median is {q.loc[2014,'p50']:g} in 2014 and {q.loc[2015,'p50']:g} in 2015, so 22 of the 32-run fall happens between 2013 and 2015. The pooled 2013–16 median is 30. The whole distribution moves: p90 {q.loc[2013,'p90']:g} → {q.loc[2025,'p90']:g}, p25 {q.loc[2013,'p25']:g} → {q.loc[2025,'p25']:g}; the share with ≤3 runs is {100*q.loc[2013,'share_le3']:.1f}% → {100*q.loc[2025,'share_le3']:.1f}% (`A0-distribution-by-release-year.csv`). | **WEAKENS the numbers, not the direction.** "48 → 16" and "28.9% → 8.7%" overstate the change; the defensible size is about a halving from the 2013–16 level |
| A1 | Splitting: one study deposited as several datasets | Datasets per unit are flat over time: {S['A1_publication_components']['components']:,} publication components (DOI + PubMed + OpenAlex, connected; {S['A1_publication_components']['datasets_with_runs_linked']:,} linked datasets; largest has {S['A1_publication_components']['largest_component']}) hold 1.10–1.30 datasets each in every year. Publication component: {f('grain: publication component')}. Group submission episode (same group, ≤30 days between submissions): {f('grain: group submission episode (<=30 d)')}. Group-year: {f('grain: group-year')}. Among linked datasets, the share in multi-dataset publications is 10% in 2013 and 18–27% in 2014–2025 (`A1-*.csv`). | **SURVIVES.** Every aggregation halves as the dataset grain does |
| A2 | Run-counting artefacts (archives, multi-file runs, Waters/Bruker folders, fractionation) | Hand audit of 30 + 30 random datasets (`A2-hand-audit-60.csv`, file lists in `Audit-file-lists.csv`): every exact count matched the injections in the file list. `.wiff.scan`/`.mtd` companions, `.d.zip`/`.d.7z` and Waters `.raw.zip` were each counted once. Bundles made 3/30 (2013–15) and 4/30 (2023–25) counts lower bounds; explicit fractionation (gel bands, SCX, Frac1–n) appeared in 5/30 and 4/30. Loose Waters/Bruker folder files appear in ≤16 datasets per year. The exact share is 88–97% in every year. Exact counts only: {f('exact counts only')}. Thermo exact only: {f('Thermo, exact counts only')}. Runs that differ only by a fraction token fell from {100*cn.loc[2013,'share_named_runs_in_fraction_family']:.0f}% of named runs (2013) to {100*cn.loc[2025,'share_named_runs_in_fraction_family']:.0f}% (2025). With fraction families collapsed into one unit (a sample proxy): {f('fraction families collapsed (sample proxy)')}. | **SURVIVES.** Less fractionation explains only a small part. Runs are still not samples: TMT datasets (rising from 0.5% to 14–17% of datasets) put up to 18 samples in each run |
| A3 | Composition: species, subfield, DIA, TMT, submission type, vendor, unpublished, PRD/PAD, very small deposits | The late period was reweighted to the early stratum mix (vendor × DIA × TMT × human × submission type × publication linked; {RW['strata_early']} early strata cover {100*RW['early_mix_covered'][0]:.1f}% of the early mix). Median {RW['median_early'][0]:g} → {RW['median_late_reweighted'][0]:g} (ratio {ci3(RW['median_ratio'],2)}); ≥100 {100*RW['ge100_early'][0]:.1f}% → {100*RW['ge100_late_reweighted'][0]:.1f}% (diff {ci3(RW['ge100_diff'],1,True)} pp). Like-for-like (Thermo, DIA and TMT not detected, human, published): {f('like-for-like: Thermo, DIA & TMT not detected, human, publication linked')}. Non-human: {f('like-for-like: Thermo, DIA & TMT not detected, non-human')}. Excluding 1–3-run datasets: {f('small_s=>=4 runs')}. Unpublished: {f('pub_s=no publication linked')}. COMPLETE: {f('submission_type=COMPLETE')}. The fall also holds in the top publication subfields except Spectroscopy (median {C_.loc['subfield=Spectroscopy','median_early']:g} → {C_.loc['subfield=Spectroscopy','median_late']:g}). **Median does not fall in:** TMT-detected (median {C_.loc['tmt_s=TMT detected','median_early']:g} → {C_.loc['tmt_s=TMT detected','median_late']:g}), SCIEX ({C_.loc['vendor4=SCIEX','median_early']:g} → {C_.loc['vendor4=SCIEX','median_late']:g}), Bruker ({C_.loc['vendor4=Bruker','median_early']:g} → {C_.loc['vendor4=Bruker','median_late']:g}; early n=56). The ≥100 share still falls for TMT ({100*C_.loc['tmt_s=TMT detected','ge100_diff']:+.1f} pp) and Bruker ({100*C_.loc['vendor4=Bruker','ge100_diff']:+.1f} pp), but not for SCIEX ({100*C_.loc['vendor4=SCIEX','ge100_diff']:+.1f} pp) or DIA-detected datasets ({f('dia_s=DIA detected', med=False)}). | **SURVIVES.** Mix changes (more TMT and DIA) work against the decline. The decline sits in Thermo deposits without TMT or DIA, the majority of datasets |
| A4 | Driven by a few countries or newly active regions | China's share rose from {100*S['A4_newly_active']['china_share_early']:.1f}% to {100*S['A4_newly_active']['china_share_late']:.1f}%, and its datasets are small (median 16 → 6). Excluding China: {f('excluding China')}. Excluding {S['A4_newly_active']['countries']} countries with <20 early datasets: {f('excluding countries with <20 datasets in 2013-16')}. Shift-share over {GS['strata']} country strata: within-country change {ci3(GS['mean_log2_runs']['within_country'],2)} of a total {ci3(GS['mean_log2_runs']['total'],2)} log2 runs; for ≥100, within {ci3(GS['share_ge100']['within_country'],1,True)} pp of {ci3(GS['share_ge100']['total'],1,True)} pp. The median falls in 14 of the 16 countries with ≥30 datasets in both periods (CI of the ratio below 1 for the US, Germany, Netherlands, Switzerland, Canada, Sweden, Belgium, Norway and China). Australia and Austria rise, and Denmark's ≥100 share does not fall (all with wide CIs). | **SURVIVES.** About 90% of the change is within countries |
| A5 | Embargo censoring: large studies of recent years are still private | Release lag is unrelated to size (submitted 2016–22: median 18–21 runs and ≥100 share 11–12% in every lag bin, `A5-runs-by-release-lag.csv`). By submission year, keeping only datasets released ≤365 days after submission (complete for 2013–2024): median {y5.loc[2013,'median_runs']:g} (2013, n={int(y5.loc[2013,'datasets'])}) → {y5.loc[2016,'median_runs']:g} (2016) → {y5.loc[2024,'median_runs']:g} (2024, CI {y5.loc[2024,'median_lo']:g}–{y5.loc[2024,'median_hi']:g}); ≥100 {100*y5.loc[2013,'share_ge100']:.1f}% → {100*y5.loc[2024,'share_ge100']:.1f}%. | **SURVIVES** |
| A6a | "Continuing groups did not shrink" | Same groups and windows as H1.3 (both 2015–18 and 2022–25), dataset level: {f('groups present in 2015-18 and 2022-25 (H1.3 set)')}. All datasets, same windows: {f('all datasets, H1.3 windows')}. Within-group slope (group fixed effects; {wi['groups']:,} groups with ≥2 release years): {ci3(wi['within_group_slope_log2_per_year'],3)} log2 runs per year (×{wi['within_fold_per_decade']:.2f} per decade). The pooled slope on the same datasets is {ci3(wi['pooled_slope_same_datasets'],3)} (×{wi['pooled_fold_per_decade']:.2f}). So about {100*frac_within:.0f}% of the decline is within groups. Excluding each group's first-year datasets gives the same split. | **WEAKENED.** Continuing groups shrank modestly; "did not shrink" is wrong. The ledger's own interval (0.87–1.00) already excluded growth |
| A6b | "New groups' datasets are smaller" (turnover) | Each entry cohort starts lower and stays roughly flat (median runs in 2022–25: groups entering 2013–15 = 24–26, 2016–18 = 23–24, 2019–21 = 20–21, 2022–25 = 12; `A6-entry-cohorts-by-year.csv`). First-year datasets vs returning groups, same year × vendor × DIA × human × continent stratum (2018–25, {NG['datasets']:,} datasets): ×{NG['fold']:.2f} runs (log2 diff {ci3(NG['mean_log2_runs_diff_new_minus_returning'],2)}); ≥100 {ci3(NG['ge100_diff'],1,True)} pp. About 42–45% of each year's datasets since 2021 come from a group's first year. | **SURVIVES** like-for-like. Caveat: a hashed PI group is a proxy, so a lab that changes contact PI looks like a new group |

## Overall verdict: WEAKENED (direction robust; magnitude and turnover wording need correcting)

The decline in raw runs per dataset is not an artefact of splitting, counting rules, embargo timing, the instrument, method or species mix, or
geography. Each of these attacks leaves the late/early ratio near 0.5. What fails is (i) the anchor: the 2013 numbers (n = 201) roughly double
the size of the change against a 2013–16 baseline, and (ii) the turnover statement. Continuing groups shrank too, by a smaller amount, and about
{100*frac_within:.0f}% of the decline is within groups. Not tested: samples per study. With TMT multiplexing, fewer runs can hold more samples. Larger runs (H1.4) kept bytes per dataset
flat (H1.5).

## Corrected wording

> Raw runs per dataset roughly halved. Median {C_.loc['all datasets with runs (claim population)','median_early']:g} runs in 2013–16
> (n = {C_.loc['all datasets with runs (claim population)','n_early']:,}) vs {C_.loc['all datasets with runs (claim population)','median_late']:g} in 2023–25
> (n = {C_.loc['all datasets with runs (claim population)','n_late']:,}; ratio {C_.loc['all datasets with runs (claim population)','median_ratio']:.2f}, 95% group-bootstrap CI
> {C_.loc['all datasets with runs (claim population)','median_ratio_lo']:.2f}–{C_.loc['all datasets with runs (claim population)','median_ratio_hi']:.2f}).
> The share with ≥100 runs fell from {100*C_.loc['all datasets with runs (claim population)','ge100_early']:.1f}% to {100*C_.loc['all datasets with runs (claim population)','ge100_late']:.1f}%
> ({100*C_.loc['all datasets with runs (claim population)','ge100_diff']:+.1f} pp, {100*C_.loc['all datasets with runs (claim population)','ge100_diff_lo']:+.1f} to {100*C_.loc['all datasets with runs (claim population)','ge100_diff_hi']:+.1f}).
> The small 2013 cohort (n = 201; median 48, 28.9% ≥100) exaggerates the change.
> The decline holds when datasets of one publication, one group-month or one group-year are summed. It also holds in exact counts, with
> fraction-numbered runs collapsed, by submission year with a fixed one-year release window, after reweighting to the 2013–16 vendor × DIA × TMT ×
> species × submission-type × publication mix, and within countries (about 90% of the change).
> The median does not fall in TMT-labelled, SCIEX or Bruker datasets (few early Bruker datasets).
> About {100*frac_within:.0f}% of the decline is within contributor groups (×{wi['within_fold_per_decade']:.2f} per decade). Groups present in
> both 2015–18 and 2022–25 shrank modestly (≥100 share {100*C_.loc['groups present in 2015-18 and 2022-25 (H1.3 set)','ge100_early']:.1f}% →
> {100*C_.loc['groups present in 2015-18 and 2022-25 (H1.3 set)','ge100_late']:.1f}%). The rest comes with contributor turnover: each later
> entry cohort of groups deposits smaller datasets, and a group's first-year datasets have ×{NG['fold']:.2f} the runs of returning groups' in the
> same year, vendor, DIA, species and continent stratum.
> Runs are not samples, and bytes per dataset did not fall (H1.5). This is descriptive: no claim about why.
"""
(OUT / 'Verdict.md').write_text(md)
print('done')
