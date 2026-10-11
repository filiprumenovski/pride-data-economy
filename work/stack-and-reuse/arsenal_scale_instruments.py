"""H1 experiment scale and H3 instrument-model life cycles in deposits (arsenal expansion, October 2026).

Writes only to outputs/Stack-and-reuse/Arsenal/H1-H3-scale-instruments/. Descriptive only.

H1 run rule (vendor_raw layer of layer_rules.RULES, non-generated files only). After stripping any
chain of compression suffixes (.zip .gz .bz2 .7z .tar .rar .xz .tgz), a file is
  run         .raw (Thermo file or zipped Waters folder), .d (zipped Bruker/Agilent folder),
              .wiff, .wiff2, .lcd, .uimf, or a split .d archive (name.d.N.dar);
              runs = distinct stems (name without run extension), so x.raw + x.raw.zip count once
  companion   .wiff.scan / .wiff.mtd / .wiff2.scan (bytes counted, not runs)
  run_loose   loose contents of a Bruker folder deposited without the folder: .baf, .tdf_bin, .yep,
              .fid (one per acquisition; .tdf/.baf_idx are their companions and fall in "other")
  maldi_spot  .t2d (one MALDI spot each; not counted as runs)
  bundle      an archive whose name carries no run extension (e.g. raw_data.7z, sample01.zip);
              its number of runs is unknown, so each bundle counts as >= 1 run (lower bound)
  other       anything else in the layer (indexes, calibration, xml, txt, ...)
Run counts are exact for datasets without bundles ("exact" subset) and a lower bound otherwise.

H3 instrument labels come from pride_instruments. Each PRIDE term is typed with the PSI-MS ontology
already in the repo (work/vendor-patterns/psi-ms.obo): a leaf under 'instrument model'
(MS:1000031) is a model; a non-leaf under it is a series or vendor-level term; anything else is not
an instrument model (technique, software, file format, chromatography system, affinity platform).
Four labels whose accession disagrees with the PRIDE name are harmonised by name (table written).
"""
import json, math
import numpy as np, pandas as pd
from functools import lru_cache
from common import R, O, SEED, cluster_boot
from layer_rules import RULES
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Arsenal/H1-H3-scale-instruments'; OUT.mkdir(parents=True, exist_ok=True)
c = connect(); c.execute("set threads=1; set memory_limit='2GB'")
S = {}  # summary
YEARS = list(range(2013, 2026))
rng_seed = SEED


def w(df, name, **kw):
    df.to_csv(OUT / name, index=kw.pop('index', False), **kw)


def gid(d):
    """Cluster id: hashed contributor group, or the accession itself when missing."""
    return d.pi_group.fillna('acc:' + d.accession)


def share_ci(mask, groups):
    m = np.asarray(mask, float)
    if len(m) == 0: return np.nan, np.nan, np.nan
    lo, hi = cluster_boot(m, np.asarray(groups), reps=2000, seed=SEED)
    return m.mean(), lo, hi


def gboot_stat(d, fn, reps=1000, seed=SEED):
    """Group bootstrap of an arbitrary weighted statistic fn(d, weights); clusters = gid(d)."""
    g = pd.factorize(gid(d))[0]; ng = g.max() + 1; rng = np.random.default_rng(seed)
    out = []
    for _ in range(reps):
        cnt = np.bincount(rng.integers(0, ng, ng), minlength=ng)
        out.append(fn(d, cnt[g].astype(float)))
    out = np.array(out, float)
    return np.nanquantile(out, .025, axis=0), np.nanquantile(out, .975, axis=0)


def wmedian(x, wt):
    o = np.argsort(x); x = np.asarray(x)[o]; wt = np.asarray(wt)[o]; cw = np.cumsum(wt)
    return x[np.searchsorted(cw, cw[-1] / 2)] if cw[-1] > 0 else np.nan


# ============================================================ dataset features
d = pd.read_parquet(O / 'Deposit-features.parquet')
lay = pd.read_parquet(O / 'Layers/Dataset-layers.parquet')
bcols = [x for x in lay.columns if x.startswith('bytes_')]
lay['bytes_all'] = lay[bcols].fillna(0).sum(axis=1)
d = d.merge(lay[['accession', 'bytes_all', 'bytes_vendor_raw']], on='accession', how='left', validate='one_to_one')
d['g'] = gid(d)
S['cohort'] = {'datasets': len(d), 'pi_group_missing': int(d.pi_group.isna().sum()),
               'prefix_counts': d.accession.str[:3].value_counts().to_dict()}

# ============================================================ H1 file-level run classification
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
RUNEXT = "('raw','d','wiff','wiff2','lcd','uimf')"
c.execute(rf"""create temp table vr as
with f as (
  select project_accession a, lower(file_name) n,
         lower(regexp_replace(file_name, '\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
         coalesce(file_category_value, '') cat, coalesce(file_size_bytes, 0) sz
  from pride_files
  where project_accession in (select accession from study)
    and coalesce(primary_public_location_value, '') not like '%/generated/%'),
g as (select *, {case} layer from f),
h as (select a, n, cat, sz, regexp_replace(n, '(\.(zip|gz|bz2|7z|tar|rar|xz|tgz))+$', '') s1
      from g where layer = 'vendor_raw'),
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
roles = c.execute("""select rrole, count(*) files, count(distinct a) datasets, sum(sz)/1e12 tb,
  median(sz)/1e9 median_file_gb from vr group by 1 order by files desc""").df()
w(roles, 'H1-run-rule-file-roles.csv')
# spot-check sample (deterministic) of every role, written for audit
spot = c.execute("""select rrole, a accession, n file_name, cat, round(sz/1e9, 3) gb from (
  select *, row_number() over (partition by rrole order by hash(n || a)) r from vr) where r <= 25 order by rrole, r""").df()
w(spot, 'H1-run-rule-spot-check.csv')
print(roles.to_string(index=False))
for r, x in spot.groupby('rrole'):
    print('--', r); print(x.head(8).to_string(index=False))

runs = c.execute("""select a accession,
  count(distinct stem) filter (where rrole = 'run') runs_named,
  count(*) filter (where rrole = 'run_loose') runs_loose,
  count(*) filter (where rrole = 'bundle') bundles,
  count(*) filter (where rrole = 'maldi_spot') maldi_spots,
  coalesce(sum(sz) filter (where rrole in ('run','companion','run_loose')), 0) bytes_resolved,
  coalesce(sum(sz) filter (where rrole = 'bundle'), 0) bytes_bundle,
  sum(sz) bytes_raw_layer,
  count(*) filter (where rrole = 'run' and ext = 'raw') f_raw, count(*) filter (where rrole = 'run' and ext = 'd') f_d,
  count(*) filter (where rrole = 'run' and ext in ('wiff','wiff2')) f_wiff,
  count(*) filter (where rrole = 'run' and sz = 0) zero_byte_runs
  from vr group by 1""").df()
d = d.merge(runs, on='accession', how='left', validate='one_to_one')
for col in runs.columns[1:]:
    d[col] = d[col].fillna(0)
d['runs_exact'] = d.runs_named + d.runs_loose
d['runs_lb'] = d.runs_exact + d.bundles                      # lower bound when bundles are present
d['has_raw'] = d.bytes_raw_layer > 0
d['status'] = np.select([~d.has_raw, d.runs_lb == 0, d.bundles == 0, d.runs_exact == 0],
                       ['no vendor raw', 'raw-layer files but no runs identified', 'exact', 'bundles only'], 'runs + bundles')
d['per_run'] = np.where((d.status == 'exact') & (d.runs_exact > 0), d.bytes_resolved / d.runs_exact.clip(lower=1), np.nan)
d['dia'] = np.where(d.dia_expanded, 'DIA detected', 'DIA not detected')
S['H1_status'] = d.status.value_counts().to_dict()
S['H1_status_by_year'] = pd.crosstab(d.release_year, d.status).to_dict(orient='index')
st = pd.crosstab(d.release_year, d.status, margins=True); w(st.reset_index(), 'H1-run-count-status-by-year.csv')

# largest run counts: printed for a plausibility check
top = d.sort_values('runs_lb', ascending=False).head(30)[['accession', 'release_year', 'runs_named', 'runs_loose', 'bundles',
       'bytes_raw_layer', 'vendor_single', 'dia', 'instruments', 'primary_subfield']]
w(top, 'H1-largest-run-counts.csv'); print(top.head(12).to_string(index=False))

w(d[['accession', 'release_year', 'pi_group', 'vendor_single', 'dia_expanded', 'status', 'runs_named', 'runs_loose', 'bundles',
   'maldi_spots', 'runs_exact', 'runs_lb', 'bytes_resolved', 'bytes_bundle', 'bytes_raw_layer', 'bytes_all', 'per_run']], 'H1-dataset-scale.csv.gz')
R_ = d[d.has_raw & (d.runs_lb > 0)].copy()  # datasets with >= 1 identified run or raw bundle
EX = R_[R_.status == 'exact'].copy()


def scale_table(df, by):
    rows = []
    for key, x in df.groupby(by, observed=True):
        ex = x[x.status == 'exact']
        r = dict(zip(by, key if isinstance(key, tuple) else (key,)))
        r.update(datasets_with_raw=len(x), groups=x.g.nunique(), exact_datasets=len(ex),
                 median_runs_lb=x.runs_lb.median(), p75_runs_lb=x.runs_lb.quantile(.75), p90_runs_lb=x.runs_lb.quantile(.9),
                 mean_runs_lb=x.runs_lb.mean(), median_runs_exact=ex.runs_exact.median(),
                 median_raw_GB=x.bytes_raw_layer.median() / 1e9, median_total_GB=x.bytes_all.median() / 1e9,
                 median_per_run_GB=ex.per_run.median() / 1e9 if len(ex) else np.nan,
                 p90_per_run_GB=ex.per_run.quantile(.9) / 1e9 if len(ex) else np.nan)
        for t in (100, 1000):
            m, lo, hi = share_ci(x.runs_lb >= t, x.g)
            r[f'share_ge{t}_runs'] = m; r[f'ge{t}_ci_lo'] = lo; r[f'ge{t}_ci_hi'] = hi
        mx, lo, hi = share_ci(ex.runs_exact >= 100, ex.g) if len(ex) else (np.nan,) * 3
        r['share_ge100_exact_only'] = mx
        rows.append(r)
    return pd.DataFrame(rows)


by_year = scale_table(R_, ['release_year']); w(by_year, 'H1-scale-by-year.csv')
print(by_year[['release_year', 'datasets_with_raw', 'median_runs_lb', 'p90_runs_lb', 'share_ge100_runs', 'ge100_ci_lo', 'ge100_ci_hi',
               'share_ge1000_runs', 'median_per_run_GB', 'median_total_GB']].round(3).to_string(index=False))
# CI for the year medians (runs, per-run GB): group bootstrap
med_ci = []
for y, x in R_.groupby('release_year'):
    lo, hi = gboot_stat(x, lambda z, wt: wmedian(z.runs_lb.to_numpy(), wt), reps=400)
    ex = x[x.status == 'exact']
    plo, phi = gboot_stat(ex, lambda z, wt: wmedian(z.per_run.to_numpy(), wt), reps=400)
    med_ci.append({'release_year': y, 'median_runs_lb': x.runs_lb.median(), 'ci_lo': lo, 'ci_hi': hi,
                   'median_per_run_GB': ex.per_run.median() / 1e9, 'per_run_ci_lo': plo / 1e9, 'per_run_ci_hi': phi / 1e9})
med_ci = pd.DataFrame(med_ci); w(med_ci, 'H1-year-medians-with-CI.csv')

R_['period'] = pd.cut(R_.release_year, [2012, 2016, 2019, 2022, 2025], labels=['2013-16', '2017-19', '2020-22', '2023-25'])
EX = R_[R_.status == 'exact'].copy()
by_dia = scale_table(R_, ['release_year', 'dia']); w(by_dia, 'H1-scale-by-year-DIA.csv')
by_vendor = scale_table(R_, ['period', 'vendor_single']); w(by_vendor, 'H1-scale-by-period-vendor.csv')
by_type = scale_table(R_, ['period', 'submission_type']); w(by_type, 'H1-scale-by-period-submission-type.csv')
by_sub = scale_table(R_[R_.release_year >= 2020], ['primary_subfield']).sort_values('datasets_with_raw', ascending=False)
w(by_sub, 'H1-scale-by-subfield-2020-25.csv')

# ---------- like-for-like DIA vs not: same year x vendor x human stratum, standardised to the DIA mix
def std_diff(z, wt, t=100):
    z = z.assign(wt=wt, y=(z.runs_lb >= t).astype(float))
    s = z.groupby(['stratum', 'dia']).apply(lambda q: pd.Series({'w': q.wt.sum(), 'yw': (q.y * q.wt).sum()}), include_groups=False).unstack('dia')
    s = s.dropna()
    s = s[(s[('w', 'DIA detected')] > 0) & (s[('w', 'DIA not detected')] > 0)]
    pd_ = s[('yw', 'DIA detected')] / s[('w', 'DIA detected')]; pn = s[('yw', 'DIA not detected')] / s[('w', 'DIA not detected')]
    wd = s[('w', 'DIA detected')]
    return np.array([(pd_ * wd).sum() / wd.sum(), (pn * wd).sum() / wd.sum(), ((pd_ - pn) * wd).sum() / wd.sum()])


ll = R_[R_.release_year >= 2018].copy()
ll['stratum'] = ll.release_year.astype(str) + '|' + ll.vendor_single + '|' + ll.human.astype(str)
like = {}
for t in (100, 1000):
    est = std_diff(ll, np.ones(len(ll)), t); lo, hi = gboot_stat(ll, lambda z, wt: std_diff(z, wt, t), reps=300)
    like[f'ge{t}'] = {'DIA_rate': est[0], 'not_DIA_rate_std': est[1], 'diff_pp': 100 * est[2], 'diff_ci_pp': [100 * lo[2], 100 * hi[2]]}
like['n'] = {'datasets_2018_25_with_raw': len(ll), 'DIA': int((ll.dia == 'DIA detected').sum()), 'groups': int(ll.g.nunique())}
# per-run size, DIA vs not, within Thermo single-model datasets (model fixed below)
S['H1_DIA_like_for_like'] = like; print(json.dumps(like, indent=1, default=float))

# ---------- byte concentration (all non-generated bytes, every cohort dataset)
def conc(x):
    x = np.sort(np.asarray(x, float))[::-1]; tot = x.sum(); n = len(x)
    k = lambda f: x[:max(1, int(round(f * n)))].sum() / tot
    lor = np.cumsum(np.sort(x)) / tot; gini = 1 - 2 * lor.mean() + 1 / n
    return {'datasets': n, 'PB': tot / 1e15, 'top0.1pct_share': k(.001), 'top1pct_share': k(.01), 'top10pct_share': k(.1),
            'top50pct_share': k(.5), 'gini': gini}


conc_all = conc(d.bytes_all.fillna(0)); conc_raw = conc(d.bytes_raw_layer)
top1 = d.sort_values('bytes_all', ascending=False).head(int(round(.01 * len(d))))
conc_all['top1pct_groups'] = int(top1.g.nunique()); conc_all['top1pct_min_TB'] = top1.bytes_all.min() / 1e12
conc_all['top1pct_DIA_share'] = float(top1.dia_expanded.mean()); conc_all['all_DIA_share'] = float(d.dia_expanded.mean())
conc_all['top1pct_vendor'] = top1.vendor_single.value_counts().to_dict()
conc_all['top1pct_release_year'] = top1.release_year.value_counts().sort_index().to_dict()
conc_all['top1pct_median_runs_lb'] = float(top1.runs_lb.median()); conc_all['top1pct_status'] = top1.status.value_counts().to_dict()
# share of top-1% bytes that are vendor raw
conc_all['top1pct_vendor_raw_byte_share'] = float(top1.bytes_raw_layer.sum() / top1.bytes_all.sum())
cy = pd.DataFrame([{'release_year': y, **conc(x.bytes_all.fillna(0))} for y, x in d.groupby('release_year')])
w(cy, 'H1-byte-concentration-by-year.csv')
# group-level concentration: one cluster per group (does a handful of groups hold the bytes?)
gb = d.groupby('g').bytes_all.sum()
conc_groups = conc(gb.to_numpy())
S['H1_bytes_concentration'] = {'all_bytes': conc_all, 'vendor_raw_bytes': conc_raw, 'by_contributor_group': conc_groups}
w(d.sort_values('bytes_all', ascending=False).head(50)[['accession', 'release_year', 'country', 'vendor_single', 'dia', 'runs_lb',
  'status', 'bytes_all', 'bytes_raw_layer', 'primary_subfield', 'instruments']], 'H1-largest-datasets-by-bytes.csv')
print(json.dumps(S['H1_bytes_concentration'], indent=1, default=float)[:2500])

# ============================================================ H3 instrument labels
txt = (R / 'work/vendor-patterns/psi-ms.obo').read_text()
terms = {}
for chunk in txt.split('[Term]')[1:]:
    L = chunk.split('\n[')[0].strip().splitlines()
    tid = next(x[4:] for x in L if x.startswith('id: '))
    terms[tid] = {'name': next(x[6:] for x in L if x.startswith('name: ')), 'parents': [x[6:].split()[0] for x in L if x.startswith('is_a: ')]}
kids = {}
for t, v in terms.items():
    for p in v['parents']: kids.setdefault(p, []).append(t)


@lru_cache(None)
def anc(t):
    s = set()
    for p in terms.get(t, {}).get('parents', []): s |= {p} | anc(p)
    return frozenset(s)


lab = c.execute("""select i.instrument_accession acc, i.instrument_name nm, count(distinct i.project_accession) datasets
  from pride_instruments i join study s on s.accession = i.project_accession group by 1, 2""").df()
vmap = pd.read_csv(R / 'outputs/Global-proteomics-study/Vendor-patterns/Instrument-vendor-map.csv')
vmap = vmap.drop_duplicates(['instrument_accession', 'instrument_name']).rename(columns={'instrument_accession': 'acc', 'instrument_name': 'nm'})
lab = lab.merge(vmap[['acc', 'nm', 'vendor']], on=['acc', 'nm'], how='left')
lab['ontology_name'] = lab.acc.map(lambda t: terms.get(t, {}).get('name'))
# Labels whose PRIDE name disagrees with the accession: harmonised by the product name shown
NAME_FIX = {('MS:100047', 'Orbitrap Elite'): ('LTQ Orbitrap Elite', 'malformed accession; PSI-MS names MS:1001910 "Orbitrap Elite"'),
            ('MS:1002626', 'Exactive Plus'): ('Exactive Plus', 'accession is "iTRAQ reagent 116" in PSI-MS; name kept'),
            ('MS:1001741', 'LTQ Orbitrap Velos'): ('LTQ Orbitrap Velos', 'accession is a ProteomeDiscoverer term in PSI-MS; name kept'),
            ('MS:1000655', 'QSTAR'): ('QSTAR', 'accession is "QSTAR Elite" in PSI-MS; PRIDE name kept')}
NAME_MODEL_ACC = {'LTQ Orbitrap Elite': 'MS:1001910', 'Exactive Plus': 'MS:1002526', 'LTQ Orbitrap Velos': 'MS:1001742', 'QSTAR': 'MS:1000190'}


def label_type(acc, nm):
    if (acc, nm) in NAME_FIX: acc = NAME_MODEL_ACC[NAME_FIX[(acc, nm)][0]]
    if acc.startswith('PRIDE:'): return 'not an MS instrument model (affinity platform)'
    if acc not in terms: return 'not in PSI-MS'
    if acc == 'MS:1000031' or 'MS:1000031' in anc(acc):
        return 'series or vendor-level term' if kids.get(acc) else 'model'
    return 'not an instrument model (technique, software, format or LC system)'


lab['label_type'] = [label_type(a, n) for a, n in zip(lab.acc, lab.nm)]
lab['model'] = [NAME_FIX.get((a, n), (n,))[0] for a, n in zip(lab.acc, lab.nm)]
lab['harmonisation_note'] = [NAME_FIX.get((a, n), (None, ''))[1] for a, n in zip(lab.acc, lab.nm)]
# vendor for harmonised rows: take the vendor of the model's canonical accession
canon_v = lab.loc[lab.label_type == 'model'].groupby('model').vendor.agg(lambda s: s[s != 'Unknown'].iloc[0] if (s != 'Unknown').any() else 'Unknown')
lab['model_vendor'] = lab.model.map(canon_v).fillna(lab.vendor)
w(lab.sort_values('datasets', ascending=False), 'H3-label-harmonisation.csv')
lt = lab.groupby('label_type').agg(labels=('acc', 'size'), dataset_mentions=('datasets', 'sum')).reset_index()
w(lt, 'H3-label-types.csv'); print(lt.to_string(index=False))

pi = c.execute("select distinct project_accession accession, instrument_accession acc, instrument_name nm from pride_instruments").df()
pi = pi[pi.accession.isin(d.accession)].merge(lab[['acc', 'nm', 'label_type', 'model', 'model_vendor']], on=['acc', 'nm'], how='left')
ds_types = pi.groupby('accession').label_type.agg(lambda s: 'model' if (s == 'model').any() else
                                                  ('series or vendor-level term' if (s == 'series or vendor-level term').any() else 'other only'))
d['label_status'] = d.accession.map(ds_types).fillna('no instrument label')
mod = pi[pi.label_type == 'model'][['accession', 'model', 'model_vendor']].drop_duplicates(['accession', 'model'])
mod = mod.merge(d[['accession', 'release_year', 'submission_date', 'country', 'pi_group', 'g', 'dia_expanded']], on='accession')
mod['credit'] = 1 / mod.groupby('accession').model.transform('size')
mod['sub_year'] = mod.submission_date.dt.year
n_models = mod.groupby('accession').size()
S['H3_label_status'] = d.label_status.value_counts().to_dict()
S['H3_models_per_dataset'] = n_models.value_counts().sort_index().to_dict()
ls_year = pd.crosstab(d.release_year, d.label_status); w(ls_year.reset_index(), 'H3-label-status-by-year.csv')
ser = pi[pi.label_type != 'model'].merge(d[['accession', 'release_year']], on='accession')
ser = pd.crosstab(ser.nm, ser.release_year); ser = ser[ser.sum(axis=1) >= 10]; w(ser.reset_index(), 'H3-non-model-labels-by-year.csv'); print(ser.to_string())

# ---------- per model x year
NY = d.groupby('release_year').size()
NYM = d[d.label_status == 'model'].groupby('release_year').size()
my = mod.groupby(['model', 'release_year']).agg(datasets=('accession', 'size'), credit=('credit', 'sum'), groups=('g', 'nunique')).reset_index()
full = pd.MultiIndex.from_product([sorted(mod.model.unique()), YEARS], names=['model', 'release_year'])
my = my.set_index(['model', 'release_year']).reindex(full, fill_value=0).reset_index()
my['all_datasets'] = my.release_year.map(NY); my['model_labelled_datasets'] = my.release_year.map(NYM)
my['share_all'] = my.datasets / my.all_datasets; my['share_labelled'] = my.datasets / my.model_labelled_datasets
my['credit_share_labelled'] = my.credit / my.model_labelled_datasets
my['vendor'] = my.model.map(mod.drop_duplicates('model').set_index('model').model_vendor)
w(my, 'H3-model-by-year.csv')

# same by submission year (sensitivity: closer to acquisition than public release)
NS = d[d.label_status == 'model'].groupby(d.submission_date.dt.year).size()
ms_ = mod.groupby(['model', 'sub_year']).size().rename('datasets').reset_index()
SY0 = int(mod.sub_year.min())


def lifecycle(x, y_col, n_col, sh_col, min_first=5):
    x = x.sort_values(y_col); yrs = x[y_col].to_numpy(); n = x[n_col].to_numpy(); sh = x[sh_col].to_numpy()
    if n.sum() == 0: return {}
    first_any = yrs[n > 0].min(); first5 = yrs[n >= min_first].min() if (n >= min_first).any() else np.nan
    ip = int(np.argmax(sh)); peak, peak_sh = yrs[ip], sh[ip]
    r = {'datasets': int(n.sum()), 'first_year_any': int(first_any), 'first_year_ge5': first5,
         'peak_year_share': int(peak), 'peak_share': peak_sh, 'peak_year_count': int(yrs[int(np.argmax(n))]),
         'last_year_any': int(yrs[n > 0].max()), 'last_year_ge5': yrs[n >= min_first].max() if (n >= min_first).any() else np.nan,
         'share_2025': sh[yrs == 2025][0] if (yrs == 2025).any() else np.nan}
    r['left_censored'] = bool(first_any == YEARS[0] or (not np.isnan(first5) and first5 == YEARS[0]))
    r['right_censored_peak'] = bool(peak == YEARS[-1])
    base = first5 if not np.isnan(first5) else first_any
    q = yrs[(yrs >= base) & (sh >= .25 * peak_sh)]
    r['year_reach_25pct_peak'] = int(q.min()) if len(q) else np.nan
    r['years_first_to_25pct_peak'] = (q.min() - base) if len(q) and not r['left_censored'] else np.nan
    r['years_first_to_peak'] = (peak - base) if not r['left_censored'] else np.nan
    after = yrs[(yrs > peak) & (sh < .5 * peak_sh)]
    r['years_peak_to_half'] = (after.min() - peak) if len(after) else np.nan
    # age-matched uptake (not dependent on a possibly censored peak): share 0, 1, 2 years after first year with >= 5
    for k in (0, 1, 2):
        hit = sh[yrs == base + k] if not np.isnan(base) else []
        r[f'share_age{k}'] = hit[0] if len(hit) and base + k <= YEARS[-1] else np.nan
    return r


life = []
for m_, x in my.groupby('model'):
    if x.datasets.sum() < 50: continue
    r = {'model': m_, 'vendor': x.vendor.iloc[0], **lifecycle(x, 'release_year', 'datasets', 'share_labelled')}
    xs = ms_[ms_.model == m_].set_index('sub_year')[['datasets']].reindex(range(SY0, 2026), fill_value=0).rename_axis('sub_year').reset_index()
    xs['share'] = xs.datasets / xs.sub_year.map(NS).fillna(np.inf)
    ls = lifecycle(xs.assign(model=m_), 'sub_year', 'datasets', 'share')
    r['sub_first_year_any'] = ls.get('first_year_any'); r['sub_first_year_ge5'] = ls.get('first_year_ge5'); r['sub_peak_year_share'] = ls.get('peak_year_share')
    r['groups_total'] = int(mod[mod.model == m_].g.nunique())
    r['share_country_resolved'] = float(mod[mod.model == m_].country.notna().mean())
    life.append(r)
life = pd.DataFrame(life).sort_values('first_year_ge5'); w(life, 'H3-model-life-cycles.csv')
life['peak_completed'] = life.peak_year_share < YEARS[-1]
w(life, 'H3-model-life-cycles.csv')
print(life[['model', 'share_age0', 'share_age1', 'share_age2', 'peak_completed']].round(3).to_string(index=False))
print(life[['model', 'datasets', 'first_year_any', 'first_year_ge5', 'sub_first_year_ge5', 'peak_year_share', 'peak_share', 'year_reach_25pct_peak',
            'years_first_to_25pct_peak', 'years_peak_to_half', 'last_year_any', 'share_2025']].round(3).to_string(index=False))

# anachronism check: datasets released >= 2 years before a model's first year with >= 5 datasets
early = mod.merge(life[['model', 'first_year_ge5']], on='model')
early = early[early.release_year <= early.first_year_ge5 - 2][['accession', 'model', 'release_year', 'submission_date', 'first_year_ge5']]
w(early.sort_values(['model', 'release_year']), 'H3-early-outliers.csv'); S['H3_early_outlier_rows'] = len(early)
print('early outliers', len(early)); print(early.head(20).to_string(index=False))

# ---------- label vs raw-file format consistency (data quality)
fmt = d[['accession', 'f_raw', 'f_d', 'f_wiff']].copy()
vs = mod.groupby('accession').model_vendor.agg(lambda s: '|'.join(sorted(set(s)))).rename('label_vendors')
fmt = fmt.join(vs, on='accession', how='inner')
one = fmt[~fmt.label_vendors.str.contains(r'\|')]
def fmt_check(r):
    has = {'raw': r.f_raw > 0, 'd': r.f_d > 0, 'wiff': r.f_wiff > 0}
    if not any(has.values()): return 'no named run files'
    exp = {'Thermo Fisher': {'raw'}, 'SCIEX': {'wiff'}, 'Bruker': {'d'}, 'Waters': {'raw'}, 'Agilent': {'d'}}.get(r.label_vendors)
    if exp is None: return 'vendor without single format'
    got = {k for k, v in has.items() if v}
    return 'consistent' if got <= exp else ('mixed' if got & exp else 'inconsistent')
one = one.assign(check=one.apply(fmt_check, axis=1))
fc = pd.crosstab(one.label_vendors, one.check); w(fc.reset_index(), 'H3-label-vs-file-format.csv'); print(fc)
S['H3_label_vs_format'] = fc.to_dict(orient='index')

# ---------- fleet diversity per year
def lcomb(a, b):  # log C(a, b)
    return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)


div = []
for y in YEARS:
    x = mod[mod.release_year == y]; N = x.accession.nunique(); cnt = x.groupby('model').size(); cr = x.groupby('model').credit.sum() / N
    nn = 1000
    rare = sum(1 - math.exp(lcomb(N - k, nn) - lcomb(N, nn)) if N - k >= nn else 1 for k in cnt) if N >= nn else np.nan
    p = cr.to_numpy()
    def eff(z, wt):
        q = z.assign(w=wt * z.credit).groupby('model').w.sum(); q = q / q.sum(); return 1 / (q ** 2).sum()
    lo, hi = gboot_stat(x, eff, reps=300)
    div.append({'release_year': y, 'model_labelled_datasets': N, 'models_ge1': len(cnt), 'models_ge5': int((cnt >= 5).sum()),
                'models_ge1pct': int((cnt / N >= .01).sum()), 'rarefied_models_at_1000': rare,
                'effective_models_inv_simpson': 1 / (p ** 2).sum(), 'eff_ci_lo': lo, 'eff_ci_hi': hi,
                'top1_model': cr.idxmax(), 'top1_credit_share': cr.max(), 'top3_credit_share': cr.nlargest(3).sum(),
                'thermo_credit_share': cr[cr.index.map(lambda m_: canon_v.get(m_, '') == 'Thermo Fisher')].sum(),
                'vendors_ge5': int((x.groupby('model_vendor').size() >= 5).sum()),
                'effective_vendors_inv_simpson': 1 / ((x.groupby('model_vendor').credit.sum() / N) ** 2).sum(),
                'effective_models_within_thermo': (lambda q: 1 / ((q / q.sum()) ** 2).sum())(x[x.model_vendor == 'Thermo Fisher'].groupby('model').credit.sum())})
div = pd.DataFrame(div); w(div, 'H3-fleet-diversity-by-year.csv'); print(div.round(3).to_string(index=False))

# ---------- regional timing of first appearance (continent; deposits, not adoption)
geo = json.loads((R / 'work/stack-and-reuse/ne_50m_admin_0_countries.geojson').read_text())
cont = {}
for f in geo['features']:
    p = f['properties']
    for k in ('ADMIN', 'NAME', 'NAME_LONG'): cont[p[k]] = p['CONTINENT']
d['continent'] = d.country.map(cont)
assert d.country.notna().sum() == d.continent.notna().sum(), 'unmapped country'
mod['continent'] = mod.accession.map(d.set_index('accession').continent)
CONTS = ['Europe', 'North America', 'Asia', 'Oceania', 'South America', 'Africa']
S['H3_continent_datasets'] = d.continent.value_counts(dropna=False).to_dict()
mg = mod[mod.continent.notna()]
NyC = d[(d.label_status == 'model') & d.continent.notna()].groupby(['release_year', 'continent']).size()
Ny = d[(d.label_status == 'model') & d.continent.notna()].groupby('release_year').size()
eligible = life[(~life.left_censored) & (life.datasets >= 100)].model.tolist()
reg = []
for m_ in eligible:
    x = mg[mg.model == m_]; Ky = x.groupby('release_year').accession.nunique()
    gfirst = int(x.release_year.min()) if len(x) else np.nan
    for ct in CONTS:
        xc = x[x.continent == ct]; obs = int(xc.release_year.min()) if len(xc) else np.nan
        # first year by which >= 2 distinct contributor groups in the continent had deposited it
        gg = xc.dropna(subset=['pi_group']).sort_values('release_year').drop_duplicates('pi_group')
        obs2 = int(gg.release_year.iloc[1]) if len(gg) >= 2 else np.nan
        # null: model labels independent of continent within each release year (hypergeometric)
        surv = 1.0; S_ = {}
        for y in YEARS:
            N_, K_, n_ = Ny.get(y, 0), Ky.get(y, 0), NyC.get((y, ct), 0)
            p0 = math.exp(lcomb(N_ - K_, n_) - lcomb(N_, n_)) if N_ - K_ >= n_ and N_ > 0 else 0.0
            surv *= p0; S_[y] = surv
        med_exp = next((y for y in YEARS if S_[y] <= .5), np.nan)
        p_late = S_.get(obs - 1, 1.0) if not np.isnan(obs) else S_[YEARS[-1]]  # P(first >= obs) under null
        p_early = 1 - S_[obs] if not np.isnan(obs) else np.nan                  # P(first <= obs) under null
        reg.append({'model': m_, 'continent': ct, 'global_first_year': gfirst, 'first_year': obs, 'lag_vs_global': obs - gfirst if not np.isnan(obs) else np.nan,
                    'first_year_2groups': obs2, 'expected_median_first_year': med_exp, 'p_null_first_by_2025': 1 - S_[YEARS[-1]],
                    'obs_minus_expected': obs - med_exp if not (np.isnan(obs) or np.isnan(med_exp)) else np.nan,
                    'p_late': p_late, 'p_early': p_early, 'continent_datasets_with_model': len(xc)})
reg = pd.DataFrame(reg); w(reg, 'H3-region-first-appearance.csv')
rs = reg.groupby('continent').agg(models=('model', 'size'), seen=('first_year', lambda s: s.notna().sum()),
                                  median_lag_vs_global=('lag_vs_global', 'median'), mean_lag_vs_global=('lag_vs_global', 'mean'),
                                  median_obs_minus_expected=('obs_minus_expected', 'median'),
                                  later_than_null_p_lt_05=('p_late', lambda s: int((s < .05).sum())),
                                  earlier_than_null_p_lt_05=('p_early', lambda s: int((s < .05).sum()))).reindex(CONTS)
rs['continent_datasets_2013_25'] = rs.index.map(d.continent.value_counts())
w(rs.reset_index(), 'H3-region-summary.csv'); print(rs.round(2).to_string())
S['H3_region_models'] = eligible

# share-based timing: year each continent's share reaches 25% of the model's global peak share
shr = []
for m_ in eligible:
    pk = life.set_index('model').loc[m_, 'peak_share']
    x = mg[mg.model == m_].groupby(['continent', 'release_year']).accession.nunique()
    for ct in CONTS:
        yy = [y for y in YEARS if NyC.get((y, ct), 0) >= 20 and x.get((ct, y), 0) / NyC.get((y, ct), 1) >= .25 * pk]
        shr.append({'model': m_, 'continent': ct, 'year_share_reaches_25pct_global_peak': min(yy) if yy else np.nan})
shr = pd.DataFrame(shr).pivot(index='model', columns='continent', values='year_share_reaches_25pct_global_peak')[CONTS]
w(shr.reset_index(), 'H3-region-share-timing.csv'); print(shr.to_string())

# ============================================================ H1 x H3: per-run size by model
mod1 = mod[mod.credit == 1][['accession', 'model', 'model_vendor']]
EXm = EX[EX.per_run > 0].merge(mod1, on='accession')
pm = EXm.groupby('model').agg(datasets=('accession', 'size'), groups=('g', 'nunique'), median_runs=('runs_exact', 'median'),
                              median_per_run_GB=('per_run', lambda s: s.median() / 1e9), p25=('per_run', lambda s: s.quantile(.25) / 1e9),
                              p75=('per_run', lambda s: s.quantile(.75) / 1e9), share_ge100=('runs_exact', lambda s: (s >= 100).mean()),
                              dia_share=('dia_expanded', 'mean'), first_year=('release_year', 'min')).reset_index()
pm = pm[pm.datasets >= 50].sort_values('median_per_run_GB'); w(pm, 'H1-per-run-size-by-model.csv'); print(pm.round(3).to_string(index=False))
MAJOR = pm[pm.datasets >= 300].model.tolist()
pmy = EXm[EXm.model.isin(MAJOR)].groupby(['model', 'release_year']).agg(datasets=('accession', 'size'),
      median_per_run_GB=('per_run', lambda s: s.median() / 1e9), median_runs=('runs_exact', 'median'),
      share_ge100=('runs_exact', lambda s: (s >= 100).mean())).reset_index()
w(pmy, 'H1-per-run-size-by-model-year.csv')
pmd = EXm[EXm.model.isin(MAJOR)].groupby(['model', 'dia']).agg(datasets=('accession', 'size'),
      median_per_run_GB=('per_run', lambda s: s.median() / 1e9), median_runs=('runs_exact', 'median')).reset_index()
w(pmd, 'H1-per-run-size-by-model-DIA.csv'); print(pmd.round(3).to_string(index=False))

# shift-share of mean log10 per-run bytes, 2014-16 vs 2023-25 (single-model exact datasets)
EXm['lp'] = np.log10(EXm.per_run.clip(lower=1))
EXm['per'] = np.where(EXm.release_year.between(2014, 2016), 'a', np.where(EXm.release_year.between(2023, 2025), 'b', None))


def shift_share(z, wt):
    z = z.assign(wt=wt); z = z[z.per.notna() & (z.wt > 0)]
    t = z.groupby(['per', 'model']).apply(lambda q: pd.Series({'w': q.wt.sum(), 'm': np.average(q.lp, weights=q.wt)}), include_groups=False)
    a, b = t.loc['a'], t.loc['b']; a = a.assign(s=a.w / a.w.sum()); b = b.assign(s=b.w / b.w.sum())
    tot = (b.s * b.m).sum() - (a.s * a.m).sum()
    both = a.index.intersection(b.index)
    within = (((a.s[both] + b.s[both]) / 2) * (b.m[both] - a.m[both])).sum()
    return np.array([tot, within, tot - within, a.s[both].sum(), b.s[both].sum()])


ss = shift_share(EXm, np.ones(len(EXm))); lo, hi = gboot_stat(EXm, shift_share, reps=300)
S['H1_per_run_shift_share_log10'] = {'total_change_log10': ss[0], 'within_model': ss[1], 'between_model_and_new_models': ss[2],
    'ci_total': [lo[0], hi[0]], 'ci_within': [lo[1], hi[1]], 'ci_between': [lo[2], hi[2]],
    'fold_total': 10 ** ss[0], 'fold_within': 10 ** ss[1], 'fold_between': 10 ** ss[2],
    'early_share_in_models_present_both': ss[3], 'late_share_in_models_present_both': ss[4],
    'n_datasets': {'2014-16': int((EXm.per == 'a').sum()), '2023-25': int((EXm.per == 'b').sum())}}
print(json.dumps(S['H1_per_run_shift_share_log10'], indent=1, default=float))

# ============================================================ turnover: returning groups' own scale
gp = R_[R_.pi_group.notna()].copy()
gp['win'] = np.where(gp.release_year.between(2015, 2018), 'early', np.where(gp.release_year.between(2022, 2025), 'late', None))
gm = gp[gp.win.notna()].groupby(['pi_group', 'win']).agg(med=('runs_lb', 'median'), n=('accession', 'size'),
                                                        ge100=('runs_lb', lambda s: (s >= 100).mean())).unstack('win').dropna()
ratio = np.log2(gm[('med', 'late')].clip(lower=1) / gm[('med', 'early')].clip(lower=1))
rng = np.random.default_rng(SEED); bs = [np.median(ratio.to_numpy()[rng.integers(0, len(ratio), len(ratio))]) for _ in range(2000)]
dd = gm[('ge100', 'late')] - gm[('ge100', 'early')]
bs2 = [dd.to_numpy()[rng.integers(0, len(dd), len(dd))].mean() for _ in range(2000)]
# compare with all datasets in the same windows
allw = gp[gp.win.notna()].groupby('win').runs_lb.median()
newg = R_[(R_.release_year.between(2022, 2025)) & (R_.group_first_year >= 2022)]
S['H1_returning_groups'] = {'groups_in_both_windows': len(gm), 'median_of_group_log2_ratio_late_vs_early': float(np.median(ratio)),
    'ci': [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], 'fold': float(2 ** np.median(ratio)),
    'share_groups_with_higher_median': float((ratio > 0).mean()), 'share_groups_lower': float((ratio < 0).mean()),
    'mean_change_share_ge100_pp': float(100 * dd.mean()), 'ci_pp': [float(100 * np.quantile(bs2, .025)), float(100 * np.quantile(bs2, .975))],
    'all_datasets_median_runs_2015_18': float(allw.get('early')), 'all_datasets_median_runs_2022_25': float(allw.get('late')),
    'returning_groups_dataset_median_early': float(gp[(gp.win == 'early') & gp.pi_group.isin(gm.index)].runs_lb.median()),
    'returning_groups_dataset_median_late': float(gp[(gp.win == 'late') & gp.pi_group.isin(gm.index)].runs_lb.median()),
    'new_groups_2022_25_dataset_median': float(newg.runs_lb.median()), 'new_groups_2022_25_datasets': len(newg)}
print(json.dumps(S['H1_returning_groups'], indent=1, default=float))

# subfield/human-standardised share >= 100 runs (2019 mix held fixed)
mix = R_.assign(sf=R_.primary_subfield.fillna('none') + '|' + R_.human.astype(str))
ref = mix[mix.release_year == 2019].sf.value_counts(normalize=True)
stdz = []
for y, x in mix.groupby('release_year'):
    r_ = x.groupby('sf').runs_lb.apply(lambda s: (s >= 100).mean())
    common = r_.index.intersection(ref.index)
    stdz.append({'release_year': y, 'crude_ge100': (x.runs_lb >= 100).mean(),
                 'std_ge100_2019_subfield_human_mix': (r_[common] * ref[common]).sum() / ref[common].sum(),
                 'ref_mix_covered': ref[common].sum()})
stdz = pd.DataFrame(stdz); w(stdz, 'H1-ge100-standardised-by-year.csv'); print(stdz.round(3).to_string(index=False))

# newcomer vs returning group datasets per year (group's first release year == this year)
R_['newcomer'] = np.where(R_.pi_group.isna(), 'no group', np.where(R_.group_first_year == R_.release_year, 'group in first year', 'returning group'))
nc = scale_table(R_, ['release_year', 'newcomer']); w(nc, 'H1-scale-by-year-newcomer.csv')
print(nc.pivot(index='release_year', columns='newcomer', values='median_runs_lb').to_string())
print(nc.pivot(index='release_year', columns='newcomer', values='share_ge100_runs').round(3).to_string())

# publication grain: datasets linked to the same OpenAlex work are summed (runs lower bound, bytes)
pub = R_[R_.work_id.notna()].groupby('work_id').agg(datasets=('accession', 'size'), runs_lb=('runs_lb', 'sum'),
        bytes_raw=('bytes_raw_layer', 'sum'), first_release=('release_year', 'min'), g=('g', 'first')).reset_index()
pr = []
for y, x in pub.groupby('first_release'):
    m, lo, hi = share_ci(x.runs_lb >= 100, x.g)
    pr.append({'first_release_year': y, 'publications': len(x), 'multi_dataset_share': (x.datasets > 1).mean(),
               'median_runs_lb': x.runs_lb.median(), 'share_ge100': m, 'ci_lo': lo, 'ci_hi': hi, 'median_raw_GB': x.bytes_raw.median() / 1e9})
pr = pd.DataFrame(pr); w(pr, 'H1-scale-by-year-publication-grain.csv'); print(pr.round(3).to_string(index=False))
S['H1_publication_grain'] = {'datasets_with_work_id': int(R_.work_id.notna().sum()), 'of_datasets_with_runs': len(R_), 'publications': len(pub)}

# ---------- summary
S['H1_year'] = by_year.set_index('release_year')[['datasets_with_raw', 'exact_datasets', 'median_runs_lb', 'median_runs_exact', 'p90_runs_lb',
    'share_ge100_runs', 'ge100_ci_lo', 'ge100_ci_hi', 'share_ge1000_runs', 'ge1000_ci_lo', 'ge1000_ci_hi', 'median_per_run_GB',
    'median_total_GB']].round(4).to_dict(orient='index')
S['H1_overall'] = {'datasets_with_raw': len(R_), 'exact': len(EX), 'share_ge100': float((R_.runs_lb >= 100).mean()),
                   'share_ge1000': float((R_.runs_lb >= 1000).mean()), 'median_runs_lb': float(R_.runs_lb.median()),
                   'zero_byte_named_runs': int(R_.zero_byte_runs.sum()), 'datasets_with_zero_byte_runs': int((R_.zero_byte_runs > 0).sum())}
S['H3_fleet'] = div.set_index('release_year')[['model_labelled_datasets', 'models_ge1', 'models_ge5', 'rarefied_models_at_1000',
                                               'effective_models_inv_simpson', 'eff_ci_lo', 'eff_ci_hi', 'top1_model', 'top1_credit_share', 'thermo_credit_share']].to_dict(orient='index')
S['notes'] = ['Descriptive only. Deposited instrument labels are not ownership; first appearance in deposits is not purchase or adoption.',
              'Run counts are lower bounds for datasets with raw bundles; exact-only figures are given alongside.',
              'Bytes are as deposited (compressed archives count at their compressed size).',
              'Group-bootstrap CIs resample hashed pi_group clusters (datasets without a group are their own cluster); seed 20261001.']
(OUT / 'Summary.json').write_text(json.dumps(S, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
print('done')
