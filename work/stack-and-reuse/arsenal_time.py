"""H2 — time structure of deposition (new experiment, October 2026).

Descriptive only. Reads the frozen databases read-only; writes only to
outputs/Stack-and-reuse/Arsenal/H2-time/.

1. Data age at deposit: a candidate acquisition date is parsed from each vendor raw file name
   (layer 'vendor_raw' of layer_rules.RULES, non-generated paths). Age at deposit =
   PRIDE submission date - latest (and median) parsed acquisition date of the dataset.
2. Submission -> public release lag and timing relative to the linked paper (OpenAlex date).
3. Acquisition span within a dataset, acquisition weekday, submission calendar.

Parsing rules (conservative; a date in a file name is a label chosen by the lab, usually the
acquisition or queue date, never verified here):
  Unambiguous tokens (tier A), digits not adjacent to further digits:
    - 14-digit timestamps YYYYMMDDhhmmss (valid clock time)
    - 8-digit YYYYMMDD
    - YYYY-MM-DD with one repeated separator (-, _ or .), two-digit month and day
    - day + English month abbreviation + 2- or 4-digit year (e.g. 8Jun14, 28OCT13)
  Order-ambiguous tokens (tier B), resolved per dataset and family, using only files with no
  tier-A date:
    - 6-digit tokens: YYMMDD, DDMMYY, MMDDYY
    - 8-digit tokens not valid as YYYYMMDD: DDMMYYYY, MMDDYYYY
    - DD-MM-YYYY / MM-DD-YYYY with one repeated separator, two-digit fields
    An order is valid for a dataset only if every token of that family in the dataset gives a
    real calendar date in 2005-2025. A family is resolved if exactly one order is valid, or if
    all valid orders give identical dates for every token; otherwise it is left unparsed.
  Every candidate date must be a real calendar date with year 2005-2025.
  A file with more than one distinct candidate date is 'conflicting' and gets no date.
  Separated two-digit-year forms (YY-MM-DD) are never used (PNNL names carry LC column IDs in
  that form).
The submission date is not used for parsing; it is used afterwards only for the plausibility
filter (latest acquisition date > submission + 2 days = impossible -> dataset discarded).
"""
from pathlib import Path
import sys, json
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import R, O, SEED, features, cluster_boot
from layer_rules import RULES
sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Arsenal/H2-time'; OUT.mkdir(parents=True, exist_ok=True)
YMIN, YMAX = 2005, 2025
REPS = 2000
COVER = 0.5          # primary: at least half of a dataset's raw files carry a date
TOL = 2              # days an acquisition date may follow the submission date
S = {}               # Summary.json


def save(df, name):
    df.to_csv(OUT / name, index=False)


# ---------------------------------------------------------------- data -------------------
c = connect(); c.execute("set threads=1; set memory_limit='2GB'")
case = 'case ' + ' '.join(f"when {cond} then '{lay}'" for lay, cond in RULES) + " else 'other' end"
raw = c.execute(rf"""
with f as (select project_accession a, file_name, lower(file_name) n,
  lower(regexp_replace(file_name, '\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b, file_category_value cat
  from pride_files where project_accession in (select accession from study)
   and coalesce(primary_public_location_value,'') not like '%/generated/%')
select a accession, file_name from f where ({case}) = 'vendor_raw'""").df()
works = c.execute("""select pp.accession, w.work_id, w.publication_date paper_date, w.work_type
  from lit.project_publications pp join lit.openalex_works w using (work_id)
  where pp.accession in (select accession from study)""").df()

d = features()
form = pd.read_csv(O / 'Habits/Dataset-form.csv.gz', usecols=['accession', 'form'])
d = d.merge(form, on='accession', how='left', validate='one_to_one')
assert len(d) == 36249
d['grp'] = d.pi_group.fillna('acc:' + d.accession)          # 5 datasets without a group hash
d['sub_year'] = d.submission_date.dt.year
d['dia'] = np.where(d.dia_expanded, 'DIA', 'not DIA')
d['human_s'] = d.human.map({True: 'human', False: 'not human'}).fillna('unknown')
d['contributor'] = np.where(d.group_seq == 1, 'first dataset of group', 'later dataset of group')
d['single_instrument'] = np.where(d.instruments.fillna('').str.contains(r'\|'), None,
                                  d.instruments.fillna('').str.strip())
raw = raw[raw.accession.isin(d.accession)].reset_index(drop=True)
S['raw_files'] = len(raw); S['datasets_with_raw_files'] = int(raw.accession.nunique())

# ---------------------------------------------------------------- parsing ----------------
s = raw.file_name.str.lower()
MON = {m: i + 1 for i, m in enumerate('jan feb mar apr may jun jul aug sep oct nov dec'.split())}


def mkdate(y, m, dd):
    y, m, dd = (pd.to_numeric(v, errors='coerce') for v in (y, m, dd))
    ok = y.between(YMIN, YMAX) & m.between(1, 12) & dd.between(1, 31)
    out = pd.Series(pd.NaT, index=y.index, dtype='datetime64[ns]')
    if ok.any():
        out[ok] = pd.to_datetime(pd.DataFrame({'year': y[ok], 'month': m[ok], 'day': dd[ok]}), errors='coerce')
    return out


def ext(pat):
    m = s.str.extractall(pat)
    m['row'] = m.index.get_level_values(0)
    return m.reset_index(drop=True)


U = []   # unambiguous: row, date, kind
m = ext(r'(?<!\d)(20\d\d)(\d\d)(\d\d)(?:[01]\d|2[0-3])[0-5]\d[0-5]\d(?!\d)')
U.append(pd.DataFrame({'row': m.row, 'date': mkdate(m[0], m[1], m[2]), 'kind': 'YYYYMMDDhhmmss'}))
m8 = ext(r'(?<!\d)(\d{4})(\d\d)(\d\d)(?!\d)')
a8 = mkdate(m8[0], m8[1], m8[2])
U.append(pd.DataFrame({'row': m8.row, 'date': a8, 'kind': 'YYYYMMDD'}))
m = ext(r'(?<!\d)(20\d\d)([-_.])(\d\d)\2(\d\d)(?!\d)')
U.append(pd.DataFrame({'row': m.row, 'date': mkdate(m[0], m[2], m[3]), 'kind': 'YYYY-MM-DD'}))
m = ext(r'(?<![a-z\d])(\d{1,2})(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)(\d{4}|\d{2})(?![a-z\d])')
yy = pd.to_numeric(m[2]); yy = np.where(yy < 100, yy + 2000, yy)
U.append(pd.DataFrame({'row': m.row, 'date': mkdate(pd.Series(yy, index=m.index), m[1].map(MON), m[0]),
                       'kind': 'DMonYY'}))
U = pd.concat(U, ignore_index=True)
tok_counts = U.groupby('kind').agg(tokens=('row', 'size'), valid=('date', lambda x: int(x.notna().sum())))
U = U.dropna(subset=['date']).drop_duplicates(['row', 'date'])

# order-ambiguous families
A = []
m = ext(r'(?<!\d)(\d\d)(\d\d)(\d\d)(?!\d)')
A.append(pd.DataFrame({'row': m.row, 'fam': '6-digit', 'tok': m[0] + m[1] + m[2],
                       'YYMMDD': mkdate(2000 + pd.to_numeric(m[0]), m[1], m[2]),
                       'DDMMYY': mkdate(2000 + pd.to_numeric(m[2]), m[1], m[0]),
                       'MMDDYY': mkdate(2000 + pd.to_numeric(m[2]), m[0], m[1])}))
yl = m8[a8.isna()]
A.append(pd.DataFrame({'row': yl.row, 'fam': '8-digit year last', 'tok': yl[0] + yl[1] + yl[2],
                       'DDMMYYYY': mkdate(yl[1] + yl[2], yl[0].str[2:], yl[0].str[:2]),
                       'MMDDYYYY': mkdate(yl[1] + yl[2], yl[0].str[:2], yl[0].str[2:])}))
m = ext(r'(?<!\d)(\d\d)([-_.])(\d\d)\2(20\d\d)(?!\d)')
A.append(pd.DataFrame({'row': m.row, 'fam': 'DD-MM-YYYY or MM-DD-YYYY', 'tok': m[0] + m[2] + m[3],
                       'DD-MM-YYYY': mkdate(m[3], m[2], m[0]), 'MM-DD-YYYY': mkdate(m[3], m[0], m[2])}))
A = pd.concat(A, ignore_index=True)
ORDERS = ['YYMMDD', 'DDMMYY', 'MMDDYY', 'DDMMYYYY', 'MMDDYYYY', 'DD-MM-YYYY', 'MM-DD-YYYY']
A['accession'] = raw.accession.to_numpy()[A.row]
A = A[~A.row.isin(U.row)]                    # only files without an unambiguous date


def resolve(g):
    cols = [o for o in ORDERS if o in g and g[o].notna().any()]
    valid = [o for o in cols if g[o].notna().all()]
    if not valid:
        return pd.Series({'order': None, 'status': 'no order valid for all tokens'})
    ref = g[valid[0]]
    if all((g[o] == ref).all() for o in valid[1:]):
        return pd.Series({'order': valid[0] if len(valid) == 1 else '+'.join(valid),
                          'status': 'resolved (one order)' if len(valid) == 1 else 'resolved (orders agree)'})
    return pd.Series({'order': None, 'status': 'ambiguous (orders disagree)'})


fam = (A.drop_duplicates(['accession', 'fam', 'tok']).groupby(['accession', 'fam'])
       .apply(resolve, include_groups=False).reset_index())
A = A.merge(fam, on=['accession', 'fam'])
A = A[A.order.notna()].copy()
A['date'] = [r[o.split('+')[0]] for o, r in zip(A.order, A[ORDERS].to_dict('records'))]
A['kind'] = A.order

cand = pd.concat([U[['row', 'date', 'kind']].assign(tier='A'), A[['row', 'date', 'kind']].assign(tier='B')])
cand = cand.drop_duplicates(['row', 'date'])
per = cand.groupby('row').agg(nd=('date', 'nunique'), date=('date', 'first'), kind=('kind', 'first'),
                              tier=('tier', 'first'))
raw['status'] = 'no date token'
raw.loc[per.index[per.nd == 1], 'status'] = 'dated'
raw.loc[per.index[per.nd > 1], 'status'] = 'conflicting dates'
ok = per[per.nd == 1]
raw.loc[ok.index, 'date'] = ok.date; raw.loc[ok.index, 'kind'] = ok.kind; raw.loc[ok.index, 'tier'] = ok.tier
raw['date'] = pd.to_datetime(raw['date'])

S['parse_rules'] = __doc__.split('Parsing rules')[1].split('"""')[0].strip()
fs = raw.status.value_counts().rename_axis('status').reset_index(name='files')
fs['share'] = fs.files / len(raw)
save(fs, 'Parse-file-status.csv')
kinds = raw[raw.status == 'dated'].groupby(['tier', 'kind']).size().rename('files').reset_index()
kinds['datasets'] = [raw[(raw.tier == t) & (raw.kind == k)].accession.nunique() for t, k in zip(kinds.tier, kinds.kind)]
save(kinds, 'Parse-by-pattern.csv')
save(tok_counts.reset_index(), 'Parse-unambiguous-token-validity.csv')
save(fam.groupby(['fam', 'status']).size().rename('dataset_families').reset_index(), 'Parse-order-resolution.csv')

# ---------------------------------------------------------------- dataset level ---------
dated = raw[raw.status == 'dated'].copy()
dated['day'] = (dated.date - pd.Timestamp('1970-01-01')).dt.days
g = dated.groupby('accession')
ds = pd.DataFrame({'n_dated': g.size(), 'latest_d': g.day.max(), 'earliest_d': g.day.min(),
                   'median_d': g.day.median(), 'n_days': g.day.nunique(),
                   'any_tier_b': g.tier.agg(lambda t: (t == 'B').any())})
ds = ds.join(raw.groupby('accession').size().rename('n_raw'), how='right')
ds = ds.join(raw[raw.status == 'conflicting dates'].groupby('accession').size().rename('n_conflict'))
ds = ds.fillna({'n_dated': 0, 'n_conflict': 0}).reset_index()
ds['coverage'] = ds.n_dated / ds.n_raw
d = d.merge(ds, on='accession', how='left')
d['n_raw'] = d.n_raw.fillna(0)
sub_d = (d.submission_date - pd.Timestamp('1970-01-01')).dt.days
d['lead'] = d.latest_d - sub_d                    # acquisition after submission if > 0
d['age_latest'] = sub_d - d.latest_d
d['age_median'] = sub_d - d.median_d
d['span'] = d.latest_d - d.earliest_d
d['latest_date'] = pd.Timestamp('1970-01-01') + pd.to_timedelta(d.latest_d, unit='D')

funnel = [('cohort', len(d)), ('with vendor raw files', int((d.n_raw > 0).sum())),
          ('with >=1 dated raw file', int((d.n_dated > 0).sum())),
          (f'dated share >= {COVER:.0%}', int((d.coverage >= COVER).sum()))]
cov = d.coverage >= COVER
neg = d.lag_days < 0
funnel.append(('... minus release before submission (submission date overwritten)', int((cov & ~neg).sum())))
imp = cov & ~neg & (d.lead > TOL)
funnel.append((f'... minus latest acquisition > submission + {TOL} d (impossible)', int((cov & ~neg & ~imp).sum())))
d['analysed'] = cov & ~neg & ~imp
save(pd.DataFrame(funnel, columns=['step', 'datasets']), 'Funnel.csv')
S['funnel'] = dict(funnel)
S['impossible_datasets'] = {'n': int(imp.sum()), 'of_covered': int((cov & ~neg).sum()),
                            'share': float(imp.sum() / (cov & ~neg).sum()),
                            'among_tier_a_only': float((imp & ~d.any_tier_b.fillna(False).astype(bool)).sum()
                                                       / (cov & ~neg & ~d.any_tier_b.fillna(False).astype(bool)).sum()),
                            'among_any_tier_b': float((imp & d.any_tier_b.fillna(False).astype(bool)).sum()
                                                      / (cov & ~neg & d.any_tier_b.fillna(False).astype(bool)).sum())}
# file-level impossibles among covered datasets
fl = dated.merge(d.loc[cov & ~neg, ['accession', 'submission_date']], on='accession')
fl_imp = (fl.date - fl.submission_date).dt.days > TOL
S['impossible_files'] = {'n': int(fl_imp.sum()), 'of_dated_files_in_covered': len(fl),
                         'by_tier': fl.assign(i=fl_imp).groupby('tier').i.mean().round(4).to_dict()}

# ---------------------------------------------------------------- bootstrap helpers -----
codes_all, groups_all = pd.factorize(d.grp)
d['gcode'] = codes_all
G = len(groups_all)
rng = np.random.default_rng(SEED)
CB = np.empty((REPS, G), dtype=np.float32)
for r in range(REPS):
    CB[r] = np.bincount(rng.integers(0, G, G), minlength=G)


def gsum(sub, v):
    return np.bincount(sub.gcode, weights=v, minlength=G)


def b_ratio(sub, v):
    """Ratio-of-sums interval with whole contributor groups resampled (all cohort groups are
    drawn, so subsets share draws and differences are paired)."""
    return (CB @ gsum(sub, v)) / np.maximum(CB @ gsum(sub, np.ones(len(sub))), 1e-12)


def b_median(sub, col, q=.5):
    x = sub[col].to_numpy(float); o = np.argsort(x); x = x[o]; gc = sub.gcode.to_numpy()[o]
    res = np.empty(REPS)
    for i in range(0, REPS, 100):
        W = CB[i:i + 100][:, gc]; cw = np.cumsum(W, axis=1)
        res[i:i + 100] = x[np.argmax(cw >= q * cw[:, -1:], axis=1)]
    return res


def ci(b):
    return np.quantile(b, .025), np.quantile(b, .975)


def share_ci(sub, flag):
    lo, hi = cluster_boot(flag.astype(float).to_numpy(), sub.grp.to_numpy())
    return flag.mean(), lo, hi


def age_table(sub, by, minn=1, col='age_latest'):
    rows = []
    for k, x in sub.groupby(by, dropna=False):
        if len(x) < minn:
            continue
        med = x[col].median(); mlo, mhi = ci(b_median(x, col))
        r = {by: k, 'datasets': len(x), 'groups': x.grp.nunique(), 'median_days': med,
             'median_lo': mlo, 'median_hi': mhi, 'p25': x[col].quantile(.25), 'p75': x[col].quantile(.75),
             'p90': x[col].quantile(.9)}
        r['share_le_30d'] = (x[col] <= 30).mean(); r['share_gt_730d'] = (x[col] > 730).mean()
        for lab, f in [('share_le_90d', x[col] <= 90), ('share_gt_365d', x[col] > 365)]:
            p, lo, hi = share_ci(x, f); r[lab] = p; r[lab + '_lo'] = lo; r[lab + '_hi'] = hi
        rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- 1. selection ---------
d['region_country'] = d.country.fillna('unknown')
top_c = d.country.value_counts().index[:15]
d['country15'] = np.where(d.country.isin(top_c), d.country, np.where(d.country.isna(), 'unknown', 'other'))
sel = []
for var in ['release_year', 'vendor_single', 'dia', 'human_s', 'family', 'contributor', 'form', 'country15']:
    for k, x in d[d.n_raw > 0].groupby(var):
        p, lo, hi = share_ci(x, x.analysed)
        sel.append({'variable': var, 'level': k, 'datasets_with_raw': len(x), 'analysed': int(x.analysed.sum()),
                    'share': p, 'lo': lo, 'hi': hi})
save(pd.DataFrame(sel), 'Coverage-selection.csv')
S['coverage_overall'] = float(d.analysed.sum() / (d.n_raw > 0).sum())

a = d[d.analysed].copy()
a['age_latest_pos'] = a.age_latest.clip(lower=0)
S['age_overall'] = {'datasets': len(a), 'groups': int(a.grp.nunique()),
                    'median_latest': float(a.age_latest.median()), 'median_latest_ci': list(map(float, ci(b_median(a, 'age_latest')))),
                    'median_median': float(a.age_median.median()), 'p25_latest': float(a.age_latest.quantile(.25)),
                    'p75_latest': float(a.age_latest.quantile(.75)), 'p90_latest': float(a.age_latest.quantile(.9)),
                    'share_gt_365': list(map(float, share_ci(a, a.age_latest > 365))),
                    'share_le_30': list(map(float, share_ci(a, a.age_latest <= 30))),
                    'share_gt_730': float((a.age_latest > 730).mean())}

# reweighting check: analysed subset reweighted to the composition of all raw-file datasets
cells = ['release_year', 'vendor_single', 'dia', 'country15']
base = d[d.n_raw > 0].groupby(cells).size().rename('N')
got = a.groupby(cells).size().rename('n')
w = (base / got).dropna(); aw = a.join(w.rename('w'), on=cells)
o = np.argsort(aw.age_latest.to_numpy()); cw = np.cumsum(aw.w.to_numpy()[o])
S['age_reweighted'] = {'cells': cells, 'datasets_in_cells_covered': float(base[w.index].sum() / base.sum()),
                       'median_latest': float(aw.age_latest.to_numpy()[o][np.searchsorted(cw, cw[-1] / 2)]),
                       'share_gt_365': float(np.average(aw.age_latest > 365, weights=aw.w))}

# ---------------------------------------------------------------- 1. distributions -----
save(age_table(a, 'release_year'), 'Age-by-release-year.csv')
save(age_table(a, 'sub_year'), 'Age-by-submission-year.csv')
save(age_table(a, 'sub_year', col='age_median'), 'Age-median-date-by-submission-year.csv')
for var, fn in [('dia', 'DIA'), ('vendor_single', 'vendor'), ('human_s', 'human'), ('family', 'software'),
                ('contributor', 'contributor'), ('form', 'results-form')]:
    save(age_table(a, var, minn=30), f'Age-by-{fn}.csv')
cc = a.country.value_counts(); big = cc[cc >= 300].index
save(age_table(a[a.country.isin(big)], 'country'), 'Age-by-country.csv')
# like-for-like: DIA vs not and Bruker vs Thermo within 2019-2025 submissions, human only
ll = a[a.sub_year.between(2019, 2025) & (a.human_s == 'human')]
save(age_table(ll, 'dia'), 'Age-by-DIA-human-2019-2025.csv')
save(age_table(ll[ll.vendor_single.isin(['Thermo Fisher', 'Bruker', 'SCIEX'])], 'vendor_single'),
     'Age-by-vendor-human-2019-2025.csv')

# lead histogram (validity: dates should sit before the submission date)
d_cov = d[cov & ~neg]
bins = [-1e9, -3650, -1825, -730, -365, -180, -90, -30, -7, 0, TOL, 30, 365, 1e9]
labs = ['>10y before', '5-10y', '2-5y', '1-2y', '180-365d', '90-180d', '30-90d', '7-30d', '0-7d before',
        f'0-{TOL}d after', f'{TOL}-30d after', '30-365d after', '>1y after']
h = pd.cut(d_cov.lead, bins, labels=labs, right=True).value_counts().reindex(labs).rename('datasets')
save(h.rename_axis('latest_acquisition_vs_submission').reset_index(), 'Latest-acquisition-vs-submission.csv')

# ---------------------------------------------------------------- trend ----------------
def slope_boot(sub, xcol, ycol, within=False):
    x = sub[xcol].to_numpy(float) - 2019; y = sub[ycol].to_numpy(float)
    if within:
        x = x + 2019 - sub.groupby('grp')[xcol].transform('mean').to_numpy()
        y = y - sub.groupby('grp')[ycol].transform('mean').to_numpy()
        sxy, sxx = gsum(sub, x * y), gsum(sub, x * x)
        est = (x * y).sum() / (x * x).sum(); b = (CB @ sxy) / (CB @ sxx)
    else:
        n, sx, sy, sxy, sxx = (gsum(sub, v) for v in (np.ones(len(x)), x, y, x * y, x * x))
        Bn, Bx, By, Bxy, Bxx = (CB @ v for v in (n, sx, sy, sxy, sxx))
        est = np.cov(x, y, bias=True)[0, 1] / x.var(); b = (Bxy - Bx * By / Bn) / (Bxx - Bx ** 2 / Bn)
    lo, hi = ci(b)
    return est, lo, hi


a['log_age'] = np.log1p(a.age_latest_pos)
a['gt365'] = (a.age_latest > 365).astype(float)
tr = []
for lab, sub, xcol in [('all, submission years 2013-2023', a[a.sub_year.between(2013, 2023)], 'sub_year'),
                       ('all, release years 2013-2025', a, 'release_year'),
                       ('not DIA, 2013-2023', a[a.sub_year.between(2013, 2023) & (a.dia == 'not DIA')], 'sub_year'),
                       ('DIA, 2013-2023', a[a.sub_year.between(2013, 2023) & (a.dia == 'DIA')], 'sub_year'),
                       ('Thermo only, 2013-2023', a[a.sub_year.between(2013, 2023) & (a.vendor_single == 'Thermo Fisher')], 'sub_year'),
                       ('human, 2013-2023', a[a.sub_year.between(2013, 2023) & (a.human_s == 'human')], 'sub_year'),
                       ('unambiguous tokens only, 2013-2023', a[a.sub_year.between(2013, 2023) & ~a.any_tier_b.astype(bool)], 'sub_year'),
                       # start-up check: a group's first datasets cannot carry data older than its start
                       ('groups first seen <=2016, submission years 2018-2023',
                        a[a.sub_year.between(2018, 2023) & (a.group_first_year <= 2016)], 'sub_year'),
                       ('excluding each group\'s first two datasets, 2013-2023',
                        a[a.sub_year.between(2013, 2023) & (a.group_seq > 2)], 'sub_year')]:
    for within in (False, True):
        sw = sub[sub.groupby('grp').grp.transform('size') >= 2] if within else sub
        for ycol, unit in [('log_age', 'log(1+days) per year'), ('gt365', 'share >365 d per year')]:
            est, lo, hi = slope_boot(sw, xcol, ycol, within)
            tr.append({'subset': lab, 'within_group': within, 'outcome': unit, 'datasets': len(sw),
                       'groups': sw.grp.nunique(), 'slope': est, 'lo': lo, 'hi': hi})
save(pd.DataFrame(tr), 'Age-trend-slopes.csv')

# like-for-like by instrument model (single-instrument datasets, submissions 2020-2025)
im = a[a.single_instrument.notna() & (a.single_instrument != '') & a.sub_year.between(2020, 2025)]
big_m = im.single_instrument.value_counts(); big_m = big_m[big_m >= 60].index
save(age_table(im[im.single_instrument.isin(big_m)], 'single_instrument'), 'Age-by-instrument-model-2020-2025.csv')
rows = []
for mdl, x in im[im.single_instrument.isin(big_m)].groupby('single_instrument'):
    if (x.dia == 'DIA').sum() >= 25 and (x.dia == 'not DIA').sum() >= 25:
        for k, y in x.groupby('dia'):
            rows.append({'single_instrument': mdl, 'dia': k, 'datasets': len(y), 'groups': y.grp.nunique(),
                         'median_days': y.age_latest.median(), 'share_gt_365d': (y.age_latest > 365).mean()})
save(pd.DataFrame(rows), 'Age-by-DIA-within-instrument-model-2020-2025.csv')
# same-generation comparison: models launched 2016-2017 (approx.), submissions 2020-2025
gen = {'timsTOF Pro': 'Bruker', 'Orbitrap Fusion Lumos': 'Thermo Fisher', 'Q Exactive HF-X': 'Thermo Fisher'}
save(age_table(im[im.single_instrument.isin(gen)], 'single_instrument'), 'Age-same-generation-models-2020-2025.csv')
# country like-for-like: not DIA, Thermo only, submissions 2016-2025
cl = a[(a.dia == 'not DIA') & (a.vendor_single == 'Thermo Fisher') & a.sub_year.between(2016, 2025)]
ccl = cl.country.value_counts(); save(age_table(cl[cl.country.isin(ccl[ccl >= 200].index)], 'country'),
                                      'Age-by-country-notDIA-Thermo-2016-2025.csv')

# sensitivities of the overall age summary
sens = []
def addsens(lab, x, col='age_latest'):
    sens.append({'specification': lab, 'datasets': len(x), 'median_days': x[col].median(),
                 'share_gt_365d': (x[col] > 365).mean(), 'share_le_30d': (x[col] <= 30).mean()})
addsens('primary (dated share >=50%, latest date)', a)
addsens('median acquisition date instead of latest', a, 'age_median')
addsens('unambiguous tokens only', a[~a.any_tier_b.astype(bool)])
addsens('only YYYYMMDD 8-digit tokens', a[a.accession.isin(set(raw[raw.kind == 'YYYYMMDD'].accession))
                                          & ~a.accession.isin(set(raw[(raw.status == 'dated') & (raw.kind != 'YYYYMMDD')].accession))])
any_cov = d[(d.n_dated > 0) & ~neg & ~(d.lead > TOL)]
addsens('any dated raw file (no coverage threshold)', any_cov)
addsens('dated share >= 90%', a[a.coverage >= .9])
addsens('dated share 100%', a[a.coverage >= 1])
addsens('single-instrument datasets', a[a.single_instrument.notna() & (a.single_instrument != '')])
addsens('one dataset per group (first analysed)', a.sort_values('submission_date').drop_duplicates('grp'))
save(pd.DataFrame(sens), 'Age-sensitivity.csv')

# ---------------------------------------------------------------- validity checks ------
# (a) instrument launch: acquisition before the model's commercial launch (approximate, -6 months)
LAUNCH = {'Orbitrap Astral': '2023-06-01', 'Orbitrap Ascend': '2022-06-01', 'Orbitrap Exploris 480': '2019-06-01',
          'Orbitrap Eclipse': '2019-06-01', 'Q Exactive HF-X': '2017-06-01', 'Orbitrap Fusion Lumos': '2016-06-01',
          'Orbitrap Fusion': '2013-06-01', 'timsTOF Pro 2': '2021-06-01', 'timsTOF Pro': '2017-06-01',
          'Orbitrap Exploris 240': '2020-01-01'}
lv = []
for mdl, dt in LAUNCH.items():
    x = d[(d.coverage >= COVER) & ~neg & (d.single_instrument == mdl)]
    early = x.earliest_d < (pd.Timestamp(dt) - pd.Timedelta(days=182) - pd.Timestamp('1970-01-01')).days
    lv.append({'instrument_model': mdl, 'approx_launch': dt, 'datasets': len(x),
               'earliest_acq_before_launch_minus_6m': int(early.sum()), 'share': early.mean() if len(x) else np.nan,
               'latest_acq_before_launch_minus_6m': int((x.latest_d < (pd.Timestamp(dt) - pd.Timedelta(days=182) - pd.Timestamp('1970-01-01')).days).sum())})
save(pd.DataFrame(lv), 'Validity-instrument-launch.csv')

# (b) weekday of acquisition dates (dataset-day level: each distinct date counted once per dataset)
dd = dated[dated.accession.isin(a.accession)].drop_duplicates(['accession', 'date'])
dd = dd.merge(a[['accession', 'country']], on='accession')
dd['weekday'] = dd.date.dt.day_name()
WD = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
wk = []
for lab, x in [('all', dd)] + [(k, dd[dd.country == k]) for k in
                                ['United States', 'Germany', 'China', 'Israel', 'Iran', 'Saudi Arabia', 'Egypt']]:
    v = x.weekday.value_counts().reindex(WD).fillna(0)
    wk.append({'subset': lab, 'dataset_days': int(v.sum()), 'datasets': x.accession.nunique(),
               **{k: v[k] / max(v.sum(), 1) for k in WD}})
save(pd.DataFrame(wk), 'Validity-acquisition-weekday.csv')

# (c) spot check: 30 random analysed datasets + 15 impossible ones, with example names and parsed dates
rs = np.random.RandomState(SEED)
spot = []
for lab, pool, n in [('analysed', a, 30), ('impossible', d[imp], 15)]:
    for acc in pool.accession.sample(min(n, len(pool)), random_state=rs).tolist():
        f = raw[raw.accession == acc]
        ex = pd.concat([f[f.status == 'dated'].sample(min(3, (f.status == 'dated').sum()), random_state=rs),
                        f[f.status != 'dated'].head(1)])
        r = d[d.accession == acc].iloc[0]
        spot.append({'set': lab, 'accession': acc, 'submission_date': r.submission_date.date(),
                     'n_raw': int(r.n_raw), 'n_dated': int(r.n_dated), 'latest': r.latest_date.date(),
                     'age_latest_days': r.age_latest, 'span_days': r.span,
                     'examples': ' || '.join(f'{fn} -> {"" if pd.isna(dt) else dt.date()} [{k}]'
                                            for fn, dt, k in zip(ex.file_name, ex.date, ex.kind.fillna(ex.status)))})
save(pd.DataFrame(spot), 'Validity-spot-check.csv')

# ---------------------------------------------------------------- 2. release lag -------
pl = d[~neg].copy()
S['negative_lag_datasets'] = int(neg.sum())
rel_count = d.release_date.value_counts()
pl['batch_day'] = pl.release_date.map(rel_count) >= 50
S['batch_days'] = {'days_with_>=50_releases': int((rel_count >= 50).sum()),
                   'share_of_datasets': float(pl.batch_day.mean())}


def lag_table(sub, by):
    rows = []
    for k, x in sub.groupby(by):
        med = x.lag_days.median(); lo, hi = ci(b_median(x, 'lag_days'))
        r = {by: k, 'datasets': len(x), 'median_lag_days': med, 'median_lo': lo, 'median_hi': hi,
             'p25': x.lag_days.quantile(.25), 'p75': x.lag_days.quantile(.75), 'p90': x.lag_days.quantile(.9),
             'share_le_7d': (x.lag_days <= 7).mean(), 'share_le_30d': (x.lag_days <= 30).mean(),
             'share_released_on_batch_day': x.batch_day.mean(), 'share_with_linked_paper': x.work_id.notna().mean()}
        for lab, f in [('share_gt_365d', x.lag_days > 365), ('share_gt_730d', x.lag_days > 730)]:
            p, lo, hi = share_ci(x, f); r[lab] = p; r[lab + '_lo'] = lo; r[lab + '_hi'] = hi
        rows.append(r)
    return pd.DataFrame(rows)


save(lag_table(pl, 'release_year'), 'Release-lag-by-release-year.csv')
save(lag_table(pl[pl.sub_year >= 2012], 'sub_year'), 'Release-lag-by-submission-year.csv')
S['release_lag_overall'] = {'datasets': len(pl), 'median': float(pl.lag_days.median()),
                            'p25': float(pl.lag_days.quantile(.25)), 'p75': float(pl.lag_days.quantile(.75)),
                            'share_gt_365': float((pl.lag_days > 365).mean())}

# relation to the linked paper (primary paper = cohort work_id; earliest linked work as sensitivity)
wd = works.dropna(subset=['paper_date']).copy(); wd['paper_date'] = pd.to_datetime(wd.paper_date)
prim = pl.merge(wd[['work_id', 'paper_date', 'work_type']].drop_duplicates('work_id'), on='work_id', how='left')
early = wd.groupby('accession').paper_date.min().rename('paper_date_earliest')
prim = prim.merge(early, on='accession', how='left')
prim['jan1'] = (prim.paper_date.dt.month == 1) & (prim.paper_date.dt.day == 1)
S['paper_link'] = {'datasets_nonneg_lag': len(prim), 'with_primary_paper_date': int(prim.paper_date.notna().sum()),
                   'jan1_paper_dates': int(prim.jan1.sum()),
                   'work_types': prim.work_type.value_counts().to_dict(),
                   'earliest_link_differs_from_primary': int((prim.paper_date_earliest < prim.paper_date).sum()),
                   'paper_date_outside_2005_2026': int(((prim.paper_date.dt.year < 2005) | (prim.paper_date.dt.year > 2026)).sum())}


def timing(x, pcol):
    rel = (x.release_date - x[pcol]).dt.days; sub = (x.submission_date - x[pcol]).dt.days
    return pd.DataFrame({'rel': rel, 'sub': sub})


pt = []
for spec, pcol, filt in [('primary paper', 'paper_date', lambda x: x.paper_date.notna()),
                         ('primary paper, excl. 1 Jan dates', 'paper_date', lambda x: x.paper_date.notna() & ~x.jan1),
                         ('primary paper, journal articles only', 'paper_date', lambda x: x.paper_date.notna() & (x.work_type == 'article')),
                         ('primary paper, preprints only', 'paper_date', lambda x: x.paper_date.notna() & (x.work_type == 'preprint')),
                         ('earliest linked paper', 'paper_date_earliest', lambda x: x.paper_date_earliest.notna())]:
    x = prim[filt(prim)].copy(); t = timing(x, pcol)
    x['rel_before'] = t.rel < -30; x['rel_around'] = t.rel.between(-30, 30); x['rel_after'] = t.rel > 30
    x['sub_after_paper'] = t['sub'] > 30; x['rel_days'] = t.rel; x['sub_to_paper'] = -t['sub']
    for k, y in [('all', x)] + list(x.groupby('release_year')):
        r = {'specification': spec, 'release_year': k, 'datasets': len(y),
             'median_release_minus_paper_days': y.rel_days.median(),
             'median_paper_minus_submission_days': y.sub_to_paper.median()}
        for lab in ['rel_before', 'rel_around', 'rel_after', 'sub_after_paper']:
            p, lo, hi = share_ci(y, y[lab]); r[lab] = p; r[lab + '_lo'] = lo; r[lab + '_hi'] = hi
        pt.append(r)
ren = {'rel_before': 'share_released_>30d_before_paper', 'rel_around': 'share_released_within_30d_of_paper',
       'rel_after': 'share_released_>30d_after_paper', 'sub_after_paper': 'share_submitted_>30d_after_paper'}
pt = pd.DataFrame(pt).rename(columns=lambda c_: next((c_.replace(k, v) for k, v in ren.items() if c_.startswith(k)), c_))
save(pt, 'Paper-timing-by-release-year.csv')
xp = prim[prim.paper_date.notna()]
xr = (xp.release_date - xp.paper_date).dt.days
bt = []
for k, y in xp.assign(r=xr).groupby('batch_day'):
    bt.append({'batch_day': k, 'datasets': len(y), 'median_release_minus_paper': y.r.median(),
               'share_released_>30d_after_paper': (y.r > 30).mean(), 'median_lag_days': y.lag_days.median()})
for (k, yr), y in xp.assign(r=xr).groupby(['batch_day', 'release_year']):
    bt.append({'batch_day': k, 'release_year': yr, 'datasets': len(y), 'median_release_minus_paper': y.r.median(),
               'share_released_>30d_after_paper': (y.r > 30).mean(), 'median_lag_days': y.lag_days.median()})
save(pd.DataFrame(bt), 'Paper-timing-by-batch-release.csv')

# ---------------------------------------------------------------- 3. other time --------
sp = a[a.n_dated >= 2].copy()
sp['raw_bin'] = pd.cut(sp.n_raw, [1, 5, 20, 50, 100, 1e9], labels=['2-5', '6-20', '21-50', '51-100', '>100'])
rows = []
for var in ['sub_year', 'dia', 'raw_bin', 'vendor_single']:
    for k, x in sp.groupby(var, observed=True):
        if len(x) < 30:
            continue
        lo, hi = ci(b_median(x, 'span'))
        rows.append({'variable': var, 'level': k, 'datasets': len(x), 'median_span_days': x.span.median(),
                     'median_lo': lo, 'median_hi': hi, 'p75': x.span.quantile(.75), 'p90': x.span.quantile(.9),
                     'share_one_day': (x.span == 0).mean(), 'share_gt_90d': (x.span > 90).mean(),
                     'share_gt_365d': (x.span > 365).mean(), 'median_distinct_days': x.n_days.median(),
                     'median_raw_files': x.n_raw.median()})
save(pd.DataFrame(rows), 'Acquisition-span.csv')
rows = []
for (k, yb), x in sp.assign(yb=pd.cut(sp.sub_year, [2012, 2016, 2019, 2022, 2025], labels=['2013-16', '2017-19', '2020-22', '2023-25'])
                            ).groupby(['raw_bin', 'yb'], observed=True):
    rows.append({'raw_files': k, 'submission_years': yb, 'datasets': len(x), 'median_span_days': x.span.median(),
                 'share_one_day': (x.span == 0).mean(), 'median_distinct_days': x.n_days.median(),
                 'files_per_distinct_day': (x.n_dated / x.n_days).median()})
save(pd.DataFrame(rows), 'Acquisition-span-by-size-and-period.csv')
S['span_overall'] = {'datasets': len(sp), 'median': float(sp.span.median()), 'p75': float(sp.span.quantile(.75)),
                     'p90': float(sp.span.quantile(.9)), 'share_one_day': float((sp.span == 0).mean()),
                     'share_gt_365': float((sp.span > 365).mean())}

# acquisition -> paper (latest acquisition to primary paper date)
ap = a.merge(prim[['accession', 'paper_date']], on='accession', how='left').dropna(subset=['paper_date'])
ap['acq_to_paper'] = (ap.paper_date - ap.latest_date).dt.days
rows = []
for k, x in [('all', ap)] + list(ap.groupby('release_year')):
    rows.append({'release_year': k, 'datasets': len(x), 'median_acq_to_paper_days': x.acq_to_paper.median(),
                 'p25': x.acq_to_paper.quantile(.25), 'p75': x.acq_to_paper.quantile(.75),
                 'median_age_at_deposit': x.age_latest.median(),
                 'median_paper_minus_submission': (x.paper_date - x.submission_date).dt.days.median()})
save(pd.DataFrame(rows), 'Acquisition-to-paper.csv')

# submission calendar (all datasets, non-negative lag)
cal = pd.DataFrame({'month': pl.submission_date.dt.month, 'weekday': pl.submission_date.dt.day_name()})
save(cal.month.value_counts(normalize=True).sort_index().rename('share').rename_axis('submission_month').reset_index(),
     'Submission-month.csv')
save(cal.weekday.value_counts(normalize=True).reindex(WD).rename('share').rename_axis('submission_weekday').reset_index(),
     'Submission-weekday.csv')

# repeat contributors: age at deposit by position in group sequence (like-for-like within groups)
w2 = a[a.groupby('grp').grp.transform('size') >= 2].copy()
w2['seq_bin'] = pd.cut(w2.group_seq, [0, 1, 3, 10, 1e9], labels=['1st', '2nd-3rd', '4th-10th', '>10th'])
save(age_table(w2, 'seq_bin'), 'Age-by-group-sequence.csv')

# dataset-level table for reuse by later scripts (no file names)
keep = ['accession', 'grp', 'release_year', 'sub_year', 'n_raw', 'n_dated', 'n_conflict', 'coverage', 'any_tier_b',
        'latest_date', 'age_latest', 'age_median', 'span', 'n_days', 'lead', 'lag_days', 'analysed']
d[keep].rename(columns={'grp': 'pi_group'}).to_csv(OUT / 'Dataset-time.csv.gz', index=False)

S['bootstrap'] = (f'{REPS} resamples of whole pi_group clusters, seed {SEED}. Shares within one subset: '
                  'common.cluster_boot. Medians and slopes: same scheme with draw counts over all cohort groups '
                  '(subsets share draws).')
S['caveats'] = ['File-name dates are labels chosen by the lab (acquisition, queue or sample date); not verified.',
                'Submission date can be overwritten on resubmission; datasets released before their recorded '
                'submission are excluded.',
                'Cohort = released 2013-2025: recent submission years miss datasets still private (right truncation).',
                'OpenAlex publication dates may be issue dates rather than online dates.']
(OUT / 'Summary.json').write_text(json.dumps(S, indent=2, default=str))
print(json.dumps(S, indent=2, default=str))
