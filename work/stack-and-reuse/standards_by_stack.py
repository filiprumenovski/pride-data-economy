"""Standard-result listing (corrected nongenerated_expanded) by software family, vendor and
country, and an exact decomposition of the post-2019 archive trend.

Descriptive. Intervals resample hashed PI-contact groups (deposits without a PI group are
their own cluster). Publication-deduplicated and equal-group versions are sensitivities.
"""
import warnings; warnings.filterwarnings('ignore', 'This pattern')
import numpy as np, pandas as pd
from common import O, SEED, features, FAMILY_ORDER, cluster_boot

OUT = O / 'Standards-by-stack'; OUT.mkdir(exist_ok=True)
d = features()
d['cluster'] = d.pi_group.fillna('solo:' + d.accession)
d['complete'] = d.submission_type.eq('COMPLETE')
d['vendor_single'] = d.vendor_single.fillna('Unknown')
d['acquisition'] = np.where(d.dia_expanded, 'DIA', 'DDA/other')
# Publication grain: earliest release per linked work, unlinked deposits kept as-is.
pub = d.sort_values(['release_date', 'accession'])
pub = pd.concat([pub[pub.work_id.notna()].drop_duplicates('work_id'), pub[pub.work_id.isna()]])

# 0. outcome anatomy: the standard flag is (almost) PRIDE's COMPLETE tier ---------------
pd.crosstab(d.submission_type, d.standards, margins=True).to_csv(OUT / 'Submission-type-by-standards.csv')

# 1. rates by stratum ---------------------------------------------------------------------
def rate_table(x, by, label, outcomes=('standards', 'complete', 'fmt_mzml', 'fmt_peaklist')):
    rows = []
    for k, g in x.groupby(by):
        row = {'stratifier': label, 'level': k, 'deposits': len(g), 'groups': g.pi_group.nunique()}
        for o in outcomes:
            v = g[o].fillna(False).astype(float).to_numpy()
            row[o] = v.mean()
            if o == 'standards':
                row['standards_ci_low'], row['standards_ci_high'] = cluster_boot(v, g.cluster.to_numpy(), reps=1000)
                w = 1 / g.groupby('cluster').accession.transform('size')
                row['standards_equal_group'] = np.average(v, weights=w)
        rows.append(row)
    return rows

rec = d[d.release_year.between(2019, 2025)]
recp = pub[pub.release_year.between(2019, 2025)]
rows = (rate_table(rec, 'family', 'software named') + rate_table(rec, 'vendor_single', 'vendor (single)')
        + rate_table(rec, 'acquisition', 'acquisition')
        + rate_table(rec[rec.dia_expanded], 'vendor_single', 'vendor within DIA')
        + rate_table(rec[~rec.dia_expanded], 'vendor_single', 'vendor within DDA/other'))
ct = rec.country.value_counts(); rows += rate_table(rec[rec.country.isin(ct[ct >= 200].index)], 'country', 'country (>=200 deposits)')
rows += [dict(r, stratifier=r['stratifier'] + ' [publication grain]') for r in rate_table(recp, 'family', 'software named')]
rt = pd.DataFrame(rows); rt.to_csv(OUT / 'Rates-2019-2025.csv', index=False)
spread = rt[~rt.stratifier.str.contains('publication')].query('deposits>=100').groupby('stratifier').standards.agg(['min', 'max', 'count'])
spread['range_pp'] = 100 * (spread['max'] - spread['min']); spread.to_csv(OUT / 'Rate-spread-by-stratifier.csv')

# 2. annual trend by acquisition, with counterfactual fixed DIA share --------------------
ann = []
for y, x in d.groupby('release_year'):
    row = {'release_year': y, 'deposits': len(x), 'dia_share': x.dia_expanded.mean(), 'standards': x.standards.mean()}
    for a, g in x.groupby('acquisition'):
        v = g.standards.astype(float).to_numpy(); row[f'standards_{a}'] = v.mean(); row[f'n_{a}'] = len(g)
        row[f'ci_low_{a}'], row[f'ci_high_{a}'] = cluster_boot(v, g.cluster.to_numpy(), reps=1000) if len(g) > 20 else (np.nan, np.nan)
    ann.append(row)
ann = pd.DataFrame(ann)
base = d[d.release_year.between(2019, 2021)].dia_expanded.mean()
ann['counterfactual_fixed_dia_share_2019_2021'] = ann['standards_DDA/other'] * (1 - base) + ann['standards_DIA'].fillna(0) * base
ann.to_csv(OUT / 'Annual-standards-by-acquisition.csv', index=False)

# 3. exact shift-share decomposition between windows, by software family ------------------
def decompose(x0, x1, cat):
    w0 = x0[cat].value_counts(normalize=True); w1 = x1[cat].value_counts(normalize=True)
    p0 = x0.groupby(cat).standards.mean(); p1 = x1.groupby(cat).standards.mean()
    # Composition terms are centred on the mean of the two overall rates, so a growing category
    # with a below-average rate counts as negative; the centring leaves every total unchanged.
    cbar = (x0.standards.mean() + x1.standards.mean()) / 2
    rows = []
    for k in sorted(set(w0.index) | set(w1.index)):
        a0, a1 = w0.get(k, 0), w1.get(k, 0); r0, r1 = p0.get(k, np.nan), p1.get(k, np.nan)
        r0f = r0 if a0 else r1; r1f = r1 if a1 else r0
        rows.append({'category': k, 'share_0': a0, 'share_1': a1, 'rate_0': r0, 'rate_1': r1,
                     'within': (a0 + a1) / 2 * (r1f - r0f), 'composition': (a1 - a0) * ((r0f + r1f) / 2 - cbar)})
    t = pd.DataFrame(rows)
    total = x1.standards.mean() - x0.standards.mean()
    assert abs(t.within.sum() + t.composition.sum() - total) < 1e-12
    return t, total

dec = []
for cat in ['family', 'acquisition', 'vendor_single']:
    for grain, src in [('deposit', d), ('publication', pub)]:
        x0 = src[src.release_year.between(2019, 2021)]; x1 = src[src.release_year.between(2023, 2025)]
        t, total = decompose(x0, x1, cat)
        t.insert(0, 'grain', grain); t.insert(0, 'partition', cat); t['total_change'] = total; dec.append(t)
dec = pd.concat(dec); dec.to_csv(OUT / 'Decomposition-2019-2021-vs-2023-2025.csv', index=False)
summ = dec.groupby(['partition', 'grain']).agg(total=('total_change', 'first'), within=('within', 'sum'), composition=('composition', 'sum'))
summ.to_csv(OUT / 'Decomposition-summary.csv')

# group-bootstrap interval for the family-composition term (deposit grain)
rng = np.random.default_rng(SEED)
x0 = d[d.release_year.between(2019, 2021)]; x1 = d[d.release_year.between(2023, 2025)]
clusters = d.cluster.unique(); idx = d.groupby('cluster').indices; boots = []
for _ in range(300):
    take = rng.choice(clusters, len(clusters))
    s = d.iloc[np.concatenate([idx[c] for c in take])]
    t, total = decompose(s[s.release_year.between(2019, 2021)], s[s.release_year.between(2023, 2025)], 'family')
    boots.append((total, t.within.sum(), t.composition.sum()))
b = np.array(boots)
pd.DataFrame({'term': ['total', 'within_family', 'family_composition'], 'estimate': summ.loc[('family', 'deposit')].to_numpy(),
              'ci_low': np.quantile(b, .025, axis=0), 'ci_high': np.quantile(b, .975, axis=0)}).to_csv(OUT / 'Decomposition-family-ci.csv', index=False)
print(rt[rt.stratifier.isin(['software named', 'acquisition', 'vendor (single)'])][['stratifier', 'level', 'deposits', 'groups', 'standards', 'standards_ci_low', 'standards_ci_high', 'standards_equal_group', 'complete']].round(3).to_string(index=False))
print(spread.round(3)); print(ann.round(3).to_string(index=False)); print(summ.round(4)); print(pd.read_csv(OUT / 'Decomposition-family-ci.csv').round(4))
