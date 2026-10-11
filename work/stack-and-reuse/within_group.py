"""Repeat contributors: does the same group include identification files for its DIA and DDA work alike?

Like-for-like designs within hashed PI-contact groups. These describe co-occurrence within
the same contributor; they are not causal effects of switching software, and a first
recorded tool is an observation event, not an adoption date.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from common import O, SEED, features, cluster_boot

OUT = O / 'Within-group'; OUT.mkdir(exist_ok=True)
d = features()
d = d[d.pi_group.notna()].sort_values(['pi_group', 'release_date', 'accession']).reset_index(drop=True)
d['std'] = d.standards.astype(int); d['dia'] = d.dia_expanded.astype(int); d['native'] = d.dia_native_tool.astype(int)
d['native_any'] = (d.dia_native_tool | d.tool_diann.fillna(False) | d.tool_spectronaut.fillna(False)).astype(int)

# A. fixed-effect linear probability models ---------------------------------------------
def fe(data, x, fe_col, label):
    v = data.groupby(fe_col)[x].agg(['min', 'max']); keep = v[(v['min'] == 0) & (v['max'] == 1)].index
    s = data[data[fe_col].isin(keep)]
    m = smf.ols(f'std ~ {x} + C(release_year) + C({fe_col})', s).fit(
        cov_type='cluster', cov_kwds={'groups': s.pi_group.astype('category').cat.codes})
    lo, hi = m.conf_int().loc[x]
    return {'design': label, 'exposure': x, 'deposits': len(s), 'groups': s.pi_group.nunique(), 'units': len(keep),
            'raw_rate_exposed': s.loc[s[x] == 1, 'std'].mean(), 'raw_rate_unexposed': s.loc[s[x] == 0, 'std'].mean(),
            'estimate': m.params[x], 'ci_low': lo, 'ci_high': hi}

d['group_instrument'] = d.pi_group + '|' + d.instruments
rows = [fe(d, 'dia', 'pi_group', 'within group'), fe(d, 'dia', 'group_instrument', 'within group and instrument string'),
        fe(d, 'native', 'pi_group', 'within group'), fe(d, 'native', 'group_instrument', 'within group and instrument string'),
        fe(d, 'native_any', 'pi_group', 'within group (text or file evidence)')]
fe_t = pd.DataFrame(rows); fe_t.to_csv(OUT / 'Fixed-effect-estimates.csv', index=False)

# B. same group before / after its first DIA dataset (any software). Primary design.
#    The same design anchored on the first DIA-NN/Spectronaut dataset is kept as a sensitivity.
d['seq'] = d.groupby('pi_group').cumcount()
def switchers(expo, suffix, PH):
    first = d[d[expo] == 1].groupby('pi_group').seq.min().rename('switch_seq')
    s = d.merge(first, left_on='pi_group', right_index=True)
    pre_n = s[s.seq < s.switch_seq].groupby('pi_group').size()
    post_other = s[(s.seq > s.switch_seq) & (s[expo] == 0)].groupby('pi_group').size()
    eligible = pre_n[pre_n >= 2].index.intersection(post_other[post_other >= 1].index)
    s = s[s.pi_group.isin(eligible)].copy()
    s['phase'] = np.select([s.seq < s.switch_seq, s[expo] == 1], PH[:2], PH[2])
    gm = s.groupby(['pi_group', 'phase'])['std'].mean().unstack()[PH]
    rng = np.random.default_rng(SEED)
    B = pd.DataFrame([gm.iloc[rng.integers(0, len(gm), len(gm))].mean() for _ in range(2000)])
    sw = pd.DataFrame({'groups': gm.notna().sum(), 'deposits': s.groupby('phase').size(), 'equal_group_rate': gm.mean(),
                       'ci_low': B.quantile(.025), 'ci_high': B.quantile(.975), 'deposit_rate': s.groupby('phase')['std'].mean()})
    dd = []
    for k, v in {'later_other_minus_earlier': gm[PH[2]] - gm[PH[0]], 'exposed_minus_earlier': gm[PH[1]] - gm[PH[0]]}.items():
        bb = [v.iloc[rng.integers(0, len(v), len(v))].mean() for _ in range(2000)]
        dd.append({'contrast': k, 'groups': len(v), 'mean_paired_difference': v.mean(), 'ci_low': np.quantile(bb, .025), 'ci_high': np.quantile(bb, .975)})
    sw.to_csv(OUT / f'Switcher-phases{suffix}.csv'); pd.DataFrame(dd).to_csv(OUT / f'Switcher-contrasts{suffix}.csv', index=False)
    strong = gm[gm[PH[0]] >= .5]
    pd.DataFrame({'groups': [len(strong)], **{c: [strong[c].mean()] for c in gm.columns}}).to_csv(OUT / f'Switchers-with-strong-prior-habit{suffix}.csv', index=False)
    s[['pi_group', 'accession', 'release_year', 'seq', 'switch_seq', expo, 'standards', 'phase']].rename(columns={expo: 'exposed'}).to_csv(
        OUT / f'Switcher-sequences{suffix}.csv', index=False)
    return sw, pd.DataFrame(dd)


sw, dd = switchers('dia', '', ['earlier datasets', 'DIA datasets (from first)', 'later DDA datasets'])
sw_n, dd_n = switchers('native', '-DIA-NN-or-Spectronaut', ['earlier datasets', 'DIA-NN/Spectronaut datasets (from first)', 'later datasets with other software'])

# C. habit x tool: rate by prior history and current stack -------------------------------
d['prior_n'] = d.seq
d['prior_pos'] = d.groupby('pi_group')['std'].transform(lambda q: q.cumsum().shift(1).fillna(0))
d['history'] = np.select([d.prior_n < 2, d.prior_pos == 0, d.prior_pos == d.prior_n], ['fewer than 2 prior', 'none', 'all'], 'some')
# DIA datasets form one category whatever software they name; DDA datasets by software named.
d['stack'] = np.where(d.dia_expanded, 'DIA (any software)', d.family.replace({
    'DIA-NN': 'Several or other tools', 'Spectronaut': 'Several or other tools',
    'Several of these': 'Several or other tools', 'Other tools only': 'Several or other tools'}))
r = d[d.release_year >= 2019]; ht = []
for (h, st), g in r.groupby(['history', 'stack']):
    lo, hi = cluster_boot(g['std'].to_numpy(float), g.pi_group.to_numpy(), reps=1000) if len(g) >= 5 else (np.nan, np.nan)
    ht.append({'history': h, 'stack': st, 'deposits': len(g), 'groups': g.pi_group.nunique(), 'rate': g['std'].mean(), 'ci_low': lo, 'ci_high': hi})
ht = pd.DataFrame(ht); ht.to_csv(OUT / 'Habit-by-stack.csv', index=False)
# additive decomposition of predictive weight: history vs stack on the same deposits
m1 = smf.logit('std ~ C(history) + C(release_year)', r).fit(disp=0)
m2 = smf.logit('std ~ C(stack) + C(release_year)', r).fit(disp=0)
m3 = smf.logit('std ~ C(history) + C(stack) + C(release_year)', r).fit(disp=0)
m0 = smf.logit('std ~ C(release_year)', r).fit(disp=0)
pd.DataFrame({'model': ['year', 'year+history', 'year+stack', 'year+history+stack'],
              'mcfadden_r2': [m0.prsquared, m1.prsquared, m2.prsquared, m3.prsquared],
              'deposits': len(r)}).to_csv(OUT / 'History-vs-stack-fit.csv', index=False)
print(fe_t.round(4).to_string(index=False)); print(sw.round(3).to_string()); print(dd.round(4).to_string(index=False))
print('DIA-NN/Spectronaut-anchored:'); print(sw_n.round(3).to_string()); print(dd_n.round(4).to_string(index=False))
print(pd.read_csv(OUT / 'Switchers-with-strong-prior-habit.csv').round(3).to_string(index=False))
print(ht.pivot(index='stack', columns='history', values='rate').round(3)); print(ht.pivot(index='stack', columns='history', values='deposits'))
print(pd.read_csv(OUT / 'History-vs-stack-fit.csv').round(4))
