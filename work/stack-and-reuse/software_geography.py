"""Do analysis tools have home countries?

Hypothesis fixed before computing (October 2026): a tool is named more often, relative to
its archive-wide share, in deposits whose PI country is the developer's country. The
developer-country map below was written before any country x tool rates were examined.
Tools with ambiguous or vendor-global development (Proteome Discoverer, ProteinPilot, PLGS)
are excluded rather than assigned.

Equal PI-group weighting within country is primary (a prolific lab is one unit); deposit
weighting is a sensitivity. Lift = country share / all-country share. Intervals resample
PI groups within country. Descriptive association; no claim about why.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from common import O, SEED, SOFTWARE, features

OUT = O / 'Software-geography'; OUT.mkdir(exist_ok=True)
HOME = {  # tool flag -> developer country/countries (institution in comment)
    'maxquant': ['Germany'],                         # MPI of Biochemistry, Martinsried
    'peptideshaker': ['Belgium', 'Norway'],          # Ghent University / VIB; University of Bergen
    'spectronaut': ['Switzerland'],                  # Biognosys
    'diann': ['United Kingdom', 'Germany'],          # Cambridge / Francis Crick; Charite Berlin
    'msfragger': ['United States'],                  # University of Michigan
    'peaks': ['Canada'],                             # Bioinformatics Solutions, Waterloo
    'mascot': ['United Kingdom'],                    # Matrix Science, London
    'scaffold': ['United States'],                   # Proteome Software, Portland
    'skyline': ['United States'],                    # University of Washington
    'pfind': ['China'],                              # Institute of Computing Technology, CAS
    'openms': ['Germany', 'Switzerland'],            # Tuebingen / FU Berlin; ETH Zurich (OpenSWATH)
    'comet': ['United States'], 'msgf': ['United States'], 'byonic': ['United States'],
    'metamorpheus': ['United States'], 'xtandem': ['Canada'], 'progenesis': ['United Kingdom'],
}
d = features()
d = d[d.release_year.between(2019, 2025) & d.country.notna() & d.pi_group.notna()].copy()
ct = d.groupby('country').pi_group.nunique(); keep = ct[ct >= 60].index
d = d[d.country.isin(keep)].copy()
d['w'] = 1 / d.groupby(['country', 'pi_group']).accession.transform('size')  # equal group within country
# Composition metric: each deposit splits one unit of credit across all named tools (as for
# vendors), so a country that simply names more tools per deposit gains no lift.
named = d[['sw_' + k for k in SOFTWARE]].sum(axis=1)
METRIC = __import__('os').environ.get('SWGEO_METRIC', 'presence')
if METRIC == 'composition':
    d = d[named > 0].copy(); named = named[named > 0]
    d['w'] = 1 / d.groupby(['country', 'pi_group']).accession.transform('size')
    OUT = OUT / 'composition'; OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(SEED)
rows = []
for tool, homes in HOME.items():
    f = d['sw_' + tool].astype(float) / (named if METRIC == 'composition' else 1)
    # global share: average of group-level rates over all groups (equal group weighting)
    gl = (f * d.w).groupby(d.pi_group).sum().sum() / d.w.groupby(d.pi_group).sum().sum()
    gl_dep = f.mean()
    for c, x in d.groupby('country'):
        gr = (f[x.index] * x.w).groupby(x.pi_group).sum()  # per-group rate (weights sum to 1 per group)
        share = gr.mean(); boots = [gr.iloc[rng.integers(0, len(gr), len(gr))].mean() for _ in range(500)]
        rows.append({'tool': tool, 'country': c, 'home': c in homes, 'groups': len(gr), 'deposits': len(x),
                     'groups_naming_tool': int((gr > 0).sum()), 'equal_group_share': share, 'global_equal_group_share': gl,
                     'lift': share / gl if gl else np.nan, 'lift_ci_low': np.quantile(boots, .025) / gl, 'lift_ci_high': np.quantile(boots, .975) / gl,
                     'deposit_share': f[x.index].mean(), 'deposit_lift': f[x.index].mean() / gl_dep})
t = pd.DataFrame(rows)
t['rank_among_countries'] = t.groupby('tool').lift.rank(ascending=False, method='min')
t['countries'] = t.groupby('tool').country.transform('size')
t.to_csv(OUT / 'Tool-country-lift.csv', index=False)
h = t[t.home].copy()
h['home_is_top3'] = h.rank_among_countries <= 3
h.to_csv(OUT / 'Home-country-lift.csv', index=False)
# permutation-style reference: probability a random country ranks in the top 3 = 3/countries
summary = {'tools': h.tool.nunique(), 'tool_home_pairs': len(h), 'countries': int(t.country.nunique()),
           'home_pairs_lift_above_1': int((h.lift > 1).sum()), 'home_pairs_ci_above_1': int((h.lift_ci_low > 1).sum()),
           'home_pairs_top3': int(h.home_is_top3.sum()), 'expected_top3_by_chance': float((3 / h.countries).sum()),
           'median_home_lift': float(h.lift.median()), 'median_home_deposit_lift': float(h.deposit_lift.median())}
# exact-ish reference: shuffle home labels across countries 10,000 times (same number of homes per tool)
L = t.pivot(index='tool', columns='country', values='lift'); obs = h.lift.apply(np.log).median(); null = []
for _ in range(10000):
    v = []
    for tool, homes in HOME.items():
        row = L.loc[tool].dropna(); v += list(np.log(rng.choice(row.to_numpy(), size=sum(c in row.index for c in homes), replace=False)))
    null.append(np.median(v))
summary['permutation_p_median_log_lift'] = float(((np.array(null) >= obs).sum() + 1) / (len(null) + 1))
summary['metric'] = METRIC
pd.Series(summary).to_csv(OUT / 'Home-country-summary.csv')
print(h[['tool', 'country', 'groups', 'groups_naming_tool', 'equal_group_share', 'global_equal_group_share', 'lift', 'lift_ci_low', 'lift_ci_high', 'rank_among_countries', 'countries', 'deposit_lift']].round(3).to_string(index=False))
print(pd.Series(summary))

# Sensitivity: the United States is the developer country for 7 of 20 pairs. Repeat the
# reference comparison for the 13 non-US pairs only.
nh = h[h.country != 'United States']; obs = np.log(nh.lift.clip(lower=1e-9)).median(); null = []
for _ in range(10000):
    null.append(np.median([np.log(max(rng.choice(L.loc[tl].dropna().to_numpy()), 1e-9)) for tl in nh.tool]))
pd.Series({'metric': METRIC, 'non_us_pairs': len(nh), 'lift_above_1': int((nh.lift > 1).sum()), 'median_lift': float(nh.lift.median()),
           'permutation_p': float(((np.array(null) >= obs).sum() + 1) / (len(null) + 1))}).to_csv(OUT / 'Home-country-summary-non-US.csv')
