"""Source tables for the dense figures (written to outputs/Stack-and-reuse/Figure-data/).

Each table is a plain aggregation of Deposit-features plus the shared software flags; the
figures read only these tables and the analysis outputs. Cells keep their denominators so
masking thresholds are applied at plotting time and remain auditable.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from common import O, features, FAMILY_ORDER

OUT = O / 'Figure-data'; OUT.mkdir(exist_ok=True)
d = features()
d['vendor_single'] = d.vendor_single.fillna('Unknown')
d['acquisition'] = np.where(d.dia_expanded, 'DIA', 'DDA/other')
# Display folding of the exclusive partition into eight classes (identity colors <= 8).
FOLD = {k: k for k in ['Mascot', 'Proteome Discoverer', 'MaxQuant', 'FragPipe', 'Spectronaut', 'DIA-NN', 'No tool named']}
FOLD.update({'Several of these': 'Several or other tools', 'Other tools only': 'Several or other tools'})
d['family8'] = d.family.map(FOLD)

# per-tool (non-exclusive: a deposit counts under every tool it names) x year standard-result rate
from common import SOFTWARE
TOOLS = {'Mascot': 'mascot', 'Scaffold': 'scaffold', 'PEAKS': 'peaks', 'Proteome Discoverer': 'pd', 'MaxQuant': 'maxquant',
         'Skyline': 'skyline', 'FragPipe': 'msfragger', 'Spectronaut': 'spectronaut', 'DIA-NN': 'diann'}
rows = []
for name, k in TOOLS.items():
    x = d[d['sw_' + k]]
    rows += [{'tool': name, 'release_year': y, 'deposits': len(g), 'standards': g.standards.mean()} for y, g in x.groupby('release_year')]
    r = x[x.release_year.between(2019, 2025)]
    rows.append({'tool': name, 'release_year': '2019-2025', 'deposits': len(r), 'standards': r.standards.mean()})
none = d[~d[['sw_' + k for k in SOFTWARE]].any(axis=1)]
rows += [{'tool': 'No tool named', 'release_year': y, 'deposits': len(g), 'standards': g.standards.mean()} for y, g in none.groupby('release_year')]
r = none[none.release_year.between(2019, 2025)]; rows.append({'tool': 'No tool named', 'release_year': '2019-2025', 'deposits': len(r), 'standards': r.standards.mean()})
pd.DataFrame(rows).to_csv(OUT / 'Tool-year.csv', index=False)
# family x year: share and standards rate with denominators
fy = d.groupby(['release_year', 'family']).agg(deposits=('accession', 'size'), standards=('standards', 'mean')).reset_index()
fy['share'] = fy.deposits / fy.groupby('release_year').deposits.transform('sum'); fy.to_csv(OUT / 'Family-year.csv', index=False)
f8 = d.groupby(['release_year', 'family8']).size().rename('deposits').reset_index()
f8['share'] = f8.deposits / f8.groupby('release_year').deposits.transform('sum'); f8.to_csv(OUT / 'Family8-year.csv', index=False)
# acquisition x year counts
d.groupby(['release_year', 'acquisition']).agg(deposits=('accession', 'size'), standards=('standards', 'sum')).reset_index().to_csv(OUT / 'Acquisition-year.csv', index=False)
# vendor x year x family8 (single-vendor Thermo/Bruker/SCIEX/Waters)
vy = d[d.vendor_single.isin(['Thermo Fisher', 'Bruker', 'SCIEX', 'Waters'])].groupby(['vendor_single', 'release_year', 'family8']).size().rename('deposits').reset_index()
vy['share'] = vy.deposits / vy.groupby(['vendor_single', 'release_year']).deposits.transform('sum'); vy.to_csv(OUT / 'Vendor-year-family8.csv', index=False)
vr = d[d.vendor_single.isin(['Thermo Fisher', 'Bruker', 'SCIEX', 'Waters'])].groupby(['vendor_single', 'release_year']).agg(
    deposits=('accession', 'size'), standards=('standards', 'mean'), mzml=('fmt_mzml', 'mean'), dia=('dia_expanded', 'mean'),
    dia_native=('dia_native_tool', 'mean')).reset_index(); vr.to_csv(OUT / 'Vendor-year-rates.csv', index=False)
# country x family8 composition and rates, 2019-2025, equal PI-group weight within country
c = d[d.release_year.between(2019, 2025) & d.country.notna()].copy()
c['cl'] = c.pi_group.fillna('solo:' + c.accession)
c['w'] = 1 / c.groupby(['country', 'cl']).accession.transform('size')
ct = c.groupby('country').agg(deposits=('accession', 'size'), groups=('cl', 'nunique'), standards=('standards', 'mean'),
                              dia_native=('dia_native_tool', 'mean'))
cf = (c.groupby(['country', 'family8']).w.sum() / c.groupby('country').w.sum()).rename('equal_group_share').reset_index()
cf = cf.merge(ct, left_on='country', right_index=True); cf.to_csv(OUT / 'Country-family8.csv', index=False)

# same-group raster: groups with >=2 datasets before their first DIA dataset and >=1 later DDA dataset
g = d[d.pi_group.notna()].sort_values(['pi_group', 'release_date', 'accession']).copy()
g['seq'] = g.groupby('pi_group').cumcount()
first = g[g.dia_expanded].groupby('pi_group').seq.min()
g = g[g.pi_group.isin(first.index)].copy(); g['rel'] = g.seq - g.pi_group.map(first)
pre = g[g.rel < 0].groupby('pi_group').size(); post = g[(g.rel > 0) & ~g.dia_expanded].groupby('pi_group').size()
elig = pre[pre >= 2].index.intersection(post.index)
g = g[g.pi_group.isin(elig)]
assert g.pi_group.nunique() == 244, g.pi_group.nunique()
g = g.assign(exposed=g.dia_expanded)
g[['pi_group', 'accession', 'release_year', 'rel', 'exposed', 'dia_native_tool', 'dia_expanded', 'standards', 'vendor_single', 'family8']].to_csv(
    OUT / 'Switcher-sequences.csv', index=False)
print('figure data written:', sorted(p.name for p in OUT.iterdir()))

# per-country values for the world maps (denominators kept; display thresholds applied at plotting)
c_all = d[d.country.notna()]
m = c_all.groupby('country').agg(datasets_2013_2025=('accession', 'size'), groups_2013_2025=('pi_group', 'nunique'))
r1 = c_all[c_all.release_year.between(2019, 2025)].groupby('country').agg(datasets_2019_2025=('accession', 'size'), with_files_2019_2025=('standards', 'mean'))
r2 = c_all[c_all.release_year.between(2023, 2025)].groupby('country').agg(datasets_2023_2025=('accession', 'size'), dia_nn_or_spectronaut_2023_2025=('dia_native_tool', 'mean'), dia_2023_2025=('dia_expanded', 'mean'))
m.join(r1).join(r2).reset_index().to_csv(OUT / 'Country-map.csv', index=False)
# overall rates in the two windows used for the step chart (Fig. 2c)
pd.DataFrame({'window': ['2019-2021', '2023-2025'],
              'deposits': [int(d.release_year.between(2019, 2021).sum()), int(d.release_year.between(2023, 2025).sum())],
              'with_files': [d[d.release_year.between(2019, 2021)].standards.mean(), d[d.release_year.between(2023, 2025)].standards.mean()]}
             ).to_csv(OUT / 'Window-rates.csv', index=False)
