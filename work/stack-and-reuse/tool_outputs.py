"""For datasets that name a tool, what did they deposit — and does it depend on what the tool
can export?

Per (dataset, tool) outcome, files outside /generated/ only, first match wins:
  standard  — mzIdentML/mzTab listed
  native    — at least one of the tool's documented native result files is listed
  archive   — only a tool-named or unlabelled results archive (contents not visible)
  tables    — only other tables
  none      — nothing above
Capability comparison: share of datasets with a standard file by the tool's export class
(native / built-in optional / none), for datasets naming exactly one registry tool, released
from the year after the export was first documented where a date is known. Descriptive only.
"""
import numpy as np, pandas as pd
from common import R, O, SEED, cluster_boot
from software_registry import compile_registry, TOOLS
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Tool-outputs'; OUT.mkdir(exist_ok=True)
REG = compile_registry()
c = connect()
files = c.execute(r"""
select project_accession accession, lower(regexp_replace(file_name, '\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
       regexp_matches(lower(file_name), '\.(zip|tgz|tar|rar|7z)$|\.tar\.gz$') is_archive
from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
  and project_accession in (select accession from study)""").df()
names = files.groupby('accession').b.apply(list); arch = files[files.is_archive].groupby('accession').b.apply(list)
L = pd.read_parquet(O / 'Layers/Dataset-layers.parquet').set_index('accession')
F = pd.read_parquet(O / 'Deposit-features.parquet', columns=['accession', 'pi_group']).set_index('accession')
txt = connect().execute('select accession, software, processing from study').df().set_index('accession')
text = (txt.software.fillna('') + ' ' + txt.processing.fillna('')).str.lower()

named = {tool: set(text[text.str.contains(t['text'], regex=True)].index) for tool, t in REG.items()}
n_named = pd.Series(0, index=text.index)
for s in named.values():
    n_named[list(s)] += 1
rows = []
for tool, t in REG.items():
    for acc in named[tool]:
        nm = names.get(acc, []); ar = arch.get(acc, [])
        native = any(t['native_re'].search(x) for x in nm)
        std = bool(L.at[acc, 'has_standard_results'])
        tool_arch = any(t['archive_re'].search(x) for x in ar)
        outcome = ('standard' if std else 'native' if native else 'archive' if (tool_arch or L.at[acc, 'has_archived_results'])
                   else 'tables' if L.at[acc, 'has_unidentified_tables'] else 'none')
        rows.append({'accession': acc, 'tool': tool, 'export': t['export'] or 'unknown', 'since': t['since'],
                     'native_visible': native, 'standard': std, 'tool_archive': tool_arch, 'outcome': outcome,
                     'tools_named': int(n_named[acc])})
d = pd.DataFrame(rows).merge(L[['release_year', 'dia_expanded']], left_on='accession', right_index=True)
d = d.merge(F, left_on='accession', right_index=True)
d.to_csv(OUT / 'Dataset-tool-outcomes.csv.gz', index=False)

r = d[d.release_year >= 2019]
order = ['standard', 'native', 'archive', 'tables', 'none']
tab = pd.crosstab(r.tool, r.outcome, normalize='index').reindex(columns=order).fillna(0).mul(100).round(1)
tab['datasets'] = r.groupby('tool').size(); tab['native_visible_any'] = r.groupby('tool').native_visible.mean().mul(100).round(1)
tab['export'] = [TOOLS[k]['export'] or 'unknown' for k in tab.index]; tab['formats'] = [TOOLS[k]['formats'] for k in tab.index]
tab['evidence'] = [TOOLS[k]['grade'] for k in tab.index]
tab = tab.sort_values('standard', ascending=False); tab.to_csv(OUT / 'Tool-outcomes-2019-2025.csv')

# capability comparison: datasets naming exactly one registry tool, after the export existed
s = d[(d.tools_named == 1) & (d.release_year >= 2019)].copy()
s = s[s.since.isna() | (s.release_year > s.since)]
s['cluster'] = s.pi_group.fillna('solo:' + s.accession)
cap = []
for cls, g in s.groupby('export'):
    lo, hi = cluster_boot(g.standard.astype(float).to_numpy(), g.cluster.to_numpy(), reps=1000)
    cap.append({'export_class': cls, 'datasets': len(g), 'groups': g.cluster.nunique(), 'tools': ', '.join(sorted(g.tool.unique())),
                'standard_pct': 100 * g.standard.mean(), 'ci_low': 100 * lo, 'ci_high': 100 * hi,
                'native_visible_pct': 100 * g.native_visible.mean(), 'archive_only_pct': 100 * (g.outcome == 'archive').mean()})
cap = pd.DataFrame(cap); cap.to_csv(OUT / 'Standard-rate-by-export-capability.csv', index=False)
per = s.groupby(['export', 'tool']).agg(datasets=('accession', 'size'), standard_pct=('standard', lambda x: 100 * x.mean()),
                                        native_visible_pct=('native_visible', lambda x: 100 * x.mean())).round(1)
per.to_csv(OUT / 'Standard-rate-by-tool-single-tool-datasets.csv')
pd.set_option('display.width', 250)
print(tab.drop(columns=['formats']).to_string()); print(cap.round(1).to_string(index=False)); print(per.to_string())
