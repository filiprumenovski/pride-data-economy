"""42-breadth: which axes of the archive broadened and which narrowed: effective number of contributor groups, journals, instrument
models and vendors, analysis tools, research subfields, countries and species behind PRIDE datasets (descriptive).

Effective number = inverse Simpson index 1 / sum(p_k^2) of the categories' shares among a release year's
datasets (dataset grain; a dataset naming several models, tools or species contributes fractional credit
1/n to each; instrument labels harmonised to models as in Arsenal H3, non-model labels excluded). Unlike a count of categories it barely depends on sample size. Datasets without the attribute
are left out of that axis (denominators in the figure data). 95% intervals: bootstrap over contributor
groups (pi_group; solo datasets are their own group), 300 replicates, seed 20261001.
Panel b: pooled 2013–15 vs 2023–25 ratios of Hill numbers of order q = 2 (inverse Simpson), q = 1 (exp Shannon) and
q = 0 (richness rarefied to the smaller period's dataset count, mean of 50 subsamples), with 95% intervals from
bootstrap resampling of contributor groups within each period (q = 1, 2; 300 replicates). Two further axes are added
there: contributor groups (each dataset credited to its group) and instrument vendors (harmonised model vendor).
Outputs: Gallery/42-breadth.{png,pdf}, Gallery/data/42-breadth*.csv.
"""
import sys, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from common import O, R, SEED, features
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3
sys.path.insert(0, str(R / 'scripts'))
from query import connect

YRS = list(range(2016, 2026))   # 2016 baseline: every axis has >900 datasets (2014: 336–483, wide intervals)
REPS = 300
rng = np.random.default_rng(SEED)
c = connect()
D = features()
D = D[D.release_year.between(2013, 2025)].copy()
D['grp'] = D.pi_group.fillna('solo:' + D.accession)


def long_frame():
    """(accession, axis, category, weight) rows; weights sum to 1 per dataset and axis."""
    rows = []
    # countries
    x = D[D.country.notna()][['accession', 'country']].rename(columns={'country': 'cat'}).assign(w=1.0, axis='Countries')
    rows.append(x)
    # instrument models (fractional, as in the fleet analysis)
    # harmonised model names, as in the fleet analysis (Arsenal H3): labels that are models only
    lab = pd.read_csv(O / 'Arsenal/H1-H3-scale-instruments/H3-label-harmonisation.csv')[['acc', 'nm', 'label_type', 'model']]
    m = c.execute("select distinct project_accession accession, instrument_accession acc, instrument_name nm from pride_instruments").df()
    m = m[m.accession.isin(D.accession)].merge(lab, on=['acc', 'nm'], how='left')
    m = m[m.label_type == 'model'][['accession', 'model']].rename(columns={'model': 'cat'}).drop_duplicates()
    m['w'] = 1 / m.groupby('accession').cat.transform('size'); rows.append(m.assign(axis='Instrument models'))
    vend = pd.read_csv(O / 'Arsenal/H1-H3-scale-instruments/H3-label-harmonisation.csv')[['model', 'model_vendor']].dropna().drop_duplicates('model')
    v = m.merge(vend, left_on='cat', right_on='model')[['accession', 'model_vendor', 'w']].rename(columns={'model_vendor': 'cat'})
    v = v.groupby(['accession', 'cat'], as_index=False).w.sum(); rows.append(v.assign(axis='Instrument vendors'))
    rows.append(D[['accession', 'grp']].rename(columns={'grp': 'cat'}).assign(w=1.0, axis='Contributor groups'))
    # analysis tools: the 22 flagged tools named in the deposit
    sw = [k for k in D.columns if k.startswith('sw_')]
    t = D[['accession'] + sw].melt('accession', var_name='cat', value_name='on')
    t = t[t.on.astype(bool)].drop(columns='on'); t['w'] = 1 / t.groupby('accession').cat.transform('size')
    rows.append(t.assign(axis='Analysis tools'))
    # species, strains collapsed to the first two words
    s = c.execute("select project_accession accession, species_name from pride_species").df()
    s = s[s.accession.isin(D.accession)]
    s['cat'] = s.species_name.str.replace(r'\(.*?\)', '', regex=True).str.strip().str.split().str[:2].str.join(' ').str.lower()
    s = s.drop(columns='species_name').drop_duplicates(); s['w'] = 1 / s.groupby('accession').cat.transform('size')
    rows.append(s.assign(axis='Species'))
    # research subfield and journal of the linked paper
    rows.append(D[D.primary_subfield.notna()][['accession', 'primary_subfield']].rename(columns={'primary_subfield': 'cat'})
                .assign(w=1.0, axis='Research subfields'))
    j = c.execute("""select w.work_id, lower(trim(e.journal)) cat from lit.openalex_works w
                     join lit.europepmc_works e on w.pmid = e.pmid where e.journal is not null""").df().drop_duplicates('work_id')
    jj = D[['accession', 'work_id']].merge(j, on='work_id')[['accession', 'cat']].assign(w=1.0, axis='Journals')
    rows.append(jj)
    L = pd.concat(rows, ignore_index=True)
    return L.merge(D[['accession', 'release_year', 'grp']], on='accession')


def hill(w_by_cat, q):
    p = w_by_cat / w_by_cat.sum(-1, keepdims=True)
    if q == 2:
        return 1 / (p ** 2).sum(-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.exp(-np.nansum(np.where(p > 0, p * np.log(p), 0), -1))


def pooled(L, early=(2013, 2015), late=(2023, 2025), reps=REPS, sub=50):
    # DISTINCT SQL rows have no guaranteed order. Fix category/group and sampling
    # order before seeded resampling; frozen figure data remain a separate input.
    L = L.sort_values(['axis', 'accession', 'cat', 'grp']).reset_index(drop=True)
    out = []
    for ax, A_ in L.groupby('axis'):
        P = {}
        for lab, (a, b) in {'early': early, 'late': late}.items():
            Y = A_[A_.release_year.between(a, b)]
            g = pd.factorize(Y.grp)[0]
            P[lab] = (Y, g)
        cats = pd.Index(pd.concat([P['early'][0].cat, P['late'][0].cat]).unique())
        mats = {}
        for lab, (Y, g) in P.items():
            k = cats.get_indexer(Y.cat)
            M = np.zeros((g.max() + 1, len(cats))); np.add.at(M, (g, k), Y.w.to_numpy()); mats[lab] = M
        row = {'axis': ax, 'datasets_early': int(P['early'][0].accession.nunique()), 'datasets_late': int(P['late'][0].accession.nunique())}
        for q in (2, 1):
            e = hill(mats['early'].sum(0), q); l = hill(mats['late'].sum(0), q)
            bs = []
            for _ in range(reps):
                we = rng.multinomial(len(mats['early']), np.full(len(mats['early']), 1 / len(mats['early'])))
                wl = rng.multinomial(len(mats['late']), np.full(len(mats['late']), 1 / len(mats['late'])))
                bs.append(hill(wl @ mats['late'], q) / hill(we @ mats['early'], q))
            row.update({f'q{q}_early': e, f'q{q}_late': l, f'q{q}_ratio': l / e,
                        f'q{q}_lo': np.percentile(bs, 2.5), f'q{q}_hi': np.percentile(bs, 97.5)})
        # q = 0: richness rarefied to the smaller period's number of datasets
        n = min(row['datasets_early'], row['datasets_late'])
        rich = {}
        for lab, (Y, g) in P.items():
            acc = Y.accession.unique(); vals = []
            for _ in range(sub):
                pick = set(rng.choice(acc, n, replace=False)); vals.append(Y[Y.accession.isin(pick)].cat.nunique())
            rich[lab] = np.mean(vals)
        row.update({'q0_early_rarefied': rich['early'], 'q0_late_rarefied': rich['late'], 'q0_ratio': rich['late'] / rich['early'], 'rarefied_to': n})
        out.append(row)
    return pd.DataFrame(out)


def inv_simpson(w_by_cat):
    p = w_by_cat / w_by_cat.sum(-1, keepdims=True)
    return 1 / (p ** 2).sum(-1)


def series(L):
    out = []
    for ax, A_ in L.groupby('axis'):
        for y in YRS:
            Y = A_[A_.release_year == y]
            if Y.accession.nunique() < 50:
                continue
            g = pd.factorize(Y.grp)[0]; k = pd.factorize(Y.cat)[0]
            M = np.zeros((g.max() + 1, k.max() + 1)); np.add.at(M, (g, k), Y.w.to_numpy())
            est = inv_simpson(M.sum(0))
            W = rng.multinomial(len(M), np.full(len(M), 1 / len(M)), size=REPS)
            bs = inv_simpson(W @ M)
            out.append({'axis': ax, 'year': y, 'datasets': int(Y.accession.nunique()), 'categories': int(k.max() + 1),
                        'eff': est, 'lo': np.percentile(bs, 2.5), 'hi': np.percentile(bs, 97.5)})
    return pd.DataFrame(out)


def concentration():
    out = []
    for lab, (a, b) in {'2013–17': (2013, 2017), '2021–25': (2021, 2025)}.items():
        P = D[D.release_year.between(a, b) & D.pi_group.notna()]
        n = P.groupby('pi_group').size().sort_values(ascending=False).to_numpy()
        for q in [.01, .10]:
            k = max(1, int(round(len(n) * q)))
            out.append({'period': lab, 'top': f'{int(q * 100)}%', 'groups': len(n), 'share': n[:k].sum() / n.sum()})
    return pd.DataFrame(out)


def fig_breadth():
    if FROZEN_DATA:
        S, C, P = [pd.read_csv(G / 'data' / f'42-breadth{suffix}.csv')
                   for suffix in ['', '-concentration', '-pooled']]
    else:
        L = long_frame(); S = series(L); C = concentration(); P = pooled(L)
    if not FROZEN_DATA:
        (G / 'data').mkdir(exist_ok=True)
        S.to_csv(G / 'data/42-breadth.csv', index=False); C.to_csv(G / 'data/42-breadth-concentration.csv', index=False)
        P.to_csv(G / 'data/42-breadth-pooled.csv', index=False)
    order = ['Contributor groups', 'Journals', 'Instrument models', 'Analysis tools', 'Research subfields', 'Countries', 'Species']
    cols = {'Contributor groups': '#1b1b1b', 'Journals': ACC, 'Instrument models': ACC2, 'Analysis tools': ACC3, 'Countries': '#5a9e6f',
            'Research subfields': '#8d6cab', 'Species': '#9aa5ad', 'Instrument vendors': '#b7a99a'}
    base = {a: S[(S.axis == a) & (S.year == 2016)].eff.iloc[0] for a in order}
    end = {a: S[(S.axis == a) & (S.year == 2025)].eff.iloc[0] for a in order}
    pr = P.set_index('axis')
    up = [a for a in pr.index if a != 'Contributor groups' and pr.q2_lo[a] > 1 and pr.q1_lo[a] > 1]
    down = [a for a in pr.index if a != 'Contributor groups' and pr.q2_hi[a] < 1 and pr.q1_hi[a] < 1]
    r = lambda a: pr.q2_ratio[a]
    fig = plt.figure(figsize=(W, 74 * MM))
    headline(fig, f'Contributor groups (×{r("Contributor groups"):.1f}) and journals (×{r("Journals"):.1f}) broadened, instrument use evened '
                  f'(×{r("Instrument models"):.1f}) and species narrowed (×{r("Species"):.2f})',
             'Effective number (Hill number) of categories among PRIDE datasets. a, order 2 (inverse Simpson) by release year relative to 2016. '
             'b, 2023–25 vs 2013–15 at orders 2, 1\nand 0 (richness rarefied to equal dataset counts); 95% intervals from bootstrap over contributor '
             'groups for axes other than contributor-group diversity, which is descriptive. Fractional credit when a dataset names several categories.')
    ax = fig.add_axes([.07, .13, .4, .6])
    pos = repel({a: np.log2(end[a] / base[a]) for a in order}, .2)
    for a in order:
        s_ = S[S.axis == a].set_index('year')
        emphasize = a == 'Contributor groups' or a in up + down
        if a != 'Contributor groups':
            ax.fill_between(s_.index, np.log2(s_.lo / base[a]), np.log2(s_.hi / base[a]), color=cols[a], alpha=.13, lw=0)
        ax.plot(s_.index, np.log2(s_.eff / base[a]), color=cols[a], lw=1.7 if emphasize else 1.0,
                alpha=1 if emphasize else .75)
        ax.text(2025.3, pos[a], a, color=cols[a], fontsize=6.4, va='center', fontweight='bold' if emphasize else 'normal')
    ax.axhline(0, color=MUTED, lw=.6, ls='--')
    ax.set_xlim(2016, 2025); ax.set_xticks([2016, 2018, 2020, 2022, 2025])
    ax.set_yticks(np.log2([.5, 1, 2, 4])); ax.set_ylim(np.log2(.55), np.log2(4.8))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'×{2 ** v:g}'))
    ax.set_ylabel('Effective number, relative to 2016'); ax.set_xlabel('Release year')
    panel(ax, 'a', x=-.1, y=1.02)
    bx = fig.add_axes([.775, .13, .21, .6])
    rows = order + ['Instrument vendors']
    rows = sorted(rows, key=lambda a: pr.q2_ratio[a])
    for i, a in enumerate(rows):
        c_ = cols[a]; sig = a == 'Contributor groups' or a in up + down
        if a != 'Contributor groups':
            bx.plot([pr.q2_lo[a], pr.q2_hi[a]], [i, i], color=c_, lw=1.4, alpha=1 if sig else .5)
        bx.plot(pr.q2_ratio[a], i, 'o', color=c_, ms=4.5, alpha=1 if sig else .6)
        bx.plot(pr.q1_ratio[a], i - .28, 'o', mfc='white', mec=c_, ms=3.2, mew=.8)
        bx.plot(pr.q0_ratio[a], i + .28, 'x', color=c_, ms=3.2, mew=.8)
    bx.axvline(1, color=MUTED, lw=.6, ls='--')
    bx.set_yticks(range(len(rows))); bx.set_yticklabels(rows, fontsize=6.4); bx.tick_params(axis='y', length=0)
    bx.set_xscale('log'); bx.set_xlim(.45, 9); bx.set_xticks([.5, 1, 2, 4, 8])
    bx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'×{v:g}')); bx.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ''))
    bx.set_xlabel('2023–25 vs 2013–15')
    bx.set_title('● order 2 (95% CI)  ○ order 1  × richness', loc='left', fontsize=6.2, color=MUTED)
    panel(bx, 'b', x=-.66, y=1.02)
    save(fig, '42-breadth')
    return S, C, P


if __name__ == '__main__':
    S, C, P = fig_breadth()
    print(S[S.year.isin([2016, 2025])].round(2).to_string()); print(P.round(2).to_string())
