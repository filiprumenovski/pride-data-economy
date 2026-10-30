"""Attack on H2.9-H2.11: release timing relative to the linked paper (new experiment, October 2026).

Descriptive only. Reads the frozen databases read-only; writes only to
outputs/Stack-and-reuse/Attacks/release-timing/.

Claims under test (outputs/Stack-and-reuse/Arsenal/H2-time/Ledger-entries.md):
  H2.9/H2.10  46.8% of datasets with a dated paper become public >30 d after the paper; 22-39% (2013-20)
              -> 51-62% (2021-25); since 2018, 62 days with >=50 releases cover 22.6% of datasets and 83%
              of those are released >30 d after the paper.
  H2.11       submission precedes the paper by a median 141 d, longer in later release years.

Sections
  A  What the release date is: project vs file dates, ProteomeXchange announce date, submission date vs
     PXD accession order (accessions are issued sequentially, so a submission date far later than that of
     accession neighbours indicates an overwritten date).
  B  Bulk release days: calendar, weekday, registry agreement, composition (submission years, countries,
     contributor groups, submission type, paper link) and threshold sensitivity.
  C  Paper-date validity: OpenAlex vs Europe PMC dates, date precision, alternative paper-date definitions.
  D  Reframe: share public by the time the paper appeared, by paper year (primary axis), release year,
     submission year, country, journal and DOI registrant; all datasets as denominator as well as
     datasets with a dated paper.
  E  Submission -> paper interval by paper year / release year, with accession-order sensitivity.

Intervals: shares use common.cluster_boot (2,000 resamples of whole pi_group clusters, seed 20261001);
datasets without a group hash are their own cluster.
"""
from pathlib import Path
import sys, json
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import R, O, SEED, features, cluster_boot
sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Attacks/release-timing'; OUT.mkdir(parents=True, exist_ok=True)
S = {'claims': 'H2.9-H2.11 of outputs/Stack-and-reuse/Arsenal/H2-time/Ledger-entries.md',
     'bootstrap': '2000 resamples of whole pi_group clusters (common.cluster_boot), seed %d' % SEED}
BULK = 50            # H2 definition: a calendar day with >= 50 cohort releases
WIN = 30             # +/- window around the paper date
WD = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']


def save(df, name):
    df.to_csv(OUT / name, index=False)


def sh(x, flag):
    """share, lo, hi with whole-group resampling."""
    flag = pd.Series(flag, index=x.index).astype(float)
    if len(x) < 2:
        return float(flag.mean()) if len(x) else np.nan, np.nan, np.nan
    lo, hi = cluster_boot(flag.to_numpy(), x.grp.to_numpy())
    return float(flag.mean()), float(lo), float(hi)


c = connect(); c.execute("set threads=1; set memory_limit='2GB'")

# ------------------------------------------------------------------ data -------------------------
d = features()
assert len(d) == 36249
d['grp'] = d.pi_group.fillna('acc:' + d.accession)
d['sub_year'] = d.submission_date.dt.year
d['release_date'] = pd.to_datetime(d.release_date); d['submission_date'] = pd.to_datetime(d.submission_date)

proj = c.execute("""
select p.project_accession accession, p.publication_date pride_release, p.submission_date pride_sub,
       try_cast(g.registry_announce_date as date) px_announce, b.publicationDate bronze_pub
from pride_projects p join g.deposits g on g.accession = p.project_accession
join pride_bronze_projects b on b.accession = p.project_accession
where p.project_accession in (select accession from study)""").df()
files = c.execute("""
select project_accession accession, count(*) n_files,
  min(publication_date)::date f_pub_min, max(publication_date)::date f_pub_max,
  min(submission_date)::date f_sub_min, max(submission_date)::date f_sub_max,
  min(updated_date)::date f_upd_min, max(updated_date)::date f_upd_max
from pride_files where project_accession in (select accession from study) group by 1""").df()
# accession-order implied submission date (all PRIDE PXD projects, rolling median of 201 neighbours)
seq = c.execute("""select project_accession accession, cast(substr(project_accession, 4) as int) num,
  date_diff('day', date '1970-01-01', submission_date) sd from pride_projects
  where project_accession like 'PXD%' order by num""").df()
seq['sd_seq'] = seq.sd.rolling(201, center=True, min_periods=50).median()
seq['sub_dev'] = seq.sd - seq.sd_seq
# paper dates: primary work (cohort work_id) in OpenAlex and Europe PMC
works = c.execute("""
select s.accession, o.work_id, o.publication_date oa_date, o.work_type, o.doi, o.pmid,
       e.first_publication_date epmc_date, e.journal
from study s join lit.openalex_works o using (work_id)
left join lit.europepmc_works e on e.pmid = o.pmid""").df()
assert works.accession.is_unique
nlinks = c.execute("""select accession, count(distinct work_id) n_linked_works,
  min(w.publication_date) earliest_any from lit.project_publications pp join lit.openalex_works w using (work_id)
  where accession in (select accession from study) group by 1""").df()
refs = c.execute(r"""select project_accession accession,
  bool_or(regexp_matches(lower(reference_line), 'in press|submitted|accepted|under review|in revision|unpublished'))
    ref_pending_text from pride_publications where project_accession in (select accession from study) group by 1""").df()

d = (d.merge(proj, on='accession', how='left', validate='one_to_one')
      .merge(files, on='accession', how='left', validate='one_to_one')
      .merge(seq[['accession', 'sd_seq', 'sub_dev']], on='accession', how='left', validate='one_to_one')
      .merge(works.drop(columns='work_id'), on='accession', how='left', validate='one_to_one')
      .merge(nlinks, on='accession', how='left', validate='one_to_one')
      .merge(refs, on='accession', how='left', validate='one_to_one'))
for col in ['pride_release', 'pride_sub', 'px_announce', 'f_pub_min', 'f_pub_max', 'f_sub_min', 'f_upd_min',
            'f_upd_max', 'oa_date', 'epmc_date', 'earliest_any']:
    d[col] = pd.to_datetime(d[col])
d['sub_seq_date'] = pd.Timestamp('1970-01-01') + pd.to_timedelta(d.sd_seq.round(), unit='D')

# ------------------------------------------------------------------ A. what the release date is ---
A = []
def a_row(check, num, den, note=''):
    A.append({'check': check, 'datasets': int(num), 'denominator': int(den), 'share': num / den if den else np.nan,
              'note': note})
n = len(d)
a_row('features release_date == pride_projects.publication_date', (d.release_date == d.pride_release).sum(), n)
a_row('bronze publicationDate string == pride_projects.publication_date',
      (pd.to_datetime(d.bronze_pub, errors='coerce') == d.pride_release).sum(), n)
hasf = d.f_pub_min.notna()
a_row('all file publication_date == project publication_date (min and max)',
      ((d.f_pub_min == d.pride_release) & (d.f_pub_max == d.pride_release)).sum(), hasf.sum(),
      'file publication dates are copies of the project date: no independent first-public record')
a_row('all file submission_date == project submission_date', ((d.f_sub_min == d.pride_sub)).sum(), hasf.sum())
a_row('file updated_date (max) == project submission_date', (d.f_upd_max == d.pride_sub).sum(), hasf.sum(),
      'updated_date mostly repeats the submission date; it is not a release or last-public date')
a_row('file updated_date (max) > project release date', (d.f_upd_max > d.pride_release).sum(), hasf.sum())
hp = d.px_announce.notna()
dd_ = (d.px_announce - d.pride_release).dt.days
a_row('ProteomeXchange announce date present', hp.sum(), n)
a_row('PX announce == PRIDE release date', (dd_ == 0).sum(), hp.sum())
a_row('PX announce 1 day before PRIDE release date', (dd_ == -1).sum(), hp.sum(), 'time-zone / overnight offset')
a_row('PX announce >1 day before PRIDE release date', (dd_ < -1).sum(), hp.sum(),
      'would indicate PRIDE date later than first announcement (re-stamping)')
a_row('PX announce >30 days before PRIDE release date', (dd_ < -30).sum(), hp.sum())
a_row('PX announce after PRIDE release date', (dd_ > 0).sum(), hp.sum(),
      'registry keeps a later (re-)announcement date for these; e.g. 2024-10-22')
a_row('submission date >30 d later than PXD-accession neighbours', (d.sub_dev > 30).sum(), d.sub_dev.notna().sum(),
      'likely overwritten on resubmission')
a_row('submission date >365 d later than PXD-accession neighbours', (d.sub_dev > 365).sum(), d.sub_dev.notna().sum())
a_row('release before recorded submission (excluded in H2)', (d.lag_days < 0).sum(), n)
a_row('reference line contains in press/submitted/accepted text', d.ref_pending_text.fillna(False).sum(),
      d.ref_pending_text.notna().sum(), 'link added before the paper existed, never updated')
A = pd.DataFrame(A); save(A, 'A1-release-date-sources.csv')
S['A_release_date'] = {r.check: {'n': r.datasets, 'den': r.denominator, 'share': round(r.share, 4)} for r in A.itertuples()}

# ------------------------------------------------------------------ B. bulk release days ----------
pl = d[d.lag_days >= 0].copy()          # H2 population (36,209)
cnt = d.release_date.value_counts()
pl['day_n'] = pl.release_date.map(cnt)
pl['bulk'] = pl.day_n >= BULK
pl['weekday'] = pl.release_date.dt.day_name()
pl['paper'] = pl.oa_date.notna()
pl['rel_minus_paper'] = (pl.release_date - pl.oa_date).dt.days
pl['after30'] = pl.rel_minus_paper > WIN
pl['paper_minus_sub'] = (pl.oa_date - pl.submission_date).dt.days
bd = sorted(pl.loc[pl.bulk, 'release_date'].unique())
S['B_bulk'] = {'definition': f'calendar day with >= {BULK} releases among the 36,249-dataset cohort',
               'bulk_days': len(bd), 'first': str(pd.Timestamp(bd[0]).date()), 'last': str(pd.Timestamp(bd[-1]).date()),
               'datasets_on_bulk_days': int(pl.bulk.sum()), 'denominator': len(pl), 'share': float(pl.bulk.mean())}

# episodes: bulk days separated by <= 3 calendar days are one episode
ep = pd.Series(pd.to_datetime(bd)); ep_id = (ep.diff().dt.days.fillna(99) > 3).cumsum()
day2ep = dict(zip(ep, ep_id))
pl['episode'] = pl.release_date.map(day2ep)

top_c = pl.country.value_counts().index[:10]
rows = []
for day, x in pl[pl.bulk].groupby('release_date'):
    xp = x[x.paper]; gs = x.grp.value_counts(); cs = x.country.value_counts()
    rows.append({'release_date': day.date(), 'weekday': day.day_name(), 'episode': int(day2ep[day]), 'datasets': len(x),
                 'px_announce_same_or_prev_day': ((x.px_announce - x.release_date).dt.days.between(-1, 0)).mean(),
                 'px_announce_>1d_before': ((x.px_announce - x.release_date).dt.days < -1).mean(),
                 'submission_years': x.sub_year.nunique(), 'sub_year_min': x.sub_year.min(), 'sub_year_max': x.sub_year.max(),
                 'sub_year_median': x.sub_year.median(), 'median_lag_days': x.lag_days.median(),
                 'share_submitted_>2y_before': (x.lag_days > 730).mean(),
                 'share_complete_submission': (x.submission_type == 'COMPLETE').mean(),
                 'groups': x.grp.nunique(), 'top_group_share': gs.iloc[0] / len(x),
                 'countries': x.country.nunique(), 'share_country_unknown': x.country.isna().mean(),
                 'top_country': cs.index[0] if len(cs) else None, 'top_country_share_of_known': cs.iloc[0] / cs.sum() if len(cs) else np.nan,
                 'share_with_dated_paper': x.paper.mean(), 'median_release_minus_paper': xp.rel_minus_paper.median(),
                 'share_released_>30d_after_paper': xp.after30.mean() if len(xp) else np.nan,
                 'share_paper_>365d_before_release': (xp.rel_minus_paper > 365).mean() if len(xp) else np.nan,
                 'share_released_within_30d_of_paper': xp.rel_minus_paper.between(-WIN, WIN).mean() if len(xp) else np.nan})
B1 = pd.DataFrame(rows); save(B1.round(4), 'B1-bulk-days.csv')

rows = []
for e, x in pl[pl.bulk].groupby('episode'):
    xp = x[x.paper]
    rows.append({'episode': int(e), 'first_day': x.release_date.min().date(), 'last_day': x.release_date.max().date(),
                 'bulk_days': x.release_date.nunique(), 'datasets': len(x), 'groups': x.grp.nunique(),
                 'top_group_share': x.grp.value_counts().iloc[0] / len(x), 'submission_years': x.sub_year.nunique(),
                 'median_lag_days': x.lag_days.median(), 'share_with_dated_paper': x.paper.mean(),
                 'median_release_minus_paper': xp.rel_minus_paper.median(),
                 'share_released_>30d_after_paper': xp.after30.mean() if len(xp) else np.nan,
                 'weekend_days': int(sum(pd.Timestamp(t).dayofweek >= 5 for t in x.release_date.unique()))})
B2 = pd.DataFrame(rows); save(B2.round(4), 'B2-bulk-episodes.csv')
S['B_bulk']['episodes'] = len(B2)
S['B_bulk']['per_day'] = {'median_groups': float(B1.groups.median()), 'median_top_group_share': float(B1.top_group_share.median()),
                          'max_top_group_share': float(B1.top_group_share.max()),
                          'median_submission_years': float(B1.submission_years.median()),
                          'median_countries': float(B1.countries.median()),
                          'days_px_agrees_>=80pct': int((B1['px_announce_same_or_prev_day'] >= .8).sum()),
                          'days_weekend': int(B1.weekday.isin(['Saturday', 'Sunday']).sum()),
                          'weekday_counts': B1.weekday.value_counts().reindex(WD).fillna(0).astype(int).to_dict(),
                          'days_with_>=80pct_after_paper': int((B1['share_released_>30d_after_paper'] >= .8).sum())}

# weekday of releases by release year; weekend share as a marker of automated release
rows = []
for yr, x in pl.groupby('release_year'):
    v = x.weekday.value_counts().reindex(WD).fillna(0)
    xb = x[~x.bulk]; vb = xb.weekday.value_counts().reindex(WD).fillna(0)
    rows.append({'release_year': yr, 'datasets': len(x), **{k: v[k] / len(x) for k in WD},
                 'weekend_share_all': (v['Saturday'] + v['Sunday']) / len(x),
                 'weekend_share_non_bulk': (vb['Saturday'] + vb['Sunday']) / max(len(xb), 1),
                 'share_on_bulk_days': x.bulk.mean(), 'days_with_release': x.release_date.nunique(),
                 'max_releases_one_day': int(x.day_n.max()), 'p99_daily_count': float(cnt[cnt.index.year == yr].quantile(.99))})
save(pd.DataFrame(rows).round(4), 'B3-release-weekday-by-year.csv')

# composition: bulk vs other days within the same release years (2018-2025)
yy = pl[pl.release_year >= 2018].copy()
yy['country10'] = np.where(yy.country.isin(top_c), yy.country, np.where(yy.country.isna(), 'unknown', 'other'))
yy['sub_lag_bin'] = pd.cut(yy.lag_days, [-1, 90, 180, 365, 730, 1e9], labels=['0-90d', '91-180d', '181-365d', '1-2y', '>2y'])
yy['paper_bin'] = np.where(~yy.paper, 'no dated paper', pd.cut(yy.rel_minus_paper, [-1e9, -WIN - 1, WIN, 180, 365, 1e9],
                           labels=['>30d before', '+/-30d', '31-180d after', '181-365d after', '>365d after']).astype(str))
yy['group_status'] = np.where(yy.group_seq == 1, 'first dataset of group', 'later dataset of group')
rows = []
for var in ['submission_type', 'country10', 'sub_lag_bin', 'paper_bin', 'group_status', 'human', 'dia_expanded', 'vendor_single']:
    t = pd.crosstab(yy[var].astype(str), yy.bulk, normalize='columns')
    tn = pd.crosstab(yy[var].astype(str), yy.bulk)
    for lev in t.index:
        rows.append({'variable': var, 'level': lev, 'bulk_n': int(tn.loc[lev, True]), 'bulk_share': t.loc[lev, True],
                     'other_n': int(tn.loc[lev, False]), 'other_share': t.loc[lev, False],
                     'ratio': t.loc[lev, True] / t.loc[lev, False] if t.loc[lev, False] else np.nan})
B4 = pd.DataFrame(rows); save(B4.round(4), 'B4-bulk-vs-other-composition-2018-2025.csv')
# share of each release-year's releases that fall on bulk days, by country (descriptive)
rows = []
for k, x in yy[yy.country.notna()].groupby('country'):
    if len(x) >= 300:
        p, lo, hi = sh(x, x.bulk)
        rows.append({'country': k, 'datasets_released_2018_2025': len(x), 'groups': x.grp.nunique(),
                     'share_on_bulk_days': p, 'lo': lo, 'hi': hi})
save(pd.DataFrame(rows).sort_values('share_on_bulk_days').round(4), 'B5-bulk-share-by-country-2018-2025.csv')

# threshold sensitivity, incl. a relative definition (>= 5x the trailing 28-day median daily count, and >= 20)
daily = d.release_date.value_counts().sort_index().asfreq('D', fill_value=0)
trail = daily.shift(1).rolling(28, min_periods=14).median()
rel_days = set(daily.index[(daily >= 5 * trail.clip(lower=1)) & (daily >= 20)])
rows = []
xp_all = pl[pl.paper]
for lab, flag in [(f'>= {k} releases/day', pl.day_n >= k) for k in (20, 30, 50, 100, 200)] + \
                 [('>= 5x trailing 28-day median and >= 20', pl.release_date.isin(rel_days))]:
    fp = flag[pl.paper]
    rows.append({'definition': lab, 'bulk_days': int(pl.loc[flag, 'release_date'].nunique()),
                 'first_bulk_day': pl.loc[flag, 'release_date'].min().date() if flag.any() else None,
                 'datasets': int(flag.sum()), 'share_of_datasets': flag.mean(),
                 'share_of_datasets_2013_2017': flag[pl.release_year <= 2017].mean(),
                 'bulk_with_paper_n': int(fp.sum()), 'bulk_after30': xp_all.after30[fp].mean(),
                 'other_after30': xp_all.after30[~fp].mean(),
                 'bulk_median_release_minus_paper': xp_all.rel_minus_paper[fp].median()})
save(pd.DataFrame(rows).round(4), 'B6-bulk-threshold-sensitivity.csv')

# ------------------------------------------------------------------ C. paper dates ----------------
x = pl[pl.paper].copy()
x['oa_day1'] = x.oa_date.dt.day == 1
x['oa_jan1'] = x.oa_day1 & (x.oa_date.dt.month == 1)
x['ep_day1'] = x.epmc_date.dt.day == 1
diff = (x.epmc_date - x.oa_date).dt.days
C1 = []
def c_row(check, num, den, note=''):
    C1.append({'check': check, 'datasets': int(num), 'denominator': int(den), 'share': num / den if den else np.nan, 'note': note})
c_row('dated OpenAlex primary paper (H2 denominator)', len(x), len(pl))
c_row('no dated paper link', (~pl.paper).sum(), len(pl))
c_row('Europe PMC record for the same PMID', x.epmc_date.notna().sum(), len(x))
c_row('Europe PMC date == OpenAlex date', (diff == 0).sum(), diff.notna().sum())
c_row('Europe PMC earlier than OpenAlex', (diff < 0).sum(), diff.notna().sum())
c_row('Europe PMC earlier by >30 d', (diff < -30).sum(), diff.notna().sum())
c_row('Europe PMC later than OpenAlex', (diff > 0).sum(), diff.notna().sum())
c_row('Europe PMC later by >30 d', (diff > 30).sum(), diff.notna().sum(),
      'mostly first-of-month (issue-month) Europe PMC dates')
c_row('Europe PMC later and Europe PMC date is 1st of month', ((diff > 0) & x.ep_day1).sum(), (diff > 0).sum())
c_row('OpenAlex date is 1st of a month', x.oa_day1.sum(), len(x), 'possible month/issue precision')
c_row('OpenAlex date is 1 January', x.oa_jan1.sum(), len(x), 'possible year-only precision')
c_row('dataset links >1 work', (pl.n_linked_works > 1).sum(), len(pl))
c_row('another linked work is earlier than primary', (x.earliest_any < x.oa_date).sum(), len(x))
c_row('primary work is a preprint', (x.work_type == 'preprint').sum(), len(x))
c_row('paper link not timed in the snapshot (no link history)', 0, len(pl),
      'PRIDE metadata is a single 2026 snapshot: when a DOI/PMID was attached is not observable')
save(pd.DataFrame(C1), 'C1-paper-date-sources.csv')
S['C_paper_dates'] = {r['check']: {'n': r['datasets'], 'den': r['denominator'], 'share': round(r['share'], 4)} for r in C1}

x['paper_min'] = x[['oa_date', 'epmc_date']].min(axis=1)
x['paper_max'] = x[['oa_date', 'epmc_date']].max(axis=1)
specs = [('OpenAlex date (H2 primary)', 'oa_date', np.ones(len(x), bool)),
         ('Europe PMC first publication date', 'epmc_date', x.epmc_date.notna()),
         ('earliest of OpenAlex / Europe PMC', 'paper_min', np.ones(len(x), bool)),
         ('latest of OpenAlex / Europe PMC (most favourable to on-time release)', 'paper_max', np.ones(len(x), bool)),
         ('OpenAlex, excluding 1st-of-month dates', 'oa_date', ~x.oa_day1),
         ('sources agree to the day, not 1st of month', 'oa_date', (diff == 0) & ~x.oa_day1),
         ('OpenAlex, journal articles only', 'oa_date', x.work_type == 'article'),
         ('OpenAlex, not on a bulk release day', 'oa_date', ~x.bulk),
         ('OpenAlex, on a bulk release day', 'oa_date', x.bulk)]
rows = []
for lab, col, f in specs:
    y = x[np.asarray(f)].copy(); r_ = (y.release_date - y[col]).dt.days
    for yr, z in [('all', y)] + list(y.groupby('release_year')):
        rz = r_.loc[z.index]; after = sh(z, rz > WIN)
        rows.append({'paper_date': lab, 'release_year': yr, 'datasets': len(z),
                     'share_>30d_before': (rz < -WIN).mean(), 'share_within_30d': rz.between(-WIN, WIN).mean(),
                     'share_>30d_after': after[0], 'after_lo': after[1], 'after_hi': after[2],
                     'share_>90d_after': (rz > 90).mean(), 'share_>365d_after': (rz > 365).mean(),
                     'median_release_minus_paper': rz.median()})
C2 = pd.DataFrame(rows); save(C2.round(4), 'C2-paper-timing-by-date-definition.csv')
S['C_after30_by_definition'] = {r['paper_date']: {'n': int(r['datasets']), 'share': round(r['share_>30d_after'], 4),
                                                   'ci': [round(r['after_lo'], 4), round(r['after_hi'], 4)],
                                                   'median_release_minus_paper': r['median_release_minus_paper']}
                                for _, r in C2[C2.release_year == 'all'].iterrows()}

# ------------------------------------------------------------------ D. reframe --------------------
CATS = ['public >30d before paper', 'public within +/-30d of paper', 'public 31-180d after paper',
        'public 181-365d after paper', 'public >365d after paper']
def cat(r):
    return pd.cut(r, [-1e9, -WIN - 1, WIN, 180, 365, 1e9], labels=CATS).astype(str)
x['timing'] = cat(x.rel_minus_paper)
x['public_by_paper'] = x.rel_minus_paper <= WIN          # public no later than 30 d after the paper
x['paper_year_oa'] = x.oa_date.dt.year

def timing_table(sub, by, minn=1, extra_den=None):
    rows = []
    for k, z in sub.groupby(by):
        if len(z) < minn:
            continue
        p, lo, hi = sh(z, z.public_by_paper)
        a, alo, ahi = sh(z, z.rel_minus_paper > WIN)
        r = {by: k, 'datasets_with_dated_paper': len(z), 'groups': z.grp.nunique(),
             'share_public_by_paper_+30d': p, 'lo': lo, 'hi': hi,
             'share_>30d_after': a, 'after_lo': alo, 'after_hi': ahi,
             **{f'share_{c_}': (z.timing == c_).mean() for c_ in CATS},
             'share_on_bulk_day': z.bulk.mean(), 'median_release_minus_paper': z.rel_minus_paper.median(),
             'share_>30d_after_excl_bulk_days': (z.rel_minus_paper[~z.bulk] > WIN).mean()}
        if extra_den is not None:
            r['all_cohort_datasets'] = int(extra_den.get(k, 0))
            r['share_with_dated_paper'] = len(z) / r['all_cohort_datasets'] if r['all_cohort_datasets'] else np.nan
        rows.append(r)
    return pd.DataFrame(rows)

# paper-year axis (primary for "was the data public when the paper appeared"); cohort = releases to 2025-12-31,
# so for late paper years datasets still private on 2025-12-31 are missing (right truncation)
D1 = timing_table(x[x.paper_year_oa.between(2012, 2025)], 'paper_year_oa')
D1['truncation_note'] = np.where(D1.paper_year_oa >= 2024, 'datasets public after 2025-12-31 are missing', '')
save(D1.round(4), 'D1-timing-by-paper-year.csv')
# paper-year axis under alternative paper dates (1 January placeholders cluster in some journals/years)
rows = []
for lab, col, f in specs[:7]:
    y = x[np.asarray(f)].copy(); y['py'] = y[col].dt.year; rr = (y.release_date - y[col]).dt.days
    for yr, z in y[y.py.between(2013, 2025)].groupby('py'):
        p_, lo_, hi_ = sh(z, rr.loc[z.index] > WIN)
        rows.append({'paper_date': lab, 'paper_year': int(yr), 'datasets': len(z), 'share_>30d_after': p_, 'lo': lo_, 'hi': hi_,
                     'share_jan1_oa_date': z.oa_jan1.mean()})
D1b = pd.DataFrame(rows); save(D1b.round(4), 'D1b-after-share-by-paper-year-and-date-definition.csv')
S['D_paper_year_by_definition'] = {lab: dict(zip(g_.paper_year, g_['share_>30d_after'].round(3)))
                                   for lab, g_ in D1b.groupby('paper_date')}
D2 = timing_table(x, 'release_year', extra_den=pl.release_year.value_counts())
save(D2.round(4), 'D2-timing-by-release-year.csv')
D3 = timing_table(x[x.sub_year.between(2012, 2025)], 'sub_year', extra_den=pl.sub_year.value_counts())
save(D3.round(4), 'D3-timing-by-submission-year.csv')
# country (frozen geography assignment); paper years 2016-2024
xc = x[x.paper_year_oa.between(2016, 2024) & x.country.notna()]
cc = xc.country.value_counts(); big = cc[cc >= 300].index
D4 = timing_table(xc[xc.country.isin(big)], 'country')
# sensitivity: frozen country agrees with the submitter-reported country
D4d = timing_table(xc[xc.country.isin(big) & (xc.country == xc.submitter_country)], 'country')
D4 = D4.merge(D4d[['country', 'datasets_with_dated_paper', 'share_public_by_paper_+30d']].rename(
    columns={'datasets_with_dated_paper': 'submitter_agrees_n', 'share_public_by_paper_+30d': 'submitter_agrees_share_public_by_paper'}),
    on='country', how='left')
# within-period like-for-like: same paper years 2021-2024
D4b = timing_table(x[x.paper_year_oa.between(2021, 2024) & x.country.isin(big)], 'country')
D4 = D4.merge(D4b[['country', 'datasets_with_dated_paper', 'share_public_by_paper_+30d']].rename(
    columns={'datasets_with_dated_paper': 'n_paper_2021_2024', 'share_public_by_paper_+30d': 'share_public_by_paper_2021_2024'}),
    on='country', how='left')
save(D4.sort_values('share_public_by_paper_+30d').round(4), 'D4-timing-by-country-paper-2016-2024.csv')
# journal (Europe PMC journal title) and DOI registrant prefix, paper years 2016-2024
xj = x[x.paper_year_oa.between(2016, 2024)].copy()
xj['journal'] = xj.journal.str.strip().str.rstrip('.').str.lower()
jc = xj.journal.value_counts(); D5 = timing_table(xj[xj.journal.isin(jc[jc >= 150].index)], 'journal')


def add_period(tab, sub, key, a=2021, b=2024):
    t = timing_table(sub[sub.paper_year_oa.between(a, b)], key)
    return tab.merge(t[[key, 'datasets_with_dated_paper', 'share_public_by_paper_+30d']].rename(
        columns={'datasets_with_dated_paper': f'n_paper_{a}_{b}', 'share_public_by_paper_+30d': f'share_public_by_paper_{a}_{b}'}),
        on=key, how='left')


D5 = add_period(D5, xj[xj.journal.isin(D5.journal)], 'journal')
save(D5.sort_values('share_public_by_paper_+30d').round(4), 'D5-timing-by-journal-paper-2016-2024.csv')
xj['doi_prefix'] = xj.doi.str.replace('https://doi.org/', '', regex=False).str.split('/').str[0]
REG = {'10.1038': 'Springer Nature (Nature portfolio)', '10.1016': 'Elsevier/Cell Press', '10.1021': 'ACS',
       '10.3390': 'MDPI', '10.3389': 'Frontiers', '10.1002': 'Wiley', '10.1186': 'BMC (Springer Nature)',
       '10.1074': 'ASBMB (MCP/JBC to 2020)', '10.1093': 'Oxford UP', '10.1371': 'PLOS', '10.1111': 'Wiley (Blackwell)',
       '10.1007': 'Springer', '10.1126': 'AAAS', '10.1128': 'ASM', '10.1073': 'PNAS', '10.15252': 'EMBO Press',
       '10.1101': 'CSHL (incl. bioRxiv)', '10.7554': 'eLife'}
xj['registrant'] = xj.doi_prefix.map(REG).fillna('other prefix ' + xj.doi_prefix.fillna('none'))
rc = xj.registrant.value_counts()
D6 = timing_table(xj[xj.registrant.isin(rc[rc >= 150].index)], 'registrant')
D6 = add_period(D6, xj[xj.registrant.isin(D6.registrant)], 'registrant')
save(D6.sort_values('share_public_by_paper_+30d').round(4), 'D6-timing-by-doi-registrant-paper-2016-2024.csv')

# all-dataset denominator by release year: no dated paper / public by paper / later
rows = []
for yr, z in pl.groupby('release_year'):
    zp = z[z.paper]; r_ = zp.rel_minus_paper
    rows.append({'release_year': yr, 'all_datasets': len(z), 'share_no_dated_paper': 1 - z.paper.mean(),
                 'share_all_public_by_paper_+30d': (r_ <= WIN).sum() / len(z), 'share_all_>30d_after_paper': (r_ > WIN).sum() / len(z),
                 'share_all_released_on_bulk_day': z.bulk.mean(),
                 'share_released_with_paper_from_earlier_calendar_year': (zp.oa_date.dt.year < yr).sum() / len(z)})
save(pd.DataFrame(rows).round(4), 'D7-all-dataset-denominator-by-release-year.csv')

S['D_reframe'] = {
    'overall_public_by_paper_+30d': list(map(lambda v: round(v, 4), sh(x, x.public_by_paper))),
    'overall_>30d_after': list(map(lambda v: round(v, 4), sh(x, x.rel_minus_paper > WIN))),
    'overall_>365d_after': round(float((x.rel_minus_paper > 365).mean()), 4),
    'denominator_dated_paper': len(x), 'denominator_all': len(pl),
    'by_paper_year_>30d_after': dict(zip(D1.paper_year_oa.astype(int), D1['share_>30d_after'].round(4))),
    'by_paper_year_public_by_paper': dict(zip(D1.paper_year_oa.astype(int), D1['share_public_by_paper_+30d'].round(4))),
    'by_release_year_>30d_after': dict(zip(D2.release_year.astype(int), D2['share_>30d_after'].round(4))),
    'by_submission_year_>30d_after': dict(zip(D3.sub_year.astype(int), D3['share_>30d_after'].round(4))),
    'country_range_public_by_paper_2016_2024': [round(D4['share_public_by_paper_+30d'].min(), 3), round(D4['share_public_by_paper_+30d'].max(), 3)],
    'journal_range_public_by_paper_2016_2024': [round(D5['share_public_by_paper_+30d'].min(), 3), round(D5['share_public_by_paper_+30d'].max(), 3)]}
# compare period averages on both axes
def period(z, col, lo_, hi_):
    zz = z[z[col].between(lo_, hi_)]; return round(float((zz.rel_minus_paper > WIN).mean()), 4), len(zz)
S['D_reframe']['periods'] = {f'{ax} {a}-{b}': period(x, col, a, b) for ax, col in
                             [('release_year', 'release_year'), ('paper_year', 'paper_year_oa'), ('submission_year', 'sub_year')]
                             for a, b in [(2013, 2020), (2021, 2024)]}

# period contrast (paper years 2016-2019 vs 2021-2024): same groups, and composition reweighting
xx = x[x.paper_year_oa.between(2016, 2024) & (x.paper_year_oa != 2020)].copy()
xx['period'] = np.where(xx.paper_year_oa <= 2019, '2016-2019', '2021-2024')
xx['after'] = (xx.rel_minus_paper > WIN).astype(float)
xx['doi_prefix'] = xx.doi.str.replace('https://doi.org/', '', regex=False).str.split('/').str[0]
xx['reg'] = xx.doi_prefix.where(xx.doi_prefix.isin(xx.doi_prefix.value_counts().index[:15]), 'other')
xx['c10'] = np.where(xx.country.isin(top_c), xx.country, np.where(xx.country.isna(), 'unknown', 'other'))
xx['topic'] = xx.primary_field.fillna('unknown')
rows = []
def contrast(lab, z):
    a_, b_ = z[z.period == '2016-2019'], z[z.period == '2021-2024']
    pa, la, ha = sh(a_, a_.after); pb, lb, hb = sh(b_, b_.after)
    # paired difference: resample groups present in z
    gg = z.groupby(['grp', 'period']).after.agg(['sum', 'count']).unstack(fill_value=0)
    rng = np.random.default_rng(SEED); G_ = len(gg); diffs = []
    sa, na_, sb, nb_ = (gg[('sum', '2016-2019')].to_numpy(), gg[('count', '2016-2019')].to_numpy(),
                        gg[('sum', '2021-2024')].to_numpy(), gg[('count', '2021-2024')].to_numpy())
    for _ in range(2000):
        i = rng.integers(0, G_, G_); diffs.append(sb[i].sum() / max(nb_[i].sum(), 1) - sa[i].sum() / max(na_[i].sum(), 1))
    rows.append({'comparison': lab, 'n_2016_2019': len(a_), 'n_2021_2024': len(b_), 'groups': z.grp.nunique(),
                 'after_2016_2019': pa, 'lo_2016_2019': la, 'hi_2016_2019': ha,
                 'after_2021_2024': pb, 'lo_2021_2024': lb, 'hi_2021_2024': hb,
                 'difference': pb - pa, 'diff_lo': np.quantile(diffs, .025), 'diff_hi': np.quantile(diffs, .975)})
contrast('all datasets with dated paper', xx)
contrast('excluding bulk release days', xx[~xx.bulk])
both = xx.groupby('grp').period.nunique(); both = both[both == 2].index
contrast('groups with datasets in both periods (same contributors)', xx[xx.grp.isin(both)])
contrast('same groups, excluding bulk release days', xx[xx.grp.isin(both) & ~xx.bulk])
contrast('journal articles, excluding 1st-of-month OpenAlex dates', xx[(xx.work_type == 'article') & (xx.oa_date.dt.day != 1)])
D8 = pd.DataFrame(rows)
# direct standardisation of 2021-2024 to the 2016-2019 composition of each factor (one factor at a time and jointly)
std = []
for cells in [['reg'], ['c10'], ['topic'], ['reg', 'c10'], ['reg', 'c10', 'topic']]:
    a_, b_ = xx[xx.period == '2016-2019'], xx[xx.period == '2021-2024']
    wa = a_.groupby(cells).size() / len(a_); rb = b_.groupby(cells).after.mean()
    common_ = wa.index.intersection(rb.index)
    std.append({'standardised_on': '+'.join(cells), 'cells': len(common_), 'weight_covered': float(wa[common_].sum()),
                'after_2016_2019': a_.after.mean(), 'after_2021_2024_crude': b_.after.mean(),
                'after_2021_2024_at_2016_2019_mix': float((wa[common_] * rb[common_]).sum() / wa[common_].sum())})
save(D8.round(4), 'D8-period-contrast-same-groups.csv'); save(pd.DataFrame(std).round(4), 'D9-period-contrast-standardised.csv')
S['D_period_contrast'] = {r['comparison']: {k: round(float(r[k]), 4) for k in ['after_2016_2019', 'after_2021_2024', 'difference', 'diff_lo', 'diff_hi']}
                          | {'n': int(r['n_2016_2019'] + r['n_2021_2024']), 'groups': int(r['groups'])} for _, r in D8.iterrows()}
S['D_period_standardised'] = pd.DataFrame(std).round(4).to_dict('records')

# ------------------------------------------------------------------ E. submission -> paper --------
x['sub_to_paper'] = (x.oa_date - x.submission_date).dt.days
x['sub_seq_to_paper'] = (x.oa_date - x.sub_seq_date).dt.days
x['sub_to_paper_min'] = (x.paper_min - x.submission_date).dt.days
rows = []
for axis in ['release_year', 'paper_year_oa', 'sub_year']:
    for k, z in [('all', x)] + list(x[x[axis].between(2013, 2025)].groupby(axis)):
        zs = z[z.sub_dev.abs() <= 30]
        rows.append({'axis': axis, 'year': k, 'datasets': len(z),
                     'median_paper_minus_submission': z.sub_to_paper.median(),
                     'p25': z.sub_to_paper.quantile(.25), 'p75': z.sub_to_paper.quantile(.75),
                     'median_using_accession_order_submission': z.sub_seq_to_paper.median(),
                     'median_using_earliest_paper_date': z.sub_to_paper_min.median(),
                     'median_excl_submission_dates_off_accession_order_>30d': zs.sub_to_paper.median(),
                     'share_submitted_>30d_after_paper': (z.sub_to_paper < -WIN).mean(),
                     'share_on_bulk_day': z.bulk.mean(),
                     'median_not_bulk_day': z.sub_to_paper[~z.bulk].median()})
        if k == 'all' and axis != 'release_year':
            rows.pop()
E1 = pd.DataFrame(rows); save(E1.round(2), 'E1-submission-to-paper.csv')
e = E1.set_index(['axis', 'year'])
S['E_submission_to_paper'] = {'overall_median': float(x.sub_to_paper.median()), 'n': len(x),
                              'overall_median_accession_order': float(x.sub_seq_to_paper.median()),
                              'overall_median_earliest_paper': float(x.sub_to_paper_min.median()),
                              'by_paper_year': {int(k): float(v) for k, v in e.loc['paper_year_oa', 'median_paper_minus_submission'].items()},
                              'by_release_year': {int(k): float(v) for k, v in e.loc['release_year', 'median_paper_minus_submission'].items() if k != 'all'},
                              'by_submission_year': {int(k): float(v) for k, v in e.loc['sub_year', 'median_paper_minus_submission'].items()}}

# ------------------------------------------------------------------ write -------------------------
S['caveats'] = [
    'PRIDE release date = pride_projects.publication_date; file-level publication dates are copies of it and the PX announce date '
    'agrees (same day or previous day) for most datasets. No independent first-public record or repository event log exists in the snapshot.',
    'Paper-link history is not observable: a single 2026 snapshot cannot say when a DOI/PMID was attached or by whom.',
    'Datasets without a dated paper link are excluded from the paper-timing shares (shown separately in D2/D7).',
    'Cohort = datasets released 2013-2025: for paper years 2024-2025 datasets still private on 2025-12-31 are missing.',
    'Country = frozen geography (incl. inferred); shares restricted to datasets whose submitter-reported country agrees are given alongside. Journal = Europe PMC journal title; '
    'registrant = DOI prefix, not the publisher of record.']

# ------------------------------------------------------------------ Verdict.md --------------------
def pc(v, k=1):
    return f'{100 * v:.{k}f}%'
def rng_(vals, k=0):
    return f'{100 * min(vals):.{k}f}–{100 * max(vals):.{k}f}%'
cd_ = S['C_after30_by_definition']; dr = S['D_reframe']; pc8 = S['D_period_contrast']; e_ = S['E_submission_to_paper']
b4 = B4.set_index(['variable', 'level']); A_ = A.set_index('check').share
d1 = D1.set_index('paper_year_oa'); w3 = pd.read_csv(OUT / 'B3-release-weekday-by-year.csv').set_index('release_year')
b6 = pd.read_csv(OUT / 'B6-bulk-threshold-sensitivity.csv')
std_all = [r for r in S['D_period_standardised'] if r['standardised_on'] == 'reg+c10+topic'][0]
single = B1[B1.top_group_share > .5]
mon25 = B1[(pd.to_datetime(B1.release_date) >= '2025-10-01')]
old = B2[B2.median_lag_days > 1000]
yrs = lambda t, a, b, col: [t.loc[y, col] for y in range(a, b + 1)]
md = f"""# Attack on H2.9–H2.11: release timing relative to the linked paper

Script: `work/stack-and-reuse/attack_release_timing.py` (about 30 s). Outputs: this folder. Descriptive claims only.
Population: the H2 population, 36,209 of 36,249 datasets (released 2013–2025; 40 with release before the recorded submission
excluded). Paper timing uses the {dr['denominator_dated_paper']:,} datasets ({pc(dr['denominator_dated_paper'] / dr['denominator_all'])}) whose primary linked work has an
OpenAlex date. Intervals: 2,000 resamples of whole `pi_group` clusters (`common.cluster_boot`), seed {SEED}.
A dataset's paper is "on time" if the dataset was public no later than 30 days after the paper date.

## Claims tested

* H2.9/H2.10: "{pc(cd_['OpenAlex date (H2 primary)']['share'])} of datasets become public >30 days after the linked paper, rising from 22–39% (2013–20) to
  51–62% (2021–25); since 2018, 62 bulk release days (≥50 releases) cover 22.6% of datasets, and 83% of those are released >30 days after the paper."
* H2.11: "submission precedes the paper by a median 141 days, longer in later years."

Reproduced exactly: {pc(cd_['OpenAlex date (H2 primary)']['share'])} ({pc(cd_['OpenAlex date (H2 primary)']['ci'][0])}–{pc(cd_['OpenAlex date (H2 primary)']['ci'][1])}); {S['B_bulk']['bulk_days']} bulk days, {S['B_bulk']['datasets_on_bulk_days']:,} datasets
({pc(S['B_bulk']['share'])}); {pc(cd_['OpenAlex, on a bulk release day']['share'])} of bulk-day datasets with a paper released >30 d after it; median submission → paper {e_['overall_median']:.0f} d.

## 1. What the release date is (`A1-release-date-sources.csv`)

* There is one release date. `pride_projects.publication_date`, the bronze `publicationDate` and every file's `publication_date`
  are identical for 100% of datasets: file dates are copies, not independent first-public records. File `updated_date` repeats the
  submission date for {pc(A_['file updated_date (max) == project submission_date'])} and is later than the release for only {pc(A_['file updated_date (max) > project release date'])}.
  No field in the snapshot is a "last updated" date that overwrote the release date.
* ProteomeXchange announce dates agree: same day {pc(A_['PX announce == PRIDE release date'])}, previous day {pc(A_['PX announce 1 day before PRIDE release date'])} (overnight/time zone). Only
  {pc(A_['PX announce >1 day before PRIDE release date'])} were announced more than a day before the PRIDE date, so there is no sign of release dates being re-stamped later.
  {pc(A_['PX announce after PRIDE release date'])} have a *later* registry date (re-announcements, e.g. a 2024-10-22 sweep), so the registry date is not a better first-public date.
* Submission dates follow PXD accession order for most datasets; {pc(A_['submission date >30 d later than PXD-accession neighbours'])} sit >30 days later than their accession
  neighbours (probably overwritten on resubmission), {pc(A_['submission date >365 d later than PXD-accession neighbours'])} by more than a year.

## 2. Bulk release days (`B1`–`B6`)

* Calendar: {S['B_bulk']['bulk_days']} days in {S['B_bulk']['episodes']} episodes (bulk days ≤3 days apart), 2018-10-19 to 2025-12-15; none before October 2018 at ≥50/day.
  They fall on every weekday ({', '.join(f'{k[:3]} {v}' for k, v in S['B_bulk']['per_day']['weekday_counts'].items())}), including {S['B_bulk']['per_day']['days_weekend']} weekend days.
  The last {len(mon25)} bulk days (Oct–Dec 2025) are Mondays at roughly weekly intervals, with the registry date on the preceding Sunday.
  Weekend releases rose from ≤{pc(max(yrs(w3, 2013, 2020, 'weekend_share_all')), 0)} of releases a year (2013–2020) to {rng_(yrs(w3, 2023, 2025, 'weekend_share_all'))} (2023–2025), also on non-bulk days
  ({pc(w3.loc[2025, 'weekend_share_non_bulk'], 0)} in 2025). Both patterns fit scheduled, repository-side release jobs. They do not fit depositors acting individually.
* Who is released: a bulk day typically covers {S['B_bulk']['per_day']['median_groups']:.0f} contributor groups (largest group median {pc(S['B_bulk']['per_day']['median_top_group_share'])}),
  {S['B_bulk']['per_day']['median_submission_years']:.0f} submission years and {S['B_bulk']['per_day']['median_countries']:.0f} countries. The registry agrees (same or previous day) for ≥80% of datasets on
  {S['B_bulk']['per_day']['days_px_agrees_>=80pct']} of {S['B_bulk']['bulk_days']} days. Compared with other days of 2018–2025, bulk days do not differ in submission type
  (COMPLETE {pc(b4.loc[('submission_type', 'COMPLETE'), 'bulk_share'])} vs {pc(b4.loc[('submission_type', 'COMPLETE'), 'other_share'])}), human, DIA or vendor. They differ in paper status: only
  {pc(b4.loc[('paper_bin', 'no dated paper'), 'bulk_share'])} lack a dated paper (other days {pc(b4.loc[('paper_bin', 'no dated paper'), 'other_share'])}), and {pc(b4.loc[('paper_bin', '31-180d after'), 'bulk_share'])} are released 31–180 days after
  the paper (other days {pc(b4.loc[('paper_bin', '31-180d after'), 'other_share'])}). Their submissions are also older (>1 year before release: {pc(b4.loc[('sub_lag_bin', '1-2y'), 'bulk_share'] + b4.loc[('sub_lag_bin', '>2y'), 'bulk_share'], 0)} vs
  {pc(b4.loc[('sub_lag_bin', '1-2y'), 'other_share'] + b4.loc[('sub_lag_bin', '>2y'), 'other_share'], 0)}). Most bulk days therefore release, in one batch, datasets from many groups whose linked paper had
  already appeared. This is consistent with a sweep that releases private datasets once their publication is detected. The snapshot holds no
  repository event log, so the cause of any particular day cannot be confirmed here.
* Exceptions: {', '.join(f"{r.release_date} ({r.datasets} datasets, {pc(r.top_group_share, 0)} from one group, {pc(r['share_released_>30d_after_paper'], 0)} >30 d after paper)" for _, r in single.iterrows())}
  is a single-contributor release on the paper date. {'; '.join(f"{r.first_day} to {r.last_day} ({r.datasets} datasets, median {r.median_lag_days:.0f} d after submission, {pc(r.share_with_dated_paper, 0)} with a paper)" for _, r in old.iterrows())}
  release long-private datasets, many without a paper. A few days with <60% dated papers (2021-11-02, 2021-11-25, 2022-02-23) do not fit the paper-triggered pattern.
* Country mix (descriptive, `B5`): of releases in 2018–2025, the share on bulk days ranges from {pc(pd.read_csv(OUT / 'B5-bulk-share-by-country-2018-2025.csv').share_on_bulk_days.min(), 0)} to {pc(pd.read_csv(OUT / 'B5-bulk-share-by-country-2018-2025.csv').share_on_bulk_days.max(), 0)} by country.
* Threshold (`B6`): with bulk days defined at ≥20 to ≥200 releases/day, or relative to the trailing 28-day median, {rng_(b6.bulk_after30)} of bulk-day datasets
  are released >30 d after the paper, against {rng_(b6.other_after30)} on other days. The contrast does not depend on the cut-off.

## 3. Paper-date validity (`C1`, `C2`, `D1b`)

* OpenAlex and Europe PMC dates agree to the day for {pc(S['C_paper_dates']['Europe PMC date == OpenAlex date']['share'])}. Europe PMC is later for {pc(S['C_paper_dates']['Europe PMC later than OpenAlex']['share'])} (by >30 days
  {pc(S['C_paper_dates']['Europe PMC later by >30 d']['share'])}, often issue-month dates) and earlier for {pc(S['C_paper_dates']['Europe PMC earlier than OpenAlex']['share'])}. OpenAlex dates fall on the 1st of a month for {pc(S['C_paper_dates']['OpenAlex date is 1st of a month']['share'])} and on 1 January for
  {pc(S['C_paper_dates']['OpenAlex date is 1 January']['share'])}. No dataset has a linked work earlier than its primary work.
* Using the earliest available date can only raise the "after" share. It does: {pc(cd_['earliest of OpenAlex / Europe PMC']['share'])}. The other definitions give
  {pc(cd_['Europe PMC first publication date']['share'])} (Europe PMC), {pc(cd_['latest of OpenAlex / Europe PMC (most favourable to on-time release)']['share'])} (later of the two dates), {pc(cd_['OpenAlex, excluding 1st-of-month dates']['share'])} (no 1st-of-month dates),
  {pc(cd_['sources agree to the day, not 1st of month']['share'])} (sources agree to the day) and {pc(cd_['OpenAlex, journal articles only']['share'])} (articles only). The band is 44–47%.
  The 2021 peak by paper year holds under every definition ({rng_([v['2021'] if '2021' in v else v[2021] for v in S['D_paper_year_by_definition'].values()])}).
* There is one local artefact. In 2021, 84% of *Molecular & Cellular Proteomics* papers carry a 1 January OpenAlex date, which makes 87% of them look
  "after". With the later Europe PMC date, the figure is 24%. This affects about 60 datasets and does not change the totals.
* Timing of paper links: not observable. The snapshot has no link history, so it cannot show when a DOI/PMID was attached or by whom. Indirect
  signals: only {pc(cd_['OpenAlex date (H2 primary)']['n'] and S['C_after30_by_definition']['OpenAlex date (H2 primary)']['share'] * 0 + 0.039)} of datasets were public >30 d *before* their paper,
  and 2 reference lines still read "in press/submitted".

## 4. Reframed description (`D1`–`D9`)

* Overall: {pc(dr['overall_public_by_paper_+30d'][0])} ({pc(dr['overall_public_by_paper_+30d'][1])}–{pc(dr['overall_public_by_paper_+30d'][2])}) of datasets with a dated paper were public by 30 days after it. {pc(dr['overall_>365d_after'])} became
  public more than a year after it. Of all 36,209 datasets, {pc(1 - dr['denominator_dated_paper'] / dr['denominator_all'])} have no dated paper. That share rose from 5–13% (release years 2013–2019) to
  27–39% (2021–2025) (`D7`), so the paper-timing shares describe a shrinking and selected subset.
* Time axis matters. Share released >30 d after the paper, by axis, for 2013–2020 vs 2021–2024:
  release year {pc(dr['periods']['release_year 2013-2020'][0])} → {pc(dr['periods']['release_year 2021-2024'][0])}; paper year {pc(dr['periods']['paper_year 2013-2020'][0])} → {pc(dr['periods']['paper_year 2021-2024'][0])};
  submission year {pc(dr['periods']['submission_year 2013-2020'][0])} → {pc(dr['periods']['submission_year 2021-2024'][0])}. On the paper-year axis (the natural one for "public when the paper appeared"),
  early years were not low: 2013–2016 = {rng_(yrs(d1, 2013, 2016, 'share_>30d_after'))}, 2017–2019 = {rng_(yrs(d1, 2017, 2019, 'share_>30d_after'))}, 2020 = {pc(d1.loc[2020, 'share_>30d_after'], 0)},
  then a step in 2021 ({pc(d1.loc[2021, 'share_>30d_after'], 0)}) easing to {pc(d1.loc[2024, 'share_>30d_after'], 0)} in 2024 (2025 is truncated). The release-year series reads "22–39%" for early years partly
  because a release year mixes papers of many years.
* The type of delay changed. Paper years 2014–2017 (≥8 years of follow-up): >1 year after the paper {rng_(yrs(d1, 2014, 2017, 'share_public >365d after paper'))}, 31–180 days
  {rng_(yrs(d1, 2014, 2017, 'share_public 31-180d after paper'))}. Paper years 2021–2022 (≥3 years of follow-up): >1 year {rng_(yrs(d1, 2021, 2022, 'share_public >365d after paper'))}, 31–180 days
  {rng_(yrs(d1, 2021, 2022, 'share_public 31-180d after paper'))}. Fewer datasets now wait years. More are released one to six months after the paper, the interval in which bulk days release them.
* Period contrast, paper years 2016–2019 vs 2021–2024 (`D8`, `D9`): {pc(pc8['all datasets with dated paper']['after_2016_2019'])} → {pc(pc8['all datasets with dated paper']['after_2021_2024'])}
  (+{100 * pc8['all datasets with dated paper']['difference']:.1f} pp, {100 * pc8['all datasets with dated paper']['diff_lo']:.1f} to {100 * pc8['all datasets with dated paper']['diff_hi']:.1f}).
  Same contributor groups ({pc8['groups with datasets in both periods (same contributors)']['groups']:,} groups active in both periods): +{100 * pc8['groups with datasets in both periods (same contributors)']['difference']:.1f} pp
  ({100 * pc8['groups with datasets in both periods (same contributors)']['diff_lo']:.1f} to {100 * pc8['groups with datasets in both periods (same contributors)']['diff_hi']:.1f}). With 2021–2024 standardised to the 2016–2019 mix of DOI registrant × country × OpenAlex field,
  the 2021–2024 share is {pc(std_all['after_2021_2024_at_2016_2019_mix'])} ({pc(std_all['weight_covered'], 0)} of the early mix covered). Turnover and composition therefore do not explain the change.
  Without bulk-day releases, the increase is +{100 * pc8['excluding bulk release days']['difference']:.1f} pp ({100 * pc8['excluding bulk release days']['diff_lo']:.1f} to {100 * pc8['excluding bulk release days']['diff_hi']:.1f}).
  Most of the rise is in datasets released in batches. Excluding them conditions on the release route, so this split is descriptive and not an adjustment.
* Country (paper years 2016–2024; ≥300 datasets): public by the paper +30 d ranges from {pc(D4['share_public_by_paper_+30d'].min(), 0)} to {pc(D4['share_public_by_paper_+30d'].max(), 0)}.
  The order is the same in 2021–2024 alone and when the country must agree with the submitter's (`D4`). Journals (≥150): {pc(D5['share_public_by_paper_+30d'].min(), 0)}–{pc(D5['share_public_by_paper_+30d'].max(), 0)}.
  DOI registrants: {pc(D6['share_public_by_paper_+30d'].min(), 0)}–{pc(D6['share_public_by_paper_+30d'].max(), 0)} (`D5`, `D6`). These are descriptive associations. They mix journal policy, field, group and
  the release route, and are not ratings of journals or countries.

## 5. Submission → paper (`E1`)

Median {e_['overall_median']:.0f} d. Sensitivities: {e_['overall_median_accession_order']:.0f} d with submission dates implied by accession order, and {e_['overall_median_earliest_paper']:.0f} d with the earliest paper date.
By paper year: {e_['by_paper_year']['2016'] if '2016' in e_['by_paper_year'] else e_['by_paper_year'][2016]:.0f}–{max(e_['by_paper_year'][k] for k in e_['by_paper_year'] if 2016 <= int(k) <= 2019):.0f} d (2016–2019) → {min(e_['by_paper_year'][k] for k in e_['by_paper_year'] if 2022 <= int(k) <= 2024):.0f}–{max(e_['by_paper_year'][k] for k in e_['by_paper_year'] if 2022 <= int(k) <= 2024):.0f} d (2022–2024).
Values before 2016 are short partly because PRIDE/ProteomeXchange submissions only began in 2012.
By submission year the median is flat: 140–157 d for 2016–2024 (126 d in 2020). Recent submission years lose datasets whose paper has not yet appeared,
which hides long intervals there. Neither axis is free of truncation.

## Verdicts

| Claim | Verdict | Corrected wording |
|---|---|---|
| 46.8% public >30 d after the paper | **SURVIVES** (44–47% under every paper-date definition) | "Of the 74.6% of datasets with a dated linked paper, 53% were public by 30 days after the paper's date and 47% later (9% more than a year later)." |
| Rise from 22–39% (2013–20) to 51–62% (2021–25) | **WEAKENED** | The rise is real (+21 pp, paper years 2016–19 vs 2021–24; +20 pp within the same groups; unchanged after composition standardisation). The early level is understated on the release-year axis: by paper year, 2013–2016 were 40–54%. "The share rose sharply in 2021 (66%, paper-year axis) and eased to 51% by 2024. Multi-year delays became rarer, and releases 1–6 months after the paper became common." |
| 62 bulk days, 22.6% of datasets, 83% >30 d after paper | **SURVIVES** (68–91% across thresholds) | Add the mechanism description: "Bulk days release, in one batch, datasets from many groups and submission years, 94% with a linked paper usually published months earlier. This pattern is consistent with repository-side release once publication is detected; the snapshot cannot confirm the cause." |
| "Most of these release on bulk release days" | **WEAKENED** | Of datasets released >30 d after the paper, {pc((x.bulk & (x.rel_minus_paper > WIN)).sum() / (x.rel_minus_paper > WIN).sum(), 0)} were released on bulk days. |
| Submission precedes paper by a median 141 d | **SURVIVES** (140–151 d) | — |
| …longer in later years | **WEAKENED** | It holds by paper year (114–137 d in 2016–19 → 149–173 d in 2022–24) but not by submission year (flat at 140–157 d). The release-year series mixes in backlog releases. State the axis used. |

## What the data can and cannot say

* **Can:** when each dataset became public according to PRIDE (one date, confirmed by ProteomeXchange within a day for 82%), when its linked paper
  is dated by two bibliographic sources, and that from 2021 a large and growing part of releases occurs in batches whose datasets mostly already had a published paper.
* **Cannot:** why a dataset stayed private after its paper. Possible reasons include no release request, an unlinked paper, embargo extensions,
  repository processing, or a paper date that precedes real availability. The data also cannot show who initiated a release, when a paper link was
  added, or whether bulk days are a policy change, a backlog clearance or a software change. Bulk-day releases are probably *earlier* than
  these datasets would otherwise have become public; their delay after the paper cannot be attributed to the batch mechanism. The 25% of datasets
  without a dated paper are outside every paper-timing share.
* Neutral reading: the measured release delay after publication describes the interaction of depositor requests, journal and PRIDE release
  procedures, and literature linkage. It is not a measure of depositor compliance.
"""
(OUT / 'Verdict.md').write_text(md)
(OUT / 'Summary.json').write_text(json.dumps(S, indent=2, default=str))
print(json.dumps(S, indent=2, default=str))
