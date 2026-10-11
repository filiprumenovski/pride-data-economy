"""Arsenal H4 (sequence databases) and H7 (storage anatomy and duplication). Descriptive only.

Run from work/stack-and-reuse/. Reads pride_files (top-level, non-generated files), the members of
listed results zips (Archive-contents/Members.parquet, Archives.csv.gz) and deposit features.
Writes only to outputs/Stack-and-reuse/Arsenal/H4-H7-databases-storage/.

H4  Sequence databases (FASTA)
  * FASTA = file name ending .fasta/.fa/.faa/.fas (one compression suffix allowed), any PRIDE file
    category. This extends the layer rule `sequence_db` (.fasta/.fa/.faa, and not filed as RAW) by
    `.fas`, which is mostly FragPipe/Philosopher's "YYYY-MM-DD-decoys-contam-<db>.fasta.fas".
    Inside archives the same extensions are used; FragPipe/Philosopher's per-experiment
    `protein.fas` (sequences of the reported proteins, a result file) is reported but not counted
    as a database.
  * Source and date are read from file names only (rules SOURCE_RULES and DATE_PATTERNS below).
    Coverage of every parse is recorded. A date in a name may be a download, release, build or
    run date; "database age at deposit" = dataset submission date minus the earliest valid date in
    the name (dates after the submission date are counted and excluded).
  * Identical FASTA: same checksum, top-level files only. PRIDE checksums are present almost
    only for datasets released 2016-2023 (see Checksum-coverage-by-year.csv); identical
    downloads of the same public release by different labs are identical files, so a shared
    checksum is not evidence that one group reused another's file. Exact-size matching is
    validated against checksums inside the window and only used as a labelled sensitivity.
  * Deposition rates: dataset-level, by release year, species (human / not), software family,
    results form, DIA; group-bootstrap 95% CIs clustered on pi_group (common.cluster_boot).

H7  Storage anatomy and duplication
  * Top-level stored bytes by layer (layer_rules.RULES), archives separated from loose files.
  * Members of listed zips: uncompressed bytes by role (MEMBER_ROLES below, first match wins);
    per-archive compression ratio by dominant role. Uncompressed and stored bytes are never summed.
  * Duplication: (a) raw runs deposited both loose and zipped in the same dataset (run-name match,
    sizes compared); (b) identical checksums among top-level files within a dataset and across
    datasets (same vs different pi_group), by layer, inside the checksum window; (c) the same
    member (base name + exact uncompressed size, >= 1 MB) in two zips of one dataset, and loose
    top-level files repeated inside a zip of the same dataset (base name + exact size, >= 1 MB).
"""
import json, re, shutil
import numpy as np, pandas as pd, duckdb
from common import R, O, SEED, cluster_boot
from layer_rules import RULES
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Arsenal/H4-H7-databases-storage'; OUT.mkdir(parents=True, exist_ok=True)
AC = O / 'Archive-contents'
MEMB = str(AC / 'Members.parquet'); ARCH = str(AC / 'Archives.csv.gz')
TMP = OUT / '.duckdb_tmp'
c = connect()
c.execute(f"set threads=2; set memory_limit='3GB'; set enable_progress_bar=false; set temp_directory='{TMP}'")
S = {}  # summary
rng = np.random.default_rng(SEED)


def ci(values, groups):
    lo, hi = cluster_boot(np.asarray(values, float), np.asarray(groups))
    return round(100 * lo, 1), round(100 * hi, 1)


def rate_table(df, by, col, min_n=30):
    rows = []
    for k, g in df.groupby(by, observed=True):
        if len(g) < min_n:
            continue
        lo, hi = ci(g[col], g.grp)
        rows.append({**dict(zip([by] if isinstance(by, str) else by, k if isinstance(k, tuple) else (k,))),
                     'datasets': len(g), 'groups': g.grp.nunique(), 'with': int(g[col].sum()),
                     'pct': round(100 * g[col].mean(), 1), 'ci_lo': lo, 'ci_hi': hi})
    return pd.DataFrame(rows)


# ---- dataset frame --------------------------------------------------------------------------
F = pd.read_parquet(O / 'Deposit-features.parquet',
                    columns=['accession', 'release_year', 'submission_date', 'release_date', 'pi_group', 'human',
                             'species_n', 'dia_expanded', 'vendor_single', 'country'])
H = pd.read_csv(O / 'Habits/Dataset-form.csv.gz', usecols=['accession', 'form', 'family'])
D = F.merge(H, on='accession', how='left', validate='one_to_one')
D['grp'] = D.pi_group.fillna('acc:' + D.accession)          # unknown group -> own cluster
D['species'] = np.select([D.human.fillna(False).astype(bool) & (D.species_n.fillna(1) == 1),
                          D.human.fillna(False).astype(bool)], ['human only', 'human + other'], 'not human')
D.loc[D.human.isna(), 'species'] = 'unknown'
DS = D[['accession', 'release_year', 'submission_date', 'release_date', 'pi_group']].copy(); c.register('ds', DS)

# ---- top-level files ------------------------------------------------------------------------
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
c.execute(f"""
create temp table tf as
with f as (
  select project_accession accession, file_name, lower(file_name) n,
         lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
         coalesce(file_category_value, '') cat, coalesce(file_size_bytes, 0) sz, coalesce(checksum, '') checksum
  from pride_files
  where project_accession in (select accession from study)
    and coalesce(primary_public_location_value, '') not like '%/generated/%')
select accession, file_name, n, b, cat, sz, checksum, {case} layer,
       regexp_matches(n, '\\.(zip|tgz|tar|rar|7z)$|\\.tar\\.gz$') is_archive,
       regexp_matches(n, '\\.(zip|tgz|tar|rar|7z|gz|bz2|xz)$') is_compressed
from f""")
S['top_level_files'] = c.execute("select count(*) files, count(distinct accession) datasets, sum(sz)/1e15 PB from tf").df().iloc[0].to_dict()

# checksum coverage (affects every checksum analysis)
cov = c.execute("""select d.release_year, count(*) files, avg((checksum <> '')::int) files_with_checksum,
    sum(case when checksum <> '' then sz else 0 end) / sum(sz) bytes_with_checksum, count(distinct tf.accession) datasets,
    count(distinct case when checksum = '' then tf.accession end) datasets_missing_any
  from tf join ds d using (accession) group by 1 order by 1""").df()
cov.round(4).to_csv(OUT / 'Checksum-coverage-by-year.csv', index=False)
WINDOW = (2016, 2023)   # release years with >= ~74% of files checksummed
S['checksum_missing'] = c.execute("""select count(*) filter (where checksum = '') files_without, count(*) files,
    sum(sz) filter (where checksum = '') / 1e15 PB_without, sum(sz) / 1e15 PB from tf""").df().iloc[0].to_dict()
S['checksum_window'] = {'release_years': WINDOW, 'note': 'files without checksum are empty strings in pride_files; '
                        '2013-2015 and 2024-2025 releases are almost entirely without checksums'}

# =============================================================================================
# H4  FASTA
# =============================================================================================
FA_RE = r"\.(fasta|fa|faa|fas)$"
fa_top = c.execute(f"""select accession, file_name, sz, checksum, cat, layer from tf where regexp_matches(b, '{FA_RE}') order by all""").df()
c.execute(f"create temp view arch as select * from read_csv('{ARCH}')")
fa_mem = c.execute(f"""
select a.accession, m.url, m.member, regexp_extract(m.member, '[^/]*$') file_name, m.usize sz
from read_parquet('{MEMB}') m join arch a using (url)
where regexp_matches(lower(m.member), '\\.(fasta|fa|faa|fas)(\\.gz)?$')
  and not regexp_matches(m.member, '(^|/)__MACOSX/|(^|/)\\._') order by all""").df()
fa_mem['tool_protein_fas'] = fa_mem.file_name.str.lower().eq('protein.fas')
S['H4_fasta_files'] = {
    'top_level_files': len(fa_top), 'top_level_datasets': int(fa_top.accession.nunique()),
    'top_level_GB': round(fa_top.sz.sum() / 1e9, 1),
    'top_level_by_layer_rule': fa_top.layer.value_counts().to_dict(),
    'top_level_dot_fas_files': int(fa_top.file_name.str.lower().str.contains(r'\.fas(?:\.gz)?$').sum()),
    'inside_archives_members': len(fa_mem), 'inside_archives_datasets': int(fa_mem.accession.nunique()),
    'inside_archives_protein_fas_members_excluded': int(fa_mem.tool_protein_fas.sum()),
    'inside_archives_protein_fas_datasets': int(fa_mem[fa_mem.tool_protein_fas].accession.nunique()),
    'inside_archives_database_members': int((~fa_mem.tool_protein_fas).sum()),
    'inside_archives_database_datasets': int(fa_mem[~fa_mem.tool_protein_fas].accession.nunique())}

# ---- source rules (file names only; flags are not exclusive) ----------------------------------
SOURCE_RULES = {
    'UniProt proteome ID (UP0...)': r'up0\d{8}',
    'Swiss-Prot / reviewed (explicit)': r'swiss.?prot|sprot|(?<!un)reviewed',
    'Swiss-Prot ("sp" token only)': r'(?:^|[_.\- ])sp(?:[_.\- ]|$)',
    'TrEMBL / unreviewed': r'trembl|unreviewed',
    'UniRef': r'uniref',
    'UniProt (any mention)': r'uniprot|uniref|uniparc|swiss.?prot|sprot|trembl|up0\d{8}|refprot|ref.?proteome|reference.?proteome',
    'RefSeq': r'refseq|gcf_\d|(?:^|[^a-z])(?:nc|nz|np|xp)_\d{6,}',
    'NCBI / GenBank (not RefSeq)': r'ncbi|genbank|gca_\d',
    'Ensembl / GENCODE': r'ensembl|gencode|\.pep\.(?:all|abinitio)|grch3\d|grcm3\d',
    'IPI (discontinued 2011)': r'(?:^|[^a-z])ipi(?:[^a-z]|$)',
    'Organism database (TAIR/Araport, SGD, FlyBase, WormBase, VEuPathDB, CGD, JGI/Phytozome...)':
        r'tair\d|araport|(?:^|[^a-z])sgd|flybase|wormbase|tritrypdb|plasmodb|toxodb|cryptodb|giardiadb|veupathdb|'
        r'(?:^|[^a-z])cgd|phytozome|(?:^|[^a-z])jgi|maizegdb|xenbase|pombase|dictybase|cyanobase',
    'Genome-, transcript- or sample-derived (custom, variant, ORF, translation, metagenome, assembly...)':
        r'custom|variant|mutant|mutat|(?:^|[^a-z])snv|(?:^|[^a-z])snp|indel|(?:^|[^a-z])orfs?(?:[^a-z]|$)|translat|'
        r'frame|transcriptom|trinity|rna.?seq|de.?novo|metagenom|metaproteom|assembl|augustus|prodigal|genemark|'
        r'personali|sample.?specific|neoantigen|fusion|proteogenom|nonchromosomal|contig|gene.?set|(?:^|[^a-z])mags?(?:[^a-z]|$)',
    'Contaminants (cRAP / contaminant list)': r'contam|crap|(?:^|[^a-z])cont(?:[^a-z]|$)',
    'Decoys included': r'(?<!no)(?<!no_)(?<!no-)decoy|revers|(?:^|[^a-z])rev(?:[^a-z]|$)|shuffl',
    'iRT / spike-in standards': r'(?:^|[^a-z])irt(?:[^a-z]|$)|biognosys|(?:^|[^a-z])ups[12]?(?:[^a-z]|$)|spike',
    'FragPipe/Philosopher decoy-annotated (.fas)': r'decoys?-.*contam|\.fas$',
    'MaxQuant contaminants.fasta': r'^contaminants\.fasta$',
}
PUBLIC_ORDER = ['UniProt (any mention)', 'RefSeq', 'NCBI / GenBank (not RefSeq)', 'Ensembl / GENCODE', 'IPI (discontinued 2011)',
                'Organism database (TAIR/Araport, SGD, FlyBase, WormBase, VEuPathDB, CGD, JGI/Phytozome...)',
                'Genome-, transcript- or sample-derived (custom, variant, ORF, translation, metagenome, assembly...)']
UP_SPECIES = {'up000005640': 'human', 'up000000589': 'mouse', 'up000002494': 'rat', 'up000002311': 'yeast S. cerevisiae',
              'up000000625': 'E. coli K-12', 'up000006548': 'Arabidopsis', 'up000000803': 'Drosophila',
              'up000001940': 'C. elegans', 'up000000437': 'zebrafish', 'up000009136': 'cattle', 'up000008227': 'pig',
              'up000000539': 'chicken', 'up000001584': 'M. tuberculosis H37Rv'}
MON = dict(jan=1, feb=2, mar=3, apr=4, may=5, jun=6, jul=7, aug=8, sep=9, oct=10, nov=11, dec=12)
MONTH = r'(jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|june?|july?|aug(?:ust)?|sep(?:t|tember)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)'
Y4 = r'((?:19|20)\d\d)'
# (format, regex, (year, month, day) group indices or callable, confidence). Earlier patterns claim spans first.
DATE_PATTERNS = [
    ('YYYY-MM-DD', Y4 + r'[-_.](0[1-9]|1[0-2])[-_.](0[1-9]|[12]\d|3[01])', (1, 2, 3), 'high'),
    ('YYYYMMDD', Y4 + r'(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])', (1, 2, 3), 'high'),
    ('DD-MM-YYYY', r'(0[1-9]|[12]\d|3[01])[-_.](0[1-9]|1[0-2])[-_.]' + Y4, (3, 2, 1), 'high'),
    ('DDMMYYYY', r'(0[1-9]|[12]\d|3[01])(0[1-9]|1[0-2])' + Y4, (3, 2, 1), 'high'),
    ('month name + YYYY', r'(?<![a-z])' + MONTH + r'(?![a-z])[-_. ]?(?:(?:0?[1-9]|[12]\d|3[01])[-_. ,]+)?' + Y4, (2, 1, None), 'high'),
    ('YYYY + month name', Y4 + r'[-_. ]?' + MONTH + r'(?![a-z])', (1, 2, None), 'high'),
    ('YYYY-MM (incl. UniProt release YYYY_MM)', Y4 + r'[-_.](0[1-9]|1[0-2])', (1, 2, None), 'high'),
    ('month name + YY', r'(?<![a-z])' + MONTH + r'(?![a-z])[-_. ]?(?:(?:0?[1-9]|[12]\d|3[01])[-_. ,]+)?([012]\d)', (2, 1, None), 'medium'),
    ('YYYYMM', r'(20[0-2]\d)(0[1-9]|1[0-2])', (1, 2, None), 'medium'),
    ('YYMMDD', r'(0[5-9]|1\d|2[0-5])(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])', (1, 2, 3), 'medium'),
    ('MM-YYYY', r'(0[1-9]|1[0-2])[-_.]' + Y4, (2, 1, None), 'medium'),
    ('MMYYYY', r'(0[1-9]|1[0-2])' + Y4, (2, 1, None), 'medium'),
    ('YYYY_M (UniProt release, one-digit month)', Y4 + r'[_.]([1-9])', (1, 2, None), 'medium'),
    ('YYYY only', r'(20[0-2]\d)', (1, None, None), 'low'),
]
DATE_RX = [(f, re.compile(('' if f.startswith('month') else r'(?<!\d)') + p + r'(?!\d)'), g, conf) for f, p, g, conf in DATE_PATTERNS]
CONF_RANK = {'high': 0, 'medium': 1, 'low': 2}
YYMMDD_RX = re.compile(r'^(0[5-9]|1\d|2[0-5])(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])$')
SNAP = pd.Timestamp('2026-04-07')


def parse_dates(name):
    """All non-overlapping dates in a lower-cased file name -> list of (date, format, confidence)."""
    s = re.sub(r'up0\d{8}|gc[af]_\d+(\.\d+)?|(?:nc|nz|np|xp|wp|nm|xm)_\d+(\.\d+)?|(?<![a-z0-9])\d{4,6}(?=_?(entries|seqs|sequences|proteins|prot\b))', ' ', name)
    taken, out = [], []
    for f, rx, (gy, gm, gd), conf in DATE_RX:
        for mt in rx.finditer(s):
            if any(mt.start() < e and mt.end() > b for b, e in taken):
                continue
            y = int(mt.group(gy)); y = y + 2000 if y < 100 else y
            m = 7 if gm is None else (int(mt.group(gm)) if mt.group(gm).isdigit() else MON[mt.group(gm)[:3]])
            d = (1 if gm is None else 15) if gd is None else int(mt.group(gd))
            try:
                dt = pd.Timestamp(year=y, month=m, day=d)
            except ValueError:
                continue
            if not (pd.Timestamp('1995-01-01') <= dt <= SNAP):
                continue
            taken.append((mt.start(), mt.end()))
            if f == 'YYYYMM' and YYMMDD_RX.match(mt.group(0)):
                continue          # six digits valid both as YYYYMM and YYMMDD (e.g. 200510): ambiguous, not used
            out.append((dt, f, conf))
    return out


def describe_fasta(df):
    n = df.file_name.str.lower()
    for k, p in SOURCE_RULES.items():
        df[k] = n.str.contains(p, regex=True)
    df['public_source'] = np.select([df[k] for k in PUBLIC_ORDER], PUBLIC_ORDER, 'no source in name')
    df['uniprot_species'] = n.str.extract(r'(up0\d{8})', expand=False).map(UP_SPECIES)
    df['version_token'] = n.str.contains(r'(?:^|[^a-z])(?:v|ver|version|release|rel)[-_. ]?\d', regex=True)
    df['size_band'] = pd.cut(df.sz, [-1, 1e4, 1e6, 1e8, 1e13], labels=['<10 kB', '10 kB-1 MB', '1-100 MB', '>100 MB'])
    parsed = n.map(parse_dates)
    df['n_dates'] = parsed.map(len)
    best = parsed.map(lambda L: min(L, key=lambda t: (CONF_RANK[t[2]], t[0])) if L else (pd.NaT, None, None))
    df['name_date'] = best.map(lambda t: t[0]); df['date_format'] = best.map(lambda t: t[1]); df['date_conf'] = best.map(lambda t: t[2])
    return df


fa_top = describe_fasta(fa_top)
fa_mem = describe_fasta(fa_mem)
fa_top = fa_top.merge(D[['accession', 'release_year', 'submission_date', 'pi_group', 'grp', 'species', 'family']], on='accession')
fa_mem = fa_mem.merge(D[['accession', 'release_year', 'submission_date', 'pi_group', 'grp', 'species', 'family']], on='accession')
fa_db_mem = fa_mem[~fa_mem.tool_protein_fas]

# spot-check sample of the parse (reviewed by hand; written for re-checking)
spot = pd.concat([fa_top.sample(80, random_state=SEED).assign(where='top level'),
                  fa_db_mem.drop_duplicates('file_name').sample(40, random_state=SEED).assign(where='inside archive')])
spot[['where', 'accession', 'file_name', 'sz', 'public_source', 'Contaminants (cRAP / contaminant list)', 'Decoys included',
      'name_date', 'date_format', 'date_conf', 'submission_date']].to_csv(OUT / 'H4-spot-check-sample.csv', index=False)

# coverage of each parse (files and datasets; top level and inside archives)
cov_rows = []
for where, df in [('top level', fa_top), ('inside archives (excl. protein.fas)', fa_db_mem)]:
    for k in list(SOURCE_RULES) + ['date (any)', 'date high-confidence', 'version token', 'no source in name', 'size <10 kB']:
        s = (df[k] if k in SOURCE_RULES else df.name_date.notna() if k == 'date (any)'
             else df.date_conf.eq('high') if k == 'date high-confidence' else df.version_token if k == 'version token'
             else df.public_source.eq('no source in name') if k == 'no source in name' else df.sz < 1e4)
        cov_rows.append({'where': where, 'parse': k, 'files': int(s.sum()), 'pct_files': round(100 * s.mean(), 1),
                         'datasets': int(df[s].accession.nunique()),
                         'pct_datasets_with_fasta': round(100 * df[s].accession.nunique() / df.accession.nunique(), 1)})
cov_rows = pd.DataFrame(cov_rows); cov_rows.to_csv(OUT / 'H4-name-parse-coverage.csv', index=False)
fa_top.public_source.value_counts().rename('files').to_frame().assign(
    datasets=fa_top.groupby('public_source').accession.nunique()).to_csv(OUT / 'H4-public-source-top-level.csv')
fa_top.date_format.value_counts(dropna=False).to_csv(OUT / 'H4-date-formats-top-level.csv')
fa_top.uniprot_species.value_counts().to_csv(OUT / 'H4-uniprot-proteome-ids-top-level.csv')

# source by release period (dataset level, top level + archives)
allfa = pd.concat([fa_top.assign(where='top'), fa_db_mem.assign(where='archive')], ignore_index=True)
dsrc = allfa.groupby('accession')[list(SOURCE_RULES)].any()
dsrc['any public source named'] = allfa.assign(p=allfa.public_source.isin(PUBLIC_ORDER[:6])).groupby('accession').p.any()
dsrc['only sources not named'] = ~allfa.assign(p=allfa.public_source.ne('no source in name')).groupby('accession').p.any()
dsrc = dsrc.join(D.set_index('accession')[['release_year', 'grp']])
dsrc['period'] = pd.cut(dsrc.release_year, [2012, 2016, 2019, 2022, 2025], labels=['2013-16', '2017-19', '2020-22', '2023-25'])
src_period = dsrc.groupby('period', observed=True)[list(SOURCE_RULES) + ['any public source named', 'only sources not named']].mean().mul(100).round(1).T
src_period.columns = src_period.columns.astype(str); src_period.index.name = 'name flag (pct of datasets with FASTA)'
src_period.loc['datasets with FASTA (denominator, n)'] = dsrc.groupby('period', observed=True).size().to_numpy()
src_period.to_csv(OUT / 'H4-sources-by-period-dataset-pct.csv')

# ---- database age at deposit ------------------------------------------------------------------
def age_table(df, label):
    df = df[df.name_date.notna()].copy()
    df['age_days'] = (df.submission_date - df.name_date).dt.days
    df['future'] = df.age_days < -31          # date more than a month after submission
    ok = df[~df.future]
    # one value per dataset: the largest dated FASTA (likely the main database)
    per_ds = ok.sort_values('sz', ascending=False).drop_duplicates('accession')
    return df, per_ds


dated_all, age_ds = age_table(pd.concat([fa_top.assign(where='top'), fa_db_mem.assign(where='archive')], ignore_index=True), 'all')
age_ds['sub_year'] = age_ds.submission_date.dt.year
age_ds['years'] = age_ds.age_days / 365.25
age_ds['gt1y'] = age_ds.years > 1; age_ds['gt2y'] = age_ds.years > 2; age_ds['gt5y'] = age_ds.years > 5
age_ds['uniprot'] = age_ds['UniProt (any mention)']


def age_summary(df, by):
    rows = []
    for k, g in df.groupby(by, observed=True):
        if len(g) < 20:
            continue
        lo, hi = ci(g.gt2y, g.grp)
        rows.append({by: k, 'datasets': len(g), 'groups': g.grp.nunique(),
                     'median_years': round(g.years.median(), 2), 'q25': round(g.years.quantile(.25), 2), 'q75': round(g.years.quantile(.75), 2),
                     'pct_gt1y': round(100 * g.gt1y.mean(), 1), 'pct_gt2y': round(100 * g.gt2y.mean(), 1), 'gt2y_ci_lo': lo, 'gt2y_ci_hi': hi,
                     'pct_gt5y': round(100 * g.gt5y.mean(), 1)})
    return pd.DataFrame(rows)


hi = age_ds[age_ds.date_conf == 'high']
# per-dataset database name date (largest dated FASTA per dataset), for cross-clock comparisons (arsenal_clocks.py)
age_ds[['accession', 'name_date', 'date_conf', 'date_format', 'where', 'file_name', 'years']].to_csv(OUT / 'H4-dataset-database-date.csv.gz', index=False)
age_by_year = age_summary(hi, 'sub_year'); age_by_year.to_csv(OUT / 'H4-database-age-by-submission-year-high-confidence.csv', index=False)
age_summary(age_ds[age_ds.date_conf.isin(['high', 'medium'])], 'sub_year').to_csv(OUT / 'H4-database-age-by-submission-year-high+medium.csv', index=False)
age_uni = age_summary(hi[hi.uniprot], 'sub_year'); age_uni.to_csv(OUT / 'H4-database-age-by-submission-year-UniProt-named.csv', index=False)
age_fmt = age_summary(age_ds, 'date_format'); age_fmt.to_csv(OUT / 'H4-database-age-by-date-format.csv', index=False)
age_src = age_summary(hi.assign(src=np.where(hi.uniprot, 'UniProt named', hi.public_source)), 'src')
age_fam = age_summary(hi[hi.sub_year >= 2019], 'family'); age_fam.to_csv(OUT / 'H4-database-age-2019-2025-by-software.csv', index=False)
age_sp = age_summary(hi[hi.sub_year >= 2019], 'species'); age_sp.to_csv(OUT / 'H4-database-age-2019-2025-by-species.csv', index=False)
# FragPipe/Philosopher .fas: the leading date is the build date; the database inside may be older
fas = dated_all[dated_all.file_name.str.lower().str.contains(r'decoys?-.*contam', regex=True) & ~dated_all.future]
S['H4_age'] = {
    'fasta_files_with_any_date': int(dated_all.shape[0]),
    'fasta_files_total': int(len(fa_top) + len(fa_db_mem)),
    'files_date_more_than_31d_after_submission_excluded': int(dated_all.future.sum()),
    'datasets_with_dated_fasta_any_conf': int(len(age_ds)), 'datasets_high_conf': int(len(hi)),
    'datasets_with_any_fasta': int(dsrc.shape[0]),
    'high_conf_median_years': round(hi.years.median(), 2), 'high_conf_pct_gt1y': round(100 * hi.gt1y.mean(), 1),
    'high_conf_pct_gt2y': round(100 * hi.gt2y.mean(), 1), 'high_conf_pct_gt5y': round(100 * hi.gt5y.mean(), 1),
    'uniprot_named_high_conf_median_years': round(hi[hi.uniprot].years.median(), 2),
    'philosopher_fas_build_date_files': int(len(fas)),
    'philosopher_fas_median_days_build_to_submission': float((fas.age_days).median()) if len(fas) else None,
    'by_source_high_conf': age_src.to_dict('records')}

# ---- identical FASTA (checksum, top level) ----------------------------------------------------
ck = fa_top[(fa_top.checksum != '') & (fa_top.sz > 0)].copy()
g = ck.groupby('checksum')
ident = pd.DataFrame({'files': g.size(), 'datasets': g.accession.nunique(), 'groups': g.grp.nunique(),
                      'size_MB': g.sz.max() / 1e6, 'first_year': g.release_year.min(), 'last_year': g.release_year.max(),
                      'names': g.file_name.agg(lambda s: ' | '.join(pd.unique(s)[:4])),
                      'public_source': g.public_source.agg(lambda s: s.mode().iat[0]),
                      'contaminants': g['Contaminants (cRAP / contaminant list)'].any(),
                      'example_accessions': g.accession.agg(lambda s: ' '.join(sorted(pd.unique(s))[:5]))})
ident.sort_values(['groups', 'datasets'], ascending=False).head(60).to_csv(OUT / 'H4-identical-fasta-most-redeposited.csv')
ck['ds_per_ck'] = ck.checksum.map(ident.datasets); ck['grp_per_ck'] = ck.checksum.map(ident.groups)
dsck = ck.groupby('accession').agg(any_shared_ds=('ds_per_ck', lambda s: (s > 1).any()), any_shared_grp=('grp_per_ck', lambda s: (s > 1).any()),
                                   grp=('grp', 'first'), release_year=('release_year', 'first'))
# one group depositing the same FASTA file (checksum) in several datasets: time between first and last submission
sg = ck[ck.pi_group.notna()].groupby(['checksum', 'pi_group']).agg(datasets=('accession', 'nunique'), first=('submission_date', 'min'),
                                                                   last=('submission_date', 'max'), size_MB=('sz', lambda x: x.max() / 1e6))
sg = sg[sg.datasets > 1].copy(); sg['span_years'] = (sg['last'] - sg['first']).dt.days / 365.25
sg.sort_values('span_years', ascending=False).head(40).to_csv(OUT / 'H4-same-group-same-fasta-longest-spans.csv')
multi = ident[ident.groups > 1]
# internal check of name dates: identical files (same checksum) whose names carry different high-confidence dates
dd = ck[ck.date_conf.eq('high')].groupby('checksum').name_date.agg(['min', 'max', 'count'])
dd = dd[dd['count'] > 1]; spread = (dd['max'] - dd['min']).dt.days
# exact-size matching as a proxy (validated inside the window): pairs of datasets with same-size FASTA (>100 kB)
big = ck[ck.sz > 1e5][['accession', 'sz', 'checksum']].drop_duplicates()
pairs = big.merge(big, on='sz'); pairs = pairs[pairs.accession_x < pairs.accession_y]
size_precision = float((pairs.checksum_x == pairs.checksum_y).mean()) if len(pairs) else None
fa_size = fa_top[fa_top.sz > 1e5][['accession', 'sz', 'grp', 'release_year']].drop_duplicates(['accession', 'sz'])
gs = fa_size.groupby('sz').agg(datasets=('accession', 'nunique'), groups=('grp', 'nunique'))
fa_size['shared_grp'] = fa_size.sz.map(gs.groups) > 1
size_by_year = fa_size.groupby(fa_size.accession).agg(shared=('shared_grp', 'any'), y=('release_year', 'first')).groupby('y').shared.mean().mul(100).round(1)
S['H4_identical'] = {
    'top_level_fasta_with_checksum': int(len(ck)), 'of_top_level_fasta': int(len(fa_top)),
    'distinct_checksums': int(len(ident)),
    'checksums_in_2plus_datasets': int((ident.datasets > 1).sum()), 'checksums_in_2plus_groups': int((ident.groups > 1).sum()),
    'files_whose_checksum_in_2plus_groups_pct': round(100 * (ck.grp_per_ck > 1).mean(), 1),
    'datasets_with_checksummed_fasta': int(len(dsck)),
    'datasets_sharing_a_fasta_with_another_dataset_pct': round(100 * dsck.any_shared_ds.mean(), 1),
    'datasets_sharing_a_fasta_with_another_group_pct': round(100 * dsck.any_shared_grp.mean(), 1),
    'datasets_sharing_with_another_group_ci': ci(dsck.any_shared_grp, dsck.grp),
    'same_group_same_fasta_pairs_(checksum x group, >=2 datasets)': int(len(sg)),
    'same_group_same_fasta_median_span_years': round(float(sg.span_years.median()), 2),
    'same_group_same_fasta_span_gt1y_pct': round(100 * float((sg.span_years > 1).mean()), 1),
    'same_group_same_fasta_span_gt3y_pct': round(100 * float((sg.span_years > 3).mean()), 1),
    'checksums_in_2plus_groups_by_name_source': multi.public_source.value_counts().to_dict(),
    'checksums_in_2plus_groups_contaminant_named': int(multi.contaminants.sum()),
    'identical_files_with_2plus_dated_names': int(len(dd)),
    'identical_files_name_date_spread_median_days': float(spread.median()) if len(dd) else None,
    'identical_files_name_date_spread_le_92_days_pct': round(100 * float((spread <= 92).mean()), 1) if len(dd) else None,
    'size_only_match_precision_vs_checksum_(>100kB, pairs)': [round(size_precision, 3) if size_precision is not None else None, int(len(pairs))],
    'size_match_shared_with_other_group_pct_by_year_(sensitivity, all years)': size_by_year.to_dict()}

# ---- FASTA deposition rates -------------------------------------------------------------------
D['fasta_top'] = D.accession.isin(fa_top.accession)
D['fasta_top_layer_rule'] = D.accession.isin(fa_top[fa_top.layer == 'sequence_db'].accession)
D['fasta_anywhere'] = D.fasta_top | D.accession.isin(fa_db_mem.accession)
D['period'] = pd.cut(D.release_year, [2012, 2016, 2019, 2022, 2025], labels=['2013-16', '2017-19', '2020-22', '2023-25'])
D['mode'] = np.where(D.dia_expanded, 'DIA', 'not DIA')
rates = []
for col in ['fasta_top', 'fasta_anywhere', 'fasta_top_layer_rule']:
    t = rate_table(D, 'release_year', col); t.insert(0, 'measure', col); rates.append(t)
pd.concat(rates).to_csv(OUT / 'H4-fasta-rate-by-year.csv', index=False)
strat = []
for by in ['species', 'family', 'form', 'mode', 'vendor_single']:
    for per, sub in [('2013-2018', D[D.release_year <= 2018]), ('2019-2025', D[D.release_year >= 2019])]:
        t = rate_table(sub, by, 'fasta_anywhere'); t.insert(0, 'period', per); t.insert(0, 'dimension', by)
        t = t.rename(columns={by: 'level'}); strat.append(t)
strat = pd.concat(strat); strat.to_csv(OUT / 'H4-fasta-rate-by-species-software-form-mode.csv', index=False)
yh = rate_table(D[D.species.isin(['human only', 'not human'])], ['release_year', 'species'], 'fasta_anywhere')
yh.to_csv(OUT / 'H4-fasta-rate-by-year-and-species.csv', index=False)
fs = rate_table(D[(D.release_year >= 2019) & D.species.isin(['human only', 'not human'])], ['family', 'species'], 'fasta_anywhere')
fs.to_csv(OUT / 'H4-fasta-rate-2019-2025-software-x-species.csv', index=False)
# year-standardised human vs not-human difference (weights = all datasets per year), group bootstrap
def std_diff(df):
    w = df.release_year.value_counts(normalize=True)
    r = df.groupby(['release_year', 'species']).fasta_anywhere.mean().unstack()
    return float(((r['not human'] - r['human only']) * w).sum())
HS = D[D.species.isin(['human only', 'not human'])]
est = std_diff(HS)
grps = HS.grp.unique(); idx = HS.groupby('grp').indices; boots = []
for _ in range(500):
    pick = rng.choice(grps, len(grps)); boots.append(std_diff(HS.iloc[np.concatenate([idx[x] for x in pick])]))
S['H4_rates'] = {
    'fasta_top_level_pct': round(100 * D.fasta_top.mean(), 1), 'fasta_anywhere_pct': round(100 * D.fasta_anywhere.mean(), 1),
    'datasets': len(D),
    'fasta_anywhere_by_period': D.groupby('period', observed=True).fasta_anywhere.mean().mul(100).round(1).astype(str).to_dict(),
    'not_human_minus_human_only_year_standardised_pp': [round(100 * est, 1), round(100 * np.quantile(boots, .025), 1), round(100 * np.quantile(boots, .975), 1)],
    'note': 'archives that could not be listed (tar.gz, rar, 7z...) are not opened; fasta_anywhere is a lower bound'}

# =============================================================================================
# H7  storage anatomy
# =============================================================================================
top = c.execute("""select d.release_year, layer, is_archive, count(*) files, sum(sz) bytes
                   from tf join ds d using (accession) group by all""").df()
lay = top.groupby(['layer', 'is_archive']).agg(files=('files', 'sum'), bytes=('bytes', 'sum')).reset_index()
lay['TB'] = (lay.bytes / 1e12).round(2); lay['pct_bytes'] = (100 * lay.bytes / lay.bytes.sum()).round(2)
lay.sort_values('bytes', ascending=False).drop(columns='bytes').to_csv(OUT / 'H7-top-level-bytes-by-layer.csv', index=False)
ly = top.assign(l=np.where(top.is_archive & ~top.layer.eq('vendor_raw'), 'archives (non-raw)', top.layer)).pivot_table(
    index='release_year', columns='l', values='bytes', aggfunc='sum', fill_value=0)
(ly.div(ly.sum(axis=1), axis=0) * 100).round(2).assign(total_TB=(ly.sum(axis=1) / 1e12).round(1)).to_csv(OUT / 'H7-top-level-byte-share-by-year.csv')

# ---- member roles ------------------------------------------------------------------------------
MQ_FINAL = r'(proteingroups|peptides|evidence|modificationspecificpeptides|summary|parameters|[a-z0-9 ()_-]*sites)\.txt$'
MQ_SCAN = r'(allpeptides|msms|msmsscans|msscans|ms3scans|pasefmsmsscans|accumulated[a-z]*scans|mzrange|matchedfeatures|librarymatch[a-z]*|msmsspectruminfo)( ?\d+)?\.txt$'
MEMBER_ROLES = [  # (role, SQL condition on n = lower path, bn = lower base name, ext) — first match wins
    ('standard result files (mzIdentML/mzTab)', r"regexp_matches(bn, '\.(mzid|mzidentml|mztab)(\.txt)?(\.gz)?$')"),
    ('open spectra (mzML, MGF...)', r"regexp_matches(bn, '\.(mzml|mzxml|mgf|ms1|ms2|mzdata|mz5|imzml|ibd|dta|pkl)(\.gz)?$')"),
    ('raw data', r"""regexp_matches(n, '\.d/') or regexp_matches(n, '\.raw/_')
        or regexp_matches(bn, '\.(raw|wiff|wiff2|scan|tdf|tdf_bin|tsf|tsf_bin|baf|baf_idx|baf_xtr|yep|lcd|t2d|uimf)$')"""),
    ('sequence databases (FASTA)', r"regexp_matches(bn, '\.(fasta|fa|faa|fas)(\.gz)?$')"),
    ('spectral libraries', r"regexp_matches(bn, '\.(speclib|blib|msp|sptxt|dlib|elib|kit)$') or bn like '%speclib%' or bn like '%spectral_lib%'"),
    ('final result tables', rf"""regexp_matches(bn, '^{MQ_FINAL}')
        or regexp_matches(bn, '(^|[_.-])report(\.pr_matrix|\.pg_matrix|\.gg_matrix|\.unique_genes_matrix|\.stats)?\.(tsv|parquet)$')
        or regexp_matches(bn, '(pr|pg|gg|unique_genes)_matrix\.(tsv|parquet)$')
        or regexp_matches(bn, '^(psm|peptide|protein|ion|combined_protein|combined_peptide|combined_ion|combined_modified_peptide|combined_site_[a-z0-9_.-]*|abundance_[a-z0-9_]*|ratio_[a-z0-9_]*)\.tsv$')
        or regexp_matches(bn, '_(proteins|psms|peptidegroups|proteingroups|peptides|modificationsites|consensusfeatures)\.(txt|xlsx|tsv|csv)$')
        or regexp_matches(bn, '^(protein-peptides|proteins|peptides|db search psm|de novo peptides|peptide|protein)\.csv$')"""),
    ('scan-level and intermediate search files', rf"""regexp_matches(bn, '{MQ_SCAN}')
        or ext in ('apl', 'aplbin', 'res', 'index', 'ind', 'sil', 'iso', 'iden_sil', 'idenall_sil', 'msmsmat', 'quant', 'pin',
                   'pepxml', 'mzbin', 'mzbin_calibrated', 'pepindex', 'pfd', 'data', 'byspec2', 'ser', 'idenall_peak', 'apar', 'spo')
        or (ext = 'out' and regexp_matches(bn, '\.\d+\.\d+\.\d\.out$'))
        or regexp_matches(n, '(^|/)[^/]+/[pn]\d+/[^/]+$')
        or ext like 'peaks%' or bn like '%.pep.xml'
        or regexp_matches(n, '(^|/)combined[^/]*/(search|andromeda|ser)/') or regexp_matches(n, '/tmps?/')
        or regexp_matches(n, '(^|/)search/a-?\d+/') or regexp_matches(n, '\.mse/') or bn like '%_ia_workflow%'
        or (ext = '' and regexp_matches(n, '(^|/)combined[^/]*/'))"""),
    ('tool project and result databases', r"""ext in ('msf', 'pdresult', 'pdstudy', 'sne', 'sf3', 'sky', 'skyd', 'group', 'byrslt', 'idpdb',
                   'psdb', 'cpsx', 'protxml', 'idxml', 'featurexml', 'consensusxml', 'p')
        or bn like '%.prot.xml' or regexp_matches(n, '/runinfo/|(^|/)seg\d+/') or regexp_matches(bn, '(^|[^a-z0-9])f\d{4,}\.dat$') or regexp_matches(bn, 'mascot.*\.dat$')"""),
    ('other tables', r"regexp_matches(bn, '\.(tsv|csv|xlsx|xls|txt|parquet|tab)$')"),
    ('parameters, documents, images and logs', r"""regexp_matches(bn, '\.(xml|json|ini|yml|yaml|pdf|docx?|pptx?|md|html?|png|jpe?g|tiff?|svg|emf|wmf|gif|bmp|log|params?|cfg|conf|method|sqlite|db)$')"""),
    ('nested archives', r"regexp_matches(bn, '\.(zip|tgz|tar|rar|7z|gz|bz2|xz)$')"),
]
role_case = 'case ' + ' '.join(f"when {cond} then '{r}'" for r, cond in MEMBER_ROLES) + " else 'other' end"
DETAIL = rf"""case when regexp_matches(bn, '{MQ_SCAN}') or regexp_matches(bn, '^{MQ_FINAL}')
                   then 'MaxQuant ' || regexp_replace(regexp_extract(bn, '([a-z0-9 ()]+)( ?\d+)?\.txt$', 1), '^.*sites$', '*sites')
                   when regexp_matches(n, '\.d/') then 'Bruker/Agilent .d folder'
                   when ext like 'peaks%' then '.peaks* (MaxQuant)'
                   when regexp_matches(n, '/runinfo/|(^|/)seg\d+/') or ext = 'p' then 'PEAKS project folder'
                   when ext = '' then '(no extension)' else '.' || ext end"""
c.execute(f"""create temp view mv as
  select a.accession, m.url, m.member, m.usize, n, bn, ext, {role_case} mrole
  from (select url, member, usize, lower(member) n, lower(regexp_extract(member, '[^/]*$')) bn,
               regexp_extract(lower(regexp_extract(member, '[^/]*$')), '\\.([a-z0-9_]+)$', 1) ext
        from read_parquet('{MEMB}')
        where not ends_with(member, '/') and not regexp_matches(member, '(^|/)__MACOSX/|(^|/)\\._')) m
  join arch a using (url)""")
# one materialised subset for the duplication queries: members >= 1 MB and every raw-data member
c.execute("create temp table mbig as select accession, url, member, n, bn, usize, mrole from mv where usize >= 1e6 or mrole = 'raw data'")
roles_arch = c.execute(f"""select url, accession, mrole as "role", {DETAIL} detail, count(*) members, sum(usize) usize
                           from mv group by all""").df()
roles_arch = roles_arch.merge(D[['accession', 'release_year']], on='accession', how='left').sort_values(['url', 'role', 'detail'], ignore_index=True)
rt = roles_arch.groupby('role').agg(members=('members', 'sum'), usize=('usize', 'sum'), archives=('url', 'nunique'), datasets=('accession', 'nunique'))
rt['TB_uncompressed'] = (rt.usize / 1e12).round(2); rt['pct'] = (100 * rt.usize / rt.usize.sum()).round(1)
rt.sort_values('usize', ascending=False).drop(columns='usize').to_csv(OUT / 'H7-inside-zips-bytes-by-role.csv')
det = roles_arch.groupby(['role', 'detail']).agg(members=('members', 'sum'), usize=('usize', 'sum'), datasets=('accession', 'nunique')).reset_index()
det['TB_uncompressed'] = (det.usize / 1e12).round(3); det['pct_of_all'] = (100 * det.usize / det.usize.sum()).round(2)
det = det.sort_values('usize', ascending=False); det.groupby('role').head(12).drop(columns='usize').to_csv(OUT / 'H7-inside-zips-bytes-by-role-detail.csv', index=False)
ry = roles_arch.pivot_table(index='release_year', columns='role', values='usize', aggfunc='sum', fill_value=0)
(ry.div(ry.sum(axis=1), axis=0) * 100).round(1).assign(total_TB=(ry.sum(axis=1) / 1e12).round(1),
    archives=roles_arch.groupby('release_year').url.nunique()).to_csv(OUT / 'H7-inside-zips-role-share-by-year.csv')
# 'other' spot check: what is left unclassified
det[det.role == 'other'].head(25).drop(columns='usize').to_csv(OUT / 'H7-inside-zips-unclassified-top.csv', index=False)
# role spot-check sample (random members per role)
c.execute(f"""copy (select mrole as "role", accession, member, usize from (
    select *, row_number() over (partition by mrole order by hash(member || url), member, url) r from mv where hash(url || member) % 1000 = 7)
    where r <= 15 order by mrole, member, url) to '{OUT / 'H7-role-spot-check-sample.csv'}' (header)""")

# per-archive compression by dominant role (stored size / uncompressed size)
A = pd.read_csv(ARCH); A = A[A.status.isin(['ok', 'partial'])]
dom = roles_arch.sort_values('usize', ascending=False).drop_duplicates('url')[['url', 'role']].rename(columns={'role': 'dominant_role'})
tot = roles_arch.groupby('url').usize.sum().rename('usize_total')
A = A.merge(dom, on='url', how='left').merge(tot, left_on='url', right_index=True, how='left')
A['ratio'] = A['size'] / A.usize_total.replace(0, np.nan)
comp = A[A.ratio.notna()].groupby('dominant_role').agg(archives=('url', 'size'), stored_TB=('size', lambda s: s.sum() / 1e12),
                                                       uncompressed_TB=('usize_total', lambda s: s.sum() / 1e12), median_ratio=('ratio', 'median'))
comp['pooled_ratio'] = comp.stored_TB / comp.uncompressed_TB
comp.round(3).sort_values('uncompressed_TB', ascending=False).to_csv(OUT / 'H7-zip-compression-by-dominant-role.csv')
S['H7_anatomy'] = {
    'top_level_TB_by_layer': lay.groupby('layer').TB.sum().sort_values(ascending=False).to_dict(),
    'top_level_archives_TB_by_layer': lay[lay.is_archive].set_index('layer').TB.to_dict(),
    'listed_zips': int(A.shape[0]), 'listed_zips_stored_TB': round(A['size'].sum() / 1e12, 1),
    'listed_zips_uncompressed_TB': round(A.usize_total.sum() / 1e12, 1),
    'inside_zips_pct_by_role': rt.pct.sort_values(ascending=False).to_dict(),
    'note': 'member bytes are uncompressed; compression ratio per dominant role in H7-zip-compression-by-dominant-role.csv'}

# =============================================================================================
# H7  duplication
# =============================================================================================
# (a) raw runs loose and zipped in the same dataset
RAWEXT = r'\.(wiff\.scan|raw|wiff|wiff2|d|baf|tdf|lcd|yep|uimf|t2d)$'
# run = name without compression/archive suffixes and without the raw extension; the raw extension is kept
# separately so that complementary files (SCIEX .wiff + .wiff.scan) are never matched to each other.
c.execute(f"""create temp table rawtop as
  with x as (select accession, file_name, sz, is_compressed is_archive,
                    regexp_replace(regexp_replace(n, '(\\.(zip|tgz|tar|rar|7z|gz|bz2|xz))+$', ''), '\\.tar$', '') k
             from tf where layer = 'vendor_raw')
  select *, regexp_replace(k, '{RAWEXT}', '') run, nullif(regexp_extract(k, '{RAWEXT}', 1), '') rawext from x""")
# zipped copy matches a loose file with the same run name and the same raw extension (or none named in the zip name)
a1 = c.execute("""
  with l as (select accession, run, rawext, sz loose_sz, file_name loose_name from rawtop where not is_archive),
       z as (select accession, run, rawext zext, sz zip_sz, file_name zip_name from rawtop where is_archive)
  select l.*, z.zip_sz, z.zip_name from l join z on l.accession = z.accession and l.run = z.run
   and (z.zext = l.rawext or (z.zext is null and l.rawext is not null and l.rawext <> 'wiff.scan'))""").df()
a1 = a1.drop_duplicates(['accession', 'run', 'zip_name']).sort_values(['accession', 'run', 'zip_name'], ignore_index=True)
a1 = a1.merge(D[['accession', 'release_year', 'grp']], on='accession')
a1['ratio'] = a1.zip_sz / a1.loose_sz
raw_mem = c.execute(f"""
  select accession, url,
         case when regexp_matches(n, '\\.d/') then regexp_extract(n, '([^/]+)\\.d/', 1)
              when regexp_matches(n, '\\.raw/_') then regexp_extract(n, '([^/]+)\\.raw/_', 1)
              else regexp_replace(bn, '{RAWEXT}', '') end run,
         sum(usize) usize, count(*) members
  from mbig where mrole = 'raw data' group by all order by all""").df()
c.register('raw_mem', raw_mem)
a2 = c.execute("""select r.accession, r.run, r.usize, r.url, t.loose_sz, t.archived_sz
   from raw_mem r join (select accession, run, sum(case when not is_archive and rawext <> 'wiff.scan' then sz end) loose_sz,
                               sum(case when is_archive then sz end) archived_sz from rawtop group by 1, 2) t using (accession, run)""").df()
a2 = a2.merge(D[['accession', 'release_year', 'grp']], on='accession').sort_values(['accession', 'run', 'url'], ignore_index=True)
a2['size_match_2pct'] = (a2.usize / a2.loose_sz).between(.98, 1.02)
raw_ds = c.execute("select accession, count(*) files, sum(sz) bytes from tf where layer='vendor_raw' group by 1").df()
D['raw_loose_and_zipped_top'] = D.accession.isin(a1.accession)
D['raw_top_and_in_results_zip'] = D.accession.isin(a2.accession)
D['has_raw'] = D.accession.isin(raw_ds.accession)
dupA = []
for col in ['raw_loose_and_zipped_top', 'raw_top_and_in_results_zip']:
    t = rate_table(D[D.has_raw], 'release_year', col, min_n=30); t.insert(0, 'measure', col); dupA.append(t)
pd.concat(dupA).to_csv(OUT / 'H7a-raw-loose-and-zipped-by-year.csv', index=False)
a1.sample(min(40, len(a1)), random_state=SEED).to_csv(OUT / 'H7a-spot-check-loose-vs-zipped-top-level.csv', index=False)
a2.sample(min(40, len(a2)), random_state=SEED).to_csv(OUT / 'H7a-spot-check-raw-inside-results-zips.csv', index=False)
a1_ds = a1.groupby('accession').agg(runs=('run', 'size'), zip_TB=('zip_sz', lambda s: s.sum() / 1e12), median_size_ratio=('ratio', 'median'))
a1_ds.sort_values('zip_TB', ascending=False).head(30).join(D.set_index('accession')[['release_year']]).to_csv(OUT / 'H7a-largest-loose-and-zipped-datasets.csv')
S['H7a_raw_loose_and_zipped'] = {
    'datasets_with_raw': int(D.has_raw.sum()),
    'top_level': {'datasets': int(a1.accession.nunique()), 'groups': int(a1.grp.nunique()), 'runs': int(len(a1)),
                  'zipped_copy_TB': round(a1.zip_sz.sum() / 1e12, 3),
                  'zipped_copies_with_size_ratio_0.5_to_1.05_TB': round(a1.zip_sz[a1.ratio.between(.5, 1.05)].sum() / 1e12, 3), 'loose_copy_TB': round(a1.loose_sz.sum() / 1e12, 2),
                  'median_zip_to_loose_size_ratio': round(float(a1.ratio.median()), 3),
                  'ratio_between_0.5_and_1.05_pct': round(100 * a1.ratio.between(.5, 1.05).mean(), 1),
                  'pct_of_datasets_with_raw': round(100 * D.raw_loose_and_zipped_top[D.has_raw].mean(), 2),
                  'ci': ci(D.raw_loose_and_zipped_top[D.has_raw], D.grp[D.has_raw])},
    'inside_results_zips': {'datasets': int(a2.accession.nunique()), 'runs': int(len(a2)),
                            'member_TB_uncompressed': round(a2.usize.sum() / 1e12, 2),
                            'runs_with_loose_copy': int(a2.loose_sz.notna().sum()),
                            'loose_copy_size_within_2pct_of_member_pct': round(100 * a2.size_match_2pct[a2.loose_sz.notna()].mean(), 1),
                            'verified_duplicate_TB_(size within 2pct)': round(a2.usize[a2.size_match_2pct].sum() / 1e12, 2),
                            'verified_duplicate_datasets': int(a2[a2.size_match_2pct].accession.nunique()),
                            'verified_duplicate_groups': int(a2[a2.size_match_2pct].grp.nunique()),
                            'runs_matched_only_to_a_top_level_archive': int((a2.loose_sz.isna()).sum()),
                            'raw_members_in_results_zips_total_TB': round(raw_mem.usize.sum() / 1e12, 2),
                            'raw_member_runs_total': int(len(raw_mem)), 'raw_member_datasets_total': int(raw_mem.accession.nunique())}}

# (b) identical checksums among top-level files (checksum window only)
c.execute(f"""create temp table ckf as
  select tf.accession, file_name, layer, sz, checksum, d.release_year, d.release_date, d.pi_group
  from tf join ds d using (accession)
  where checksum <> '' and sz > 0 and d.release_year between {WINDOW[0]} and {WINDOW[1]}""")
den = c.execute("select layer, count(*) files, sum(sz) bytes from ckf group by 1").df().set_index('layer')
# within dataset: k copies of one checksum -> (k-1) redundant copies
within = c.execute("""select layer, count(*) dup_groups, sum(k - 1) redundant_files, sum((k - 1) * sz) redundant_bytes,
                             count(distinct accession) datasets
  from (select accession, checksum, min(layer) layer, max(sz) sz, count(*) k from ckf group by 1, 2 having count(*) > 1) group by 1""").df().set_index('layer')
# across datasets: first dataset (release date, accession) holds the original; later datasets' copies are redundant
across = c.execute("""
  with per as (select checksum, accession, min(layer) layer, max(sz) sz, min(pi_group) g, min(release_date) rdate
               from ckf group by 1, 2),
  o as (select *, first_value(g) over w g0, first_value(accession) over w a0, row_number() over w rn
        from per window w as (partition by checksum order by rdate, accession))
  select layer, case when g is null or g0 is null then 'unknown group' when g = g0 then 'same group' else 'different group' end relation,
         count(*) redundant_copies, sum(sz) redundant_bytes, count(distinct accession) later_datasets, count(distinct checksum) checksums
  from o where rn > 1 group by 1, 2""").df()
dupB = den.sort_index().join(within.add_prefix('within_'), how='left').fillna(0)
for rel in ['same group', 'different group', 'unknown group']:
    s = across[across.relation == rel].set_index('layer')
    dupB['across_' + rel.replace(' ', '_') + '_bytes'] = s.redundant_bytes.reindex(dupB.index).fillna(0)
    dupB['across_' + rel.replace(' ', '_') + '_copies'] = s.redundant_copies.reindex(dupB.index).fillna(0)
dupB['TB_checksummed'] = dupB.bytes / 1e12
for col in [x for x in dupB.columns if x.endswith('_bytes')]:
    dupB['pct_' + col] = (100 * dupB[col] / dupB.bytes).round(2)
dupB.loc['ALL'] = dupB.sum(numeric_only=True)
for col in [x for x in dupB.columns if x.startswith('pct_')]:
    dupB.loc['ALL', col] = round(100 * dupB.loc['ALL', col[4:]] / dupB.loc['ALL', 'bytes'], 2)
dupB.to_csv(OUT / 'H7b-checksum-duplicates-by-layer.csv')
# near-duplicate datasets: share of a later dataset's checksummed bytes already present in an earlier dataset
nd = c.execute("""
  with per as (select checksum, accession, max(sz) sz, min(pi_group) g, min(release_date) rdate from ckf group by 1, 2),
  o as (select *, first_value(g) over w g0, row_number() over w rn from per window w as (partition by checksum order by rdate, accession))
  select accession, sum(sz) bytes, sum(case when rn > 1 then sz else 0 end) seen_bytes,
         sum(case when rn > 1 and g = g0 then sz else 0 end) seen_same, sum(case when rn > 1 and g <> g0 then sz else 0 end) seen_other
  from o group by 1""").df()
nd['share_seen'] = nd.seen_bytes / nd.bytes
nd = nd.merge(D[['accession', 'release_year', 'grp']], on='accession')
nd['near_dup_90'] = nd.share_seen >= .9
nd['near_dup_90_same'] = nd.seen_same / nd.bytes >= .9
nd['near_dup_90_other'] = nd.seen_other / nd.bytes >= .9
ndy = rate_table(nd, 'release_year', 'near_dup_90', min_n=30); ndy.to_csv(OUT / 'H7b-near-duplicate-datasets-by-year.csv', index=False)
nd[nd.near_dup_90].sort_values(['bytes', 'accession'], ascending=[False, True]).head(40).to_csv(OUT / 'H7b-near-duplicate-datasets-largest.csv', index=False)
topck = c.execute("""select checksum, min(layer) layer, max(sz) sz, count(*) files, count(distinct accession) datasets,
                     count(distinct pi_group) n_groups, (count(*) - 1) * max(sz) / 1e9 redundant_GB, min(file_name) example
                     from ckf group by 1 having count(*) > 1 order by redundant_GB desc, checksum limit 40""").df()
topck.to_csv(OUT / 'H7b-checksums-most-redundant-bytes.csv', index=False)
ay = c.execute("""
  with per as (select checksum, accession, max(sz) sz, min(pi_group) g, min(release_date) rdate, min(release_year) y
               from ckf group by 1, 2),
  o as (select *, row_number() over (partition by checksum order by rdate, accession) rn from per)
  select o.y release_year, sum(sz) across_denominator_bytes, sum(case when rn > 1 then sz else 0 end) across_bytes
  from o group by 1 order by 1""").df()
wy = c.execute("""select release_year, sum(sz) bytes, sum(wb) within_bytes from (
    select accession, checksum, min(release_year) release_year, max(sz) * count(*) sz, (count(*) - 1) * max(sz) wb
    from ckf group by 1, 2) group by 1""").df()
ay = ay.merge(wy, on='release_year')
ay['pct_across'] = (100 * ay.across_bytes / ay.bytes).round(2); ay['pct_within'] = (100 * ay.within_bytes / ay.bytes).round(2)
ay.to_csv(OUT / 'H7b-checksum-duplicate-bytes-by-year.csv', index=False)
S['H7b_checksum_duplicates'] = {
    'window_release_years': WINDOW, 'checksummed_files': int(den.files.sum()), 'checksummed_TB': round(den.bytes.sum() / 1e12, 1),
    'caveat_trend': 'a dataset can only match earlier checksummed datasets (window starts 2016), so across-dataset shares rise with year partly by construction',
    'within_dataset_redundant_TB': round(within.redundant_bytes.sum() / 1e12, 2), 'within_dataset_datasets': int(within.datasets.sum()),
    'across_redundant_TB_by_relation': (across.groupby('relation').redundant_bytes.sum() / 1e12).round(2).to_dict(),
    'across_later_datasets_by_relation': across.groupby('relation').later_datasets.sum().to_dict(),
    'near_duplicate_datasets_90pct_bytes_seen': {'all': int(nd.near_dup_90.sum()), 'same_group': int(nd.near_dup_90_same.sum()),
                                                 'other_group': int(nd.near_dup_90_other.sum()), 'of_datasets': int(len(nd)),
                                                 'TB': round(nd[nd.near_dup_90].bytes.sum() / 1e12, 2)}}

# (c1) same member (base name + exact size >= 1 MB) in two or more zips of one dataset
c1 = c.execute("""select accession, bn, usize, mrole as "role", count(distinct url) zips from mbig where usize >= 1e6
                  group by 1, 2, 3, 4 having count(distinct url) > 1""").df()
c1['redundant'] = (c1.zips - 1) * c1.usize
c1r = c1.groupby('role').agg(member_sets=('bn', 'size'), datasets=('accession', 'nunique'), redundant_TB=('redundant', lambda s: s.sum() / 1e12))
c1r.round(3).sort_values('redundant_TB', ascending=False).to_csv(OUT / 'H7c-same-member-in-several-zips-of-a-dataset.csv')
# (c2) loose top-level file repeated inside a zip of the same dataset (base name + exact size, >= 1 MB)
c2 = c.execute("""
  with t as (select accession, lower(file_name) bn, sz, layer from tf where not is_archive and sz >= 1e6),
       m as (select distinct accession, bn, usize, mrole from mbig where usize >= 1e6)
  select t.accession, t.layer, m.mrole as "role", t.bn, t.sz from t join m on t.accession = m.accession and t.bn = m.bn and t.sz = m.usize order by all""").df()
c2r = c2.groupby('layer').agg(files=('bn', 'size'), datasets=('accession', 'nunique'), TB=('sz', lambda s: s.sum() / 1e12))
c2r.round(3).sort_values('TB', ascending=False).to_csv(OUT / 'H7c-loose-file-repeated-inside-zip.csv')
c2.sample(min(30, len(c2)), random_state=SEED).to_csv(OUT / 'H7c-spot-check-loose-repeated-inside-zip.csv', index=False)
S['H7c_other_redundancy'] = {
    'same_member_in_2plus_zips': {'datasets': int(c1.accession.nunique()), 'redundant_TB_uncompressed': round(c1.redundant.sum() / 1e12, 2)},
    'loose_file_repeated_inside_zip': {'datasets': int(c2.accession.nunique()), 'files': int(len(c2)),
                                       'TB': round(c2.sz.sum() / 1e12, 2)}}

(OUT / 'Summary.json').write_text(json.dumps(S, indent=2, default=str))
shutil.rmtree(TMP, ignore_errors=True)
print(json.dumps(S, indent=1, default=str))
