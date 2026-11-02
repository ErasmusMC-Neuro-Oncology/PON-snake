# ---------------------------------------------------------------------------------
# nas.smk -- read raw data from the NAS. Include after `configfile:`:
#     include: "nas.smk"
#
#   nas("/abs/path")   file on the NAS (absolute SFTP path)      -> input: sections
#   nas_local(path)    the local path nas(path) is retrieved to  -> for samplesheets
#   nas_fastqs         input function: this sample's fastqs on the NAS
#   PATIENTS           the samples this project selected
#
# This pipeline writes NOTHING to the NAS: every output stays local. The fastqs
# are declared as nas() inputs of the Sarek rule, so Snakemake retrieves them
# when that job starts and deletes the local copies when it finishes.
#
# Needs in config:
#   all:
#     samplesheet_nas: /path/to/Samplesheet_NAS.tsv   # the shared NAS samplesheet
#     sample_regex: 'Normal'                          # which of its samples are ours
#   nas_storage:
#     host: neuro-genomic-1.erasmusmc.nl
#     user: <account>                                 # needs READ on the raw data
#     password_file: ~/.passwd/<file>                 # chmod 600
# ---------------------------------------------------------------------------------
import os
import sys
import pandas as pd

#++++++++++++++++++++++++++++++++++ 1 PATCH THE SFTP PLUGIN +++++++++++++++++++++++++++++++++++++++++
# The plugin caches ONE pysftp connection per host and never validates or refreshes
# it, and each StorageObject caches it again in self.conn. After a long job the NAS
# has dropped the idle session, and then retrieving fails with OSError while storing
# fails with RecursionError (pysftp's isdir() swallows IOError -> makedirs() recurses
# on '/'). Fix: TCP keepalives + revalidate before every operation.
import snakemake_storage_plugin_sftp as _sftp_plugin

if not getattr(_sftp_plugin, "_nas_patched", False):
    _orig_get_conn = _sftp_plugin.StorageProvider.get_conn

    def _get_conn(self, hostname, port=None):
        port = port or 22                      # plugin passes None and crashes on it
        key = (hostname, port)
        conn = self.conn_pool.get(key)
        if conn is not None:
            try:
                conn.pwd                       # cheap liveness probe
            except Exception:
                print(f"[nas] sftp connection to {hostname} went stale, reconnecting",
                      file=sys.stderr, flush=True)
                try:
                    conn.close()
                except Exception:
                    pass
                self.conn_pool.pop(key, None)
                conn = None
        if conn is None:
            conn = _orig_get_conn(self, hostname, port)
            try:
                conn._transport.set_keepalive(30)
            except Exception:
                pass
        return conn

    _sftp_plugin.StorageProvider.get_conn = _get_conn

    def _refresh_conn_first(orig):
        def wrapper(self, *args, **kwargs):
            self.conn = self.provider.get_conn(self.parsed_query.hostname,
                                               self.parsed_query.port or 22)
            return orig(self, *args, **kwargs)
        return wrapper

    for _m in ("store_object", "retrieve_object", "exists", "mtime", "size", "remove"):
        if hasattr(_sftp_plugin.StorageObject, _m):
            setattr(_sftp_plugin.StorageObject, _m,
                    _refresh_conn_first(getattr(_sftp_plugin.StorageObject, _m)))

    _sftp_plugin._nas_patched = True

#++++++++++++++++++++++++++++++++++ 2 NAS STORAGE ++++++++++++++++++++++++++++++++++++++++++++++++++++

def _read_nas_password(path):
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"NAS password file not found: {path}")
    if os.stat(path).st_mode & 0o077:
        raise PermissionError(f"{path} is readable by others; run: chmod 600 {path}")
    with open(path) as fh:
        return fh.read().strip()

_NAS_HOST     = config["nas_storage"]["host"]
_NAS_USER     = config["nas_storage"]["user"]
_NAS_PASSWORD = _read_nas_password(config["nas_storage"]["password_file"])

storage nas:
    provider="sftp",
    username=_NAS_USER,
    password=_NAS_PASSWORD,
    max_requests_per_second=20


def nas(path):
    """File on the NAS, by its absolute path.
    NB: no port in the URL -- a ':' in the path breaks singularity -B binds."""
    return storage.nas(f"sftp://{_NAS_HOST}/{str(path).lstrip('/')}")


# Where Snakemake retrieves nas(path) to: local_prefix / <host><path>
# (StorageObject.local_path() = provider.local_prefix / local_suffix()).
# ABSOLUTE: local_storage_prefix defaults to the relative '.snakemake/storage',
# and a relative path here would end up in relative (dangling) symlinks.
NAS_LOCAL_PREFIX = os.path.abspath(
    os.path.join(str(workflow.storage_settings.local_storage_prefix), "nas", _NAS_HOST))


def nas_local(path):
    """Local path that nas(path) is retrieved to, as a plain string."""
    return os.path.join(NAS_LOCAL_PREFIX, str(path).lstrip('/'))

#++++++++++++++++++++++++++++++++++ 3 SAMPLES FROM THE NAS SAMPLESHEET +++++++++++++++++++++++++++++++
# The shared samplesheet lists every project; `sample_regex` picks this one's
# normals. R1_DNA/R2_DNA hold ';'-separated paths relative to the NAS root (one
# pair per lane), which is all we need -- no directory listing.

_NAS_TSV   = config['all']['samplesheet_nas']
_NAS_REGEX = config['all']['sample_regex']


def _split_list(value):
    """';'-separated NAS paths -> list of clean, root-relative paths ('NA' -> [])."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [p.strip().lstrip('/') for p in str(value).split(';')
            if p.strip() and p.strip() != 'NA']


NAS_FASTQ = {}       # patient -> files to retrieve (fastq pairs, or bam + bai)
NAS_KIND  = {}       # patient -> 'fastq' | 'bam'

if os.path.isfile(_NAS_TSV):
    _tbl = pd.read_csv(_NAS_TSV, sep='\t')
    _keep = _tbl['Sample'].astype(str).str.contains(_NAS_REGEX, regex=True, na=False)
    _tbl = _tbl[_keep]
    if _tbl.empty:
        raise ValueError(f"sample_regex {_NAS_REGEX!r} matched no samples in {_NAS_TSV}")
    for _, _row in _tbl.iterrows():
        _p  = str(_row['Sample'])
        _r1 = _split_list(_row.get('R1_DNA'))
        _r2 = _split_list(_row.get('R2_DNA'))

        if not _r1:
            continue                            # no sequencing data for this sample

        if _r1[0].lower().endswith('.bam'):
            # pre-aligned BAM -> sarek --step markduplicates; no R2, index alongside
            if len(_r1) != 1:
                raise ValueError(f"expected exactly 1 BAM for {_p}, found {len(_r1)}")
            if _r2:
                raise ValueError(f"R2_DNA must be empty for BAM input ({_p})")
            NAS_FASTQ[_p] = _r1 + [p + '.bai' for p in _r1]
            NAS_KIND[_p]  = 'bam'
        else:
            if len(_r1) != len(_r2):
                raise ValueError(f"R1/R2 fastq counts differ for {_p}: {len(_r1)} vs {len(_r2)}")
            NAS_FASTQ[_p] = _r1 + _r2
            NAS_KIND[_p]  = 'fastq'

PATIENTS = sorted(NAS_FASTQ)


def nas_fastqs(wildcards):
    """This sample's raw files on the NAS (fastqs, or bam + bai). Declaring these
    makes Snakemake retrieve them before the job and delete the local copies once
    it has finished."""
    paths = NAS_FASTQ.get(wildcards.sample)
    if not paths:
        raise ValueError(f"no NAS input known for sample '{wildcards.sample}'")
    return [nas(p) for p in paths]


def sarek_step(wildcards):
    """Each sample is its own nextflow run, so fastq and BAM samples can be mixed:
    fastq starts at mapping, a pre-aligned BAM at markduplicates."""
    return 'markduplicates' if NAS_KIND[wildcards.sample] == 'bam' else 'mapping'


def sarek_save_mapped(wildcards):
    """--save_mapped only applies when sarek does the mapping."""
    return '--save_mapped' if NAS_KIND[wildcards.sample] == 'fastq' else ''

