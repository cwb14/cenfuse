"""The run manifest (OUT/run_info.json): what made these outputs, so a result can be traced and a rerun compared.

cenfuse's version and git commit, the command, host, threads and times, Python and the package versions, every external tool's
path and version string, and a content fingerprint of every input (runtime.file_fingerprint). Written when a run starts (status
'running') and again when it ends ('done', or 'failed' with the error), always atomically. Nothing in it decides anything."""
from __future__ import annotations

import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path

from . import __version__
from .runtime import atomic_write_text, fasta_md5, file_fingerprint

PACKAGES = ("numpy", "scipy", "pandas", "pysam", "pyhmmer", "pyabpoa", "matplotlib", "h5py")
# external tools and the arguments that make each print its version (the first line with a version-like token is kept)
TOOLS = {"FasTAN": [], "FAtoGDB": [], "ONEview": [], "trf": [], "makeblastdb": ["-version"], "blastn": ["-version"],
         "dustmasker": ["-version"], "mafft": ["--version"], "cd-hit-est": ["-h"], "minimap2": ["--version"],
         "cmsearch": ["-h"], "samtools": ["--version"]}
VERSION_RE = re.compile(r"\d+\.\d+[\w.\-+()/ ]*")
NAME = "run_info.json"


_CONDA = {}


def conda_package(path: str) -> str | None:
    """'name version build' of the conda package that installed the file at path (from the environment's conda-meta
    records), or None outside a conda environment: FasTAN and TRF print no version of their own."""
    env = Path(path).resolve().parent.parent
    meta = env / "conda-meta"
    if not meta.is_dir():
        return None
    if env not in _CONDA:
        owner = {}
        for f in meta.glob("*.json"):
            try:
                rec = json.loads(f.read_text())
            except (OSError, ValueError):
                continue
            for x in rec.get("files", ()):
                if x.startswith("bin/"):
                    owner[x[4:]] = f"{rec.get('name')} {rec.get('version')} {rec.get('build')}"
        _CONDA[env] = owner
    return _CONDA[env].get(Path(path).name)


def tool_version(tool: str, args: list) -> dict:
    """{'path', 'version', 'conda'} of one tool on PATH ('version' the first output line naming one, 'unknown' when none;
    'conda' the package that installed it), or {'path': None} when the tool is absent. Each probe times out after 10 s."""
    path = shutil.which(tool)
    if path is None:
        return dict(path=None, version=None, conda=None)
    try:
        r = subprocess.run([path, *args], capture_output=True, text=True, timeout=10, stdin=subprocess.DEVNULL)
        text = (r.stdout or "") + "\n" + (r.stderr or "")
    except (OSError, subprocess.SubprocessError):
        text = ""
    line = next((x.strip() for x in text.splitlines() if VERSION_RE.search(x) and
                 re.search(r"version|v\d|\d+\.\d+", x, re.I)), "")
    return dict(path=path, version=line[:160] or "unknown", conda=conda_package(path))


def git_state() -> dict:
    """The commit and dirty state of the cenfuse source tree, when it is a git checkout (an editable install); {} otherwise."""
    root = Path(__file__).resolve().parent.parent
    if not (root / ".git").exists():
        return {}
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0")    # read only: never take .git/index.lock beside a running checkout
    try:
        head = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=10, env=env).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--untracked-files=normal", "--", "cenfuse"],
                               capture_output=True, text=True, timeout=10, env=env).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return {}
    return dict(commit=head or None, dirty=bool(dirty), tree=str(root))


def packages() -> dict:
    out = {}
    for p in PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    return out


def _fp(path) -> dict | None:
    if not path:
        return None
    p = Path(path)
    try:
        return dict(path=str(p.resolve()), fingerprint=file_fingerprint(p), bytes=p.stat().st_size)
    except OSError as e:
        return dict(path=str(p), error=str(e))


def _meth_files(table) -> list:
    """The data files a methylation samples table names (CG, CHG, CHH columns), as the pipeline resolves them."""
    try:
        import pandas as pd
        t = pd.read_csv(table, sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")   # as meth.io
    except (OSError, ValueError):
        return []
    return [str(Path(v) if Path(v).is_absolute() else Path(table).parent / v)
            for col in ("CG", "CHG", "CHH") if col in t for v in t[col] if v]


def inputs(genomes) -> dict:
    """Per genome: the FASTA (its whole md5, as the stage keys use it), its index and every optional input, each with its
    content fingerprint."""
    out = {}
    for g in genomes:
        try:
            md5 = fasta_md5(g.fasta)
        except OSError as e:
            md5 = f"error: {e}"
        d = dict(fasta=str(Path(g.fasta).resolve()), fasta_md5=md5, fai=_fp(f"{g.fasta}.fai"), chromosomes=len(g.chroms),
                 chromosome_set=g.chrom_source, tes=_fp(g.tes), genes=_fp(g.genes), meth=_fp(g.meth),
                 meth_files=[_fp(f) for f in _meth_files(g.meth)] if g.meth else [],
                 hic=[dict(name=lib.name, files=[_fp(f) for f in lib.files]) for lib in g.hic])
        out[g.name] = d
    return out


def start(out: Path, run) -> dict:
    """The manifest at the start of a run (written at once)."""
    info = dict(cenfuse=__version__, git=git_state(), status="running", command=sys.argv, cwd=os.getcwd(),
                host=socket.gethostname(), pid=os.getpid(), threads=run.threads, started=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                genome_dir=str(run.gdir), out=str(run.out), satellite_dir=str(run.sat_import) if run.sat_import else None,
                satellite_dir_cenfuse=(read(run.sat_import) or {}).get("cenfuse") if run.sat_import else None,
                force=run.force, python=sys.version.split()[0], platform=platform.platform(), packages=packages(),
                tools={t: tool_version(t, a) for t, a in TOOLS.items()}, inputs=inputs(run.genomes))
    write(out, info)
    return info


def finish(out: Path, info: dict, status: str, t0: float, error: str | None = None, incomplete=(), cost: dict | None = None) -> None:
    """The manifest at the end of a run; cost: runtime.Usage.cost() (wall time, CPU hours, peak memory of every process)."""
    info = dict(info, status=status, ended=time.strftime("%Y-%m-%dT%H:%M:%S%z"), minutes=round((time.time() - t0) / 60, 2))
    if cost:
        info["cost"] = cost
    if error:
        info["error"] = error
    if incomplete:
        info["incomplete"] = list(incomplete)
    write(out, info)


def write(out: Path, info: dict) -> None:
    atomic_write_text(Path(out) / NAME, json.dumps(info, indent=1, default=str) + "\n")


def read(folder) -> dict | None:
    """The manifest of an earlier run's output folder (or of the folder above a satellite folder), or None; never this run's
    own (a --satellite-dir inside the output folder)."""
    for d in (Path(folder), Path(folder).parent):
        f = d / NAME
        if f.is_file():
            try:
                info = json.loads(f.read_text())
            except (OSError, ValueError):
                return None
            if info.get("host") == socket.gethostname() and info.get("pid") == os.getpid():
                return None
            return info
    return None
