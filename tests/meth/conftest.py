import sys

from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

import meth_synth as synth  # noqa: E402

@pytest.fixture(scope="session")
def syn(tmp_path_factory):
    """The synthetic genome of tests/synth.py (built once per test session)."""
    return synth.make(tmp_path_factory.mktemp("syn"))

@pytest.fixture(scope="session")
def prep(syn):
    """Synthetic inputs read, classed, binned and depth-filtered once: dict(samples, chroms, cmap, bins, geo, D)."""
    from cenfuse.meth import bins, classes, io, track
    fai = io.read_fai(syn["fai"])
    names, G, D = track.chromosome_set(fai)
    chroms = {c: fai[c] for c in names}
    samples = io.read_samples(syn["samples"])
    cmap = classes.build([io.read_satellites(syn["satellites"]), io.read_ltrquest(syn["ltrquest"]),
                          io.read_genes(syn["genes"])], names)
    raw = bins.collect(samples, chroms, cmap, "test", None, 1)
    filt, stats = bins.depth_filter(raw, samples, chroms, io.read_exclusions(syn["exclude"]))
    return dict(samples=samples, chroms=chroms, cmap=cmap, bins=filt, stats=stats, geo=track.geometry(D), D=D)

from meth_common import run_scan  # noqa: E402,F401
