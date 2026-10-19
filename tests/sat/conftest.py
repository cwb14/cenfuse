import shutil

import pytest

from cenfuse.sat.params import Params

from sat_synth import build_complex

@pytest.fixture(scope="session")
def test_params():
    return Params(**TEST_PARAMS)

@pytest.fixture(scope="session")
def synth_complex(tmp_path_factory):
    cdir = tmp_path_factory.mktemp("synth") / "Synth"
    truth = build_complex(cdir)
    return cdir, truth

@pytest.fixture(scope="session")
def synth_seed_cache(synth_complex, test_params, tmp_path_factory):
    from cenfuse.sat.genomes import discover_genomes
    from cenfuse.sat.seeds import seed_genomes
    cache = tmp_path_factory.mktemp("seedcache")
    seed_genomes(discover_genomes(synth_complex[0]), cache, 4, test_params)
    return cache

def _load(cdir, cache, p):
    from cenfuse.sat.genomes import discover_genomes
    from cenfuse.sat.seeds import load_seeds, seed_genomes
    return {g: load_seeds(path) for g, path in seed_genomes(discover_genomes(cdir), cache, 4, p).items()}

@pytest.fixture(scope="session")
def synth_seeds(synth_complex, synth_seed_cache, test_params):
    """Seeds of the default engine (FasTAN)."""
    return _load(synth_complex[0], synth_seed_cache, test_params)

@pytest.fixture(scope="session")
def synth_trf(synth_complex, test_params, tmp_path_factory):
    return _load(synth_complex[0], tmp_path_factory.mktemp("trfcache"), Params(**dict(TEST_PARAMS, seed_engine="trf")))

@pytest.fixture(scope="session")
def synth_families(synth_complex, synth_seeds, test_params, tmp_path_factory):
    from cenfuse.sat.families import build_families
    from cenfuse.sat.genomes import discover_genomes
    wd = tmp_path_factory.mktemp("families")
    fam = build_families(synth_seeds, discover_genomes(synth_complex[0]), wd, 4, test_params)
    return fam, wd

@pytest.fixture(scope="session")
def synth_annot(synth_complex, synth_families, test_params, tmp_path_factory):
    import pyhmmer
    from cenfuse.sat.annotate import annotate_genome
    from cenfuse.sat.genomes import discover_genomes
    fam, wd = synth_families
    with pyhmmer.plan7.HMMFile(wd / "families.hmm") as f:
        hmms = list(f)
    g = [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]
    out = tmp_path_factory.mktemp("annot")
    return annotate_genome(g, hmms, dict(zip(fam.fam_id, fam.consensus_len)), 4, test_params,
                           out / "Pa_A.monomers.tsv.gz", out / "Pa_A.arrays.tsv")

from sat_common import TEST_PARAMS, have  # noqa: E402,F401
