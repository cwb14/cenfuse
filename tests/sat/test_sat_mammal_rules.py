"""Rules added after the T2T mammal benchmark (docs/design.md, 'Satellites and TEs'): dispersed repeats
are not satellites, dominance is judged per satellite kind, acrocentric centromeres are not subtelomeres, classify's
per-chromosome choice makes the call, and a competing cluster is reported and lowers the call to low."""
import pandas as pd
import pysam

from cenfuse.sat.annotate import ARR_COLS
from cenfuse.sat.centromere import finalize, satellite_call
from cenfuse.sat.classify import binom_cdf, choose_centromeric, family_metrics, kinds, tandem_support
from cenfuse.sat.params import Params

P = Params()


def arr(fam, chrom, start, n, mlen=180):
    return (0, chrom, start, start + n * mlen, "+", fam, n, 0.9, mlen)


def mono(fam, chrom, start, n, mlen, cov):
    """n monomer hits of a mlen-bp model, each covering a fraction cov of the model, head to tail."""
    w = int(mlen * cov)
    return [dict(family=fam, chrom=chrom, start=start + i * w, end=start + (i + 1) * w, hmm_from=1, hmm_to=w,
                 model_len=mlen) for i in range(n)]


def test_tandem_support_separates_satellites_from_interspersed_chains():
    arrays = pd.DataFrame([arr("SAT", "c1", 1_000_000, 500), arr("INT", "c1", 5_000_000, 20, 1500)],
                          columns=ARR_COLS)
    seeds = pd.DataFrame(dict(chrom=["c1", "c1"], start=[1_000_000, 5_000_000], end=[1_090_000, 5_001_000],
                              period=[180, 40]))
    mons = pd.DataFrame(mono("SAT", "c1", 1_000_000, 500, 180, 1.0) + mono("INT", "c1", 5_000_000, 200, 1500, 0.1))
    t = tandem_support(arrays, mons, seeds, P)
    assert t["SAT"] == (1.0, 1.0)
    assert round(t["INT"][0], 3) == round(1000 / 30_000, 3) and t["INT"][1] == 0.0


def test_dispersed_family_is_never_eligible_and_none_is_without_inputs():
    rows = [arr("INT", f"c{i}", j * 400_000, 3, 1500) for i in range(1, 5) for j in range(25)]   # chained along
    a = pd.DataFrame(rows, columns=ARR_COLS)
    chroms = {f"c{i}": 10_000_000 for i in range(1, 5)}
    fams, _ = family_metrics(a, chroms, None, P, tandem={"INT": (0.05, 0.003)})
    assert not fams.eligible.iloc[0] and "dispersed repeat" in fams.reason.iloc[0]
    fams, _ = family_metrics(a, chroms, None, P, tandem={"INT": (0.05, 0.95)})          # whole units: tandem
    assert fams.eligible.iloc[0]
    fams, _ = family_metrics(a, chroms, None, P)                                         # no seeds: no test
    assert fams.eligible.iloc[0]
    conc = pd.DataFrame([arr("INT", f"c{i}", 3_000_000 + j * 30_000, 3, 1500) for i in range(1, 5)
                         for j in range(10)], columns=ARR_COLS)        # not tandem, but at one locus (TE-derived)
    fams, _ = family_metrics(conc, chroms, None, P, tandem={"INT": (0.05, 0.003)})
    assert fams.eligible.iloc[0] and fams.locus_fold.iloc[0] > 2


def test_dominance_is_judged_per_satellite_kind():
    """Alpha split into 3 subfamilies (one kind) is the largest kind on 6/8 chromosomes; a pericentric 43-bp
    satellite has the largest single-family cluster on 5/8. Per family, the 43-bp family is dominant; per kind,
    alpha is."""
    chroms = {f"c{i}": 50_000_000 for i in range(1, 9)}
    rows = []
    for i in range(1, 9):
        c = f"c{i}"
        for j, f in enumerate(("A1", "A2", "A3")):                   # alpha: three adjacent subfamily arrays
            rows.append(arr(f, c, 20_000_000 + j * 200_000, 1000 if i <= 6 else 300, 171))
        rows.append(arr("H", c, 22_000_000, 1600 if i <= 5 else 200, 43 * 4))   # pericentric, 1.6 Mb away
    info = {f: dict(superfamily="A", length=171, rdna=".", chain_group=f) for f in ("A1", "A2", "A3")}
    info["H"] = dict(superfamily="H", length=43, rdna=".", chain_group="H")
    a = pd.DataFrame(rows, columns=ARR_COLS)
    fams, mains = family_metrics(a, chroms, None, P, info)
    old, _ = choose_centromeric(fams, mains, chroms, P, info)
    new, per_chrom = choose_centromeric(fams, mains, chroms, P, info, a)
    assert set(old.family[old.dominant]) == {"H"}
    assert set(new.family[new.dominant]) == {"A1", "A2", "A3"}
    assert all(per_chrom[c][0] in ("A1", "A2", "A3") for c in chroms)
    assert kinds(["A1", "H", "A3", "A2"], info) == [["H"], ["A1", "A2", "A3"]] or \
        sorted(map(tuple, kinds(["A1", "H", "A3", "A2"], info))) == [("A1", "A2", "A3"), ("H",)]


def acro_genome(both_ends: int, interior: bool):
    """20 chromosomes of 100 Mb (5-Mb end zones), a strong terminal family at the left end of each (and at both ends
    of `both_ends`). With 20 one-ended chromosomes P = 0.75^20 = 0.003; 10 would give 0.056 (too few to test)."""
    chroms = {f"c{i}": 100_000_000 for i in range(20)}
    rows = []
    for i in range(20):
        c = f"c{i}"
        rows.append(arr("T", c, 0, 3000, 817))                        # 2.45 Mb at the left end
        if i < both_ends:
            rows.append(arr("T", c, 100_000_000 - 1000 * 817, 1000, 817))
        if interior:
            rows.append(arr("C", c, 50_000_000, 2000, 180))
    return pd.DataFrame(rows, columns=ARR_COLS), chroms


def test_one_ended_terminal_family_is_the_acrocentric_centromere():
    a, chroms = acro_genome(0, False)
    f = family_metrics(a, chroms, None, P)[0].set_index("family")
    assert f.loc["T", "eligible"] and f.loc["T", "acrocentric"] and not f.loc["T", "subtelomeric"]
    assert "acrocentric" in f.loc["T", "reason"]


def test_ends_occupied_independently_stay_subtelomeric():
    a, chroms = acro_genome(5, False)                  # both ends on 5/20 (independent ends: 20 x 0.39 = 7.8)
    f = family_metrics(a, chroms, None, P)[0].set_index("family")
    assert not f.loc["T", "eligible"] and f.loc["T", "subtelomeric"] and not f.loc["T", "acrocentric"]


def test_a_strong_interior_family_keeps_a_terminal_family_subtelomeric():
    a, chroms = acro_genome(0, True)
    f = family_metrics(a, chroms, None, P)[0].set_index("family")
    assert f.loc["C", "eligible"] and not f.loc["T", "eligible"] and not f.loc["T", "acrocentric"]


def test_binom_cdf():
    assert abs(binom_cdf(0, 10, 0.2) - 0.8 ** 10) < 1e-12
    assert abs(binom_cdf(10, 10, 0.3) - 1.0) < 1e-9


def test_satellite_call_takes_the_cluster_of_the_chosen_family():
    a = pd.DataFrame([arr("A", "c", 5_000_000, 300), arr("H", "c", 8_000_000, 3000),
                      arr("L", "c", 5_100_000, 200, 66)], columns=ARR_COLS)          # L: a layer next to A
    assert satellite_call(a, {"A", "H", "L"}, 20_000_000, P)["main_family"] == "H"
    call = satellite_call(a, {"A", "H", "L"}, 20_000_000, P, chosen={"A"})
    assert (call["cen_start"], call["cen_end"]) == (5_000_000, 5_100_000 + 200 * 66)        # A with its layer
    assert (call["alt_start"], call["alt_end"]) == (8_000_000, 8_000_000 + 3000 * 180)
    near_h = dict(start=8_100_000, end=8_400_000)                      # a decisive CRM core at the largest cluster
    assert satellite_call(a, {"A", "H", "L"}, 20_000_000, P, chosen={"A"}, core=near_h)["main_family"] == "H"
    b = pd.DataFrame([arr("A", "c", 5_000_000, 300), arr("H", "c", 8_000_000, 3000), arr("X", "c", 15_000_000, 200)],
                     columns=ARR_COLS)                                  # core at X, neither chosen nor largest
    assert satellite_call(b, {"A", "H", "X"}, 20_000_000, P, chosen={"A"},
                          core=dict(start=15_000_000, end=15_100_000))["main_family"] == "A"


def test_competing_cluster_is_reported_and_low(tmp_path):
    a = pd.DataFrame([arr("A", "c", 5_000_000, 1000), arr("A", "c", 9_000_000, 800)], columns=ARR_COLS)
    call = satellite_call(a, {"A"}, 20_000_000, P)
    assert (call["alt_start"], call["alt_end"]) == (9_000_000, 9_000_000 + 800 * 180)
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + "A" * 20_000_000 + "\n")
    pysam.faidx(str(fa))
    with pysam.FastaFile(str(fa)) as h:
        row = finalize("c", 20_000_000, call, h, None, None, P)
    assert "competing_cluster" in row["flags"] and row["confidence"] == "low"
    one = satellite_call(a.iloc[:1], {"A"}, 20_000_000, P)
    assert "alt_start" not in one
    with pysam.FastaFile(str(fa)) as h:
        assert finalize("c", 20_000_000, one, h, None, None, P)["confidence"] == "high"
