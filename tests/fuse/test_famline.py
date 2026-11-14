"""The repeat-family line (cenfuse.famline, 0.16) and its rules in fuse.decide. Synthetic FASTA only."""
import numpy as np
import pysam

from cenfuse import famline, fuse

RNG = np.random.default_rng(11)
NT = np.array(list("ACGT"))
L = 2_000_000
N_CHROMS = 12


def rnd(n):
    return "".join(RNG.choice(NT, size=n))


def genome(tmp_path, clustered=True, scattered=True, arrays_at=None):
    """N_CHROMS chromosomes of random sequence. clustered: on each, 20 copies of one 5-kb element packed into a 100-kb stretch
    at an interior place that differs between chromosomes (a family concentrated at one locus per chromosome). scattered: 20
    copies of another 5-kb element spread over each chromosome. arrays_at: {chrom: start} of a 200-copy tandem array of a
    170-bp unit (the satellite line's)."""
    cen, disp, unit = rnd(5_000), rnd(5_000), rnd(170)
    f = tmp_path / "g.fa"
    where = {}
    with open(f, "w") as o:
        for i in range(N_CHROMS):
            c = f"c{i}"
            s = list(rnd(L))
            if clustered:
                at = 500_000 + 60_000 * i
                for j in range(20):
                    p = at + j * 5_000
                    s[p:p + 5_000] = cen
                where[c] = (at, at + 100_000)
            if scattered:
                for p in RNG.choice(np.arange(10, 390) * 5_000, size=20, replace=False):
                    s[p:p + 5_000] = disp
            if arrays_at and c in arrays_at:
                a = arrays_at[c]
                s[a:a + 200 * 170] = unit * 200
            o.write(f">{c}\n{''.join(s)}\n")
    pysam.faidx(str(f))
    return f, {f"c{i}": L for i in range(N_CHROMS)}, where


def test_periodic_kmers_are_found():
    from cenfuse import repline
    code = lambda x: repline.LUT[np.frombuffer(x.encode(), dtype=np.uint8)]
    tandem, _ = repline.canonical_kmers(code("ACGTTACGTTACGTTACGTTACGTTACGTTA"), famline.K)    # period 5
    plain, _ = repline.canonical_kmers(code(rnd(31)), famline.K)
    assert famline.periodic(tandem).all() and not famline.periodic(plain).any()


def test_a_family_clustered_once_per_chromosome_passes_and_is_placed(tmp_path):
    f, chroms, where = genome(tmp_path)
    fam = famline.sample(f, chroms, chroms)
    r = famline.scan(fam, chroms)
    assert r["status"] == "pass" and r["p"] <= famline.FAM_ALPHA
    calls = r["scales"][r["best"]]["calls"]
    for c, (a, b) in where.items():
        assert calls[c]["start"] < b and calls[c]["end"] > a and calls[c]["gc"] > 0      # on the cluster, with votes


def test_scattered_copies_alone_do_not_pass(tmp_path):
    f, chroms, _ = genome(tmp_path, clustered=False)
    assert famline.scan(famline.sample(f, chroms, chroms), chroms)["status"] != "pass"


def test_satellite_arrays_are_left_out(tmp_path):
    """A tandem array is the satellite line's: its K-mers are dropped, so an array on every chromosome makes no family."""
    arrays = {f"c{i}": 1_200_000 + 20_000 * i for i in range(N_CHROMS)}
    f, chroms, _ = genome(tmp_path, clustered=False, arrays_at=arrays)
    masked = famline.sample(f, chroms, chroms, {c: [(a, a + 200 * 170)] for c, a in arrays.items()})
    pos = masked["pos"][masked["chrom"] == 0]
    assert not ((pos >= arrays["c0"] - famline.K) & (pos < arrays["c0"] + 200 * 170)).any()
    assert famline.scan(masked, chroms)["status"] != "pass"


# ---------------------------------------------------------------- fuse
FAM = dict(start=400_000, end=600_000, density=310.0, status="pass", window=100_000)


def _ev(seq=None, fam=None, rep=None, hic=None):
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence="low", flags=".", main_family="A") if seq else \
        dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 1_000_000, s, None, hic, None, lambda a, b: 0, rep=rep, fam=fam))


def test_a_passing_family_line_places_only_what_nothing_else_places():
    d = _ev(fam=FAM)
    assert d["call"] == (400_000, 600_000) and d["confidence"] == "low" and d["rule"] == "fam_only"
    assert "families" in d["support"] and d["source"] == "composition"
    assert _ev(fam=dict(FAM, status="fail"))["call"] is None
    d = _ev(seq=(100_000, 200_000), fam=FAM)                                   # a sequence call stands
    assert d["call"] == (100_000, 200_000) and "fam_only" not in d["flags"]
    rep = dict(start=700_000, end=720_000, density=6.0, status="pass", window=20_000)
    d = _ev(fam=FAM, rep=rep)                                                   # the repeat-density peak first
    assert d["call"] == (700_000, 720_000) and "rep_only" in d["flags"] and "fam_only" not in d["flags"]


def test_the_row_reports_the_family_line():
    e = fuse.Evidence("c", 1_000_000, dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family="."),
                      None, None, None, lambda a, b: 0, fam=FAM)
    r = fuse.row(e, fuse.decide(e))
    assert (r["fam_gate"], r["fam_start"], r["fam_end"], r["fam_votes"]) == ("pass", 400_001, 600_000, 310.0)
    assert set(fuse.FUSED_COLS) >= {"fam_gate", "fam_start", "fam_end", "fam_votes"}


HB = 10_000


def _hic_alone(peak):
    return dict(callable=True, bin_size=HB, gate_pass=False, gate_support=False, gate_genome_untested=True, anch_strong=False,
                anch_level="none", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
                denovo_peak_start=peak, denovo_peak_end=peak + HB, denovo_core_start=peak, denovo_core_end=peak + HB,
                denovo_confidence="high", denovo_flags="")


def test_a_family_peak_never_pairs_with_a_hic_call():
    """The repeat-family test is not calibrated (a family that clusters once per chromosome passes wherever its clusters lie),
    so even a passing, Hi-C-tested family peak at a Hi-C call is not a second tested line: the call keeps its place, its
    interval and its `low` confidence (0.16.0 raised it to medium as fam_pair; it never fired in the release)."""
    tested = dict(FAM, validated=True, agrees=True)
    for fam in (tested, dict(FAM, agrees=True), dict(tested, status="untested")):
        d = _ev(hic=_hic_alone(500_000), fam=fam)
        assert d["call"] == (500_000, 510_000) and d["confidence"] == "low" and d["rule"] == "hic_denovo_only"
        assert "families" not in d["support"] and not {"fam_pair", "fam_only"} & set(d["flags"])
    d = _ev(hic=_hic_alone(800_000), fam=tested)                                 # the peak lies elsewhere
    assert d["call"] == (800_000, 810_000) and d["confidence"] == "low"
    rep = dict(start=480_000, end=520_000, density=6.0, status="pass", window=20_000, validated=True, agrees=True)
    d = _ev(hic=_hic_alone(500_000), fam=tested, rep=rep)                       # the repeat-density peak pairs first
    assert d["call"] == (480_000, 520_000) and "rep_pair" in d["flags"] and "fam_pair" not in d["flags"]


def test_a_low_call_the_family_line_contradicts_is_flagged_and_stays():
    """0.16 (fam_elsewhere): where the passing repeat-family line clusters elsewhere on a chromosome with a low sequence call, the
    call stays and the caveat says so; not for a call in a terminal eighth (the line reads no window there), a confident call,
    a call the peak touches, or a line that did not pass. (Moving such calls onto the family window was tested and not
    adopted: the polyploid panel's homoeologous calls came no closer.)"""
    d = _ev(seq=(150_000, 250_000), fam=dict(FAM, start=600_000, end=800_000))    # 1-Mb chromosome: ends < 125 kb
    assert d["call"] == (150_000, 250_000) and d["confidence"] == "low" and "fam_elsewhere" in d["flags"]
    assert "families" in d["against"] and d["rule"] == "sequence"
    for seq, fam in (((500_000, 600_000), FAM), ((20_000, 60_000), FAM), ((150_000, 250_000), dict(FAM, status="fail"))):
        d = _ev(seq=seq, fam=fam)                                       # touches; terminal eighth; no pass
        assert d["call"] == seq and "fam_elsewhere" not in d["flags"]
    s = dict(start=150_000, end=250_000, evidence="satellite", confidence="medium", flags=".", main_family="A")
    d = fuse.decide(fuse.Evidence("c", 1_000_000, s, None, None, None, lambda a, b: 0, fam=dict(FAM, start=600_000, end=800_000)))
    assert d["call"] == (150_000, 250_000) and "fam_elsewhere" not in d["flags"]          # a confident call: not flagged
