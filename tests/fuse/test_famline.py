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


def mixed(tmp_path):
    """Three chromosomes and three unplaced sequences (in the FASTA, not in the chromosome set; one shorter than K),
    interleaved, the chromosome set in another order: random sequence with copies of one 3-kb element on both strands of
    four sequences, a tandem array of a 171-bp unit on c1 and u1, a microsatellite, N runs (one across several piece
    borders), IUPAC codes and a soft-masked stretch; arrays (unsorted, nested, one across a piece border) on c0, c1, the
    unplaced u1 and a name in no FASTA."""
    rng = np.random.default_rng(23)
    te, unit = "".join(rng.choice(NT, size=3_000)), "".join(rng.choice(NT, size=171))
    s = {c: list(rng.choice(NT, size=n)) for c, n in (("c0", 150_000), ("u1", 60_000), ("c1", 120_000), ("u2", 700),
                                                     ("c2", 90_000), ("u3", 15))}
    for c, at in (("c0", (30_000, 70_000, 141_000)), ("u1", (5_000,)), ("c1", (30_000, 95_000)), ("c2", (10_000, 45_000))):
        for i, p in enumerate(at):
            s[c][p:p + 3_000] = te if i % 2 else te.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for c, p, n in (("c1", 60_000, 60), ("u1", 20_000, 40)):
        s[c][p:p + 171 * n] = unit * n
    s["c2"][70_000:70_600] = "AC" * 300
    s["c0"][14_000:26_000] = "N" * 12_000
    s["c2"][5_000] = "N"
    for p in (1_234, 50_001, 99_999):
        s["c1"][p] = "R"
    s["c2"][40_000:52_000] = [x.lower() for x in s["c2"][40_000:52_000]]
    f = tmp_path / "m.fa"
    f.write_text("".join(f">{c}\n{''.join(x)}\n" for c, x in s.items()))
    pysam.faidx(str(f))
    fai = {c: len(x) for c, x in s.items()}
    arrays = {"c1": [(62_000, 63_000), (60_000, 60_000 + 171 * 60)], "c0": [(0, 500), (39_990, 40_010)],
              "u1": [(20_000, 20_000 + 171 * 40)], "zz": [(0, 10)]}
    return f, fai, {c: fai[c] for c in ("c2", "c0", "c1")}, arrays


def same(a, b) -> bool:
    """Equal down to the bytes: dict keys in the same order, arrays of the same dtype and shape, scalars of the same type."""
    if isinstance(a, dict):
        return type(b) is dict and list(a) == list(b) and all(same(a[k], b[k]) for k in a)
    if isinstance(a, np.ndarray):
        return type(b) is np.ndarray and a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()
    return type(a) is type(b) and a == b


def test_any_thread_count_and_piece_size_give_the_same_sample(tmp_path, monkeypatch):
    """Worker processes read the FASTA's pieces (threads); the sampled copies are the same at any thread count and piece
    size, down to the bytes: unplaced sequences, arrays, N runs and pieces ending inside a bin included."""
    from cenfuse import repline
    f, fai, chroms, arrays = mixed(tmp_path)
    ref = famline.sample(f, fai, chroms, arrays)                             # one piece per sequence, in this process
    assert ref["nk"] > 0 and np.unique(ref["chrom"]).tolist() == [0, 1, 2] and list(ref["called"]) == list(chroms)
    assert (ref["called"]["c0"][14:26] == 0).all() and (ref["called"]["c0"][:14] == 1000).all()   # the N run
    pools, pool = [], repline.process_pool

    def spy(n):
        pools.append(n)
        return pool(n)
    monkeypatch.setattr(repline, "process_pool", spy)
    for piece in (10_000, 7_919):                                             # 7,919 bp: pieces end inside bins
        monkeypatch.setattr(famline, "PIECE", piece)
        for threads in (1, 3):
            assert same(famline.sample(f, fai, chroms, arrays, threads), ref), (piece, threads)
    assert pools == [3, 3]                                                    # worker processes read the pieces


# ---------------------------------------------------------------- fuse
FAM = dict(start=400_000, end=600_000, density=310.0, status="pass", window=100_000)


def _ev(seq=None, fam=None, rep=None, hic=None, flags="."):
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence="low", flags=flags, main_family="A") if seq else \
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


def test_a_low_call_the_family_line_contradicts_beyond_resolution_moves_and_stays_low():
    """0.18 (fam_over_low, future.md #76): where the passing repeat-family window lies beyond the resolution of both readouts (the
    gap exceeds the wider interval, sat/te.agrees), a low sequence call that no other rule changed moves to the window and stays
    low; closer than that the call stays, flagged (fam_elsewhere). Not for a call in a terminal eighth (the line reads no window
    there), a confident call, a call the window touches, a line that did not pass, or a chromosome where validated Hi-C has a
    peak (it outranks an untested line)."""
    nop = "not_one_per_chromosome"                                              # the satellite's own test failed
    d = _ev(seq=(150_000, 250_000), fam=dict(FAM, start=600_000, end=800_000), flags=nop)   # gap 350 kb > 200 kb: another locus
    assert d["call"] == (600_000, 800_000) and d["confidence"] == "low" and d["rule"] == "fam_over_low"
    assert d["source"] == "composition" and "families" in d["support"] and "fam_elsewhere" not in d["flags"]
    for fl in ("no_dominant_satellite", "near_end,no_dominant_satellite"):
        assert _ev(seq=(150_000, 250_000), fam=dict(FAM, start=600_000, end=800_000), flags=fl)["rule"] == "fam_over_low"
    for fl in (".", "competing_cluster", "low_satellite,near_end"):              # low, but the satellite's own tests passed
        d = _ev(seq=(150_000, 250_000), fam=dict(FAM, start=600_000, end=800_000), flags=fl)
        assert d["call"] == (150_000, 250_000) and "fam_elsewhere" in d["flags"] and "fam_over_low" not in d["flags"], fl
    d = _ev(seq=(150_000, 250_000), fam=dict(FAM, start=300_000, end=500_000))    # gap 50 kb <= 200 kb: one region
    assert d["call"] == (150_000, 250_000) and d["confidence"] == "low" and "fam_elsewhere" in d["flags"]
    assert "families" in d["against"] and d["rule"] == "sequence" and "fam_over_low" not in d["flags"]
    for seq, fam in (((500_000, 600_000), FAM), ((20_000, 60_000), FAM), ((150_000, 250_000), dict(FAM, status="fail"))):
        d = _ev(seq=seq, fam=fam)                                       # touches; terminal eighth; no pass
        assert d["call"] == seq and not {"fam_elsewhere", "fam_over_low"} & set(d["flags"])
    s = dict(start=150_000, end=250_000, evidence="satellite", confidence="medium", flags=".", main_family="A")
    d = fuse.decide(fuse.Evidence("c", 1_000_000, s, None, None, None, lambda a, b: 0, fam=dict(FAM, start=600_000, end=800_000)))
    assert d["call"] == (150_000, 250_000) and not {"fam_elsewhere", "fam_over_low"} & set(d["flags"])   # a confident call
    hic = dict(callable=True, bin_size=HB, gate_pass=True, gate_support=False, anch_strong=False, anch_level="none",
               anch_flags="", anch_peak_start=900_000, anch_peak_end=910_000, anch_core_start=900_000, anch_core_end=910_000,
               masked=(), denovo_peak_start=-1, denovo_peak_end=-1, denovo_core_start=-1, denovo_core_end=-1,
               denovo_confidence="none", denovo_flags="")
    d = _ev(seq=(150_000, 250_000), fam=dict(FAM, start=600_000, end=800_000), hic=hic, flags=nop)   # validated Hi-C peak here
    assert d["call"] == (150_000, 250_000) and "fam_elsewhere" in d["flags"] and "fam_over_low" not in d["flags"]


def test_an_unsupported_long_unit_array_keeps_the_order_of_the_lines_that_place_a_chromosome():
    """An unsupported long-unit array is no call (fuse.decide: long_unit_unsupported), so the repeat-density peak places the
    chromosome before the family window does, as without the array: the family move does not take that call."""
    lu = dict(start=150_000, end=170_000, evidence="long_unit_tandem", confidence="low", flags=".", main_family=".")
    rep = dict(start=420_000, end=440_000, density=6.0, status="pass", window=20_000)
    for fam in (dict(FAM, start=600_000, end=800_000), dict(FAM, start=400_000, end=600_000)):
        d = fuse.decide(fuse.Evidence("c", 1_000_000, lu, None, None, None, lambda a, b: 0, rep=rep, fam=fam))
        assert d["call"] == (420_000, 440_000) and d["rule"] == "rep_only" and "long_unit_unsupported" in d["flags"]
        assert "fam_over_low" not in d["flags"]
    d = fuse.decide(fuse.Evidence("c", 1_000_000, lu, None, None, None, lambda a, b: 0, fam=dict(FAM, start=600_000, end=800_000)))
    assert d["call"] == (600_000, 800_000) and d["rule"] == "fam_only"                # no repeat-density peak: the window


def _failed_hic(**kw):
    """A Hi-C row whose gate failed on this chromosome, with no anchored peak: Hi-C says nothing unless a keyword adds a reading."""
    h = dict(callable=True, bin_size=HB, gate_pass=False, gate_support=False, gate_genome_untested=False, anch_strong=False,
             anch_level="none", anch_flags="", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1,
             masked=(), denovo_peak_start=-1, denovo_peak_end=-1, denovo_core_start=-1, denovo_core_end=-1,
             denovo_confidence="none", denovo_flags="")
    return dict(h, **kw)


def test_the_family_window_moves_only_a_satellite_call_that_no_other_line_speaks_to():
    """0.18 (fam_over_low): the window takes the place only of the satellite line's own unsure array, and only where no tested
    line speaks to that call on the chromosome (a tested reading outranks an untested one): not an LTR core held at low
    (minority_core, weak_core), not a call with centromere-like methylation, an agreeing LTR core, an agreeing tested AT trough
    or repeat-density peak, Hi-C read on the LTR cores, the cis signal or a Rabl-axis line acting here, or a CpG line; each such
    call stays, flagged fam_elsewhere. A Hi-C row whose gate failed without a reading does not block the move."""
    far = dict(FAM, start=600_000, end=800_000)
    sat = dict(start=150_000, end=250_000, evidence="satellite", confidence="low", flags="not_one_per_chromosome", main_family="A")

    def dec(seq=sat, te=None, hic=None, meth=None, **kw):
        return fuse.decide(fuse.Evidence("c", 1_000_000, seq, te, hic, meth, lambda a, b: 0, fam=far, **kw))

    assert dec(hic=_failed_hic())["rule"] == "fam_over_low"                                  # Hi-C says nothing here
    assert dec(seq=dict(sat, evidence="weak_satellite"))["rule"] == "fam_over_low"
    blocked = {
        "LTR core held at low": dict(seq=dict(sat, evidence="retrotransposon", flags="minority_core")),
        "methylation": dict(meth={"sequence_call": {"label": "centromere_like", "score": 0.9}}),
        "agreeing LTR core": dict(te=dict(clade="Tekay", start=140_000, end=260_000, decisive=False, strong=False)),
        "AT trough": dict(at=dict(start=180_000, end=220_000, gc=0.3, status="pass")),
        "repeat-density peak": dict(rep=dict(start=160_000, end=240_000, density=6.0, status="pass", window=20_000)),
        "CpG line": dict(cpg=lambda a, b: None),
        "cis signal": dict(hic=_failed_hic(cis=dict(peak=(270_000, 280_000), status="pass", tested_by="ltr", interior=False,
                                                     pair=False, trans_pair=False))),
        "Rabl axis": dict(hic=_failed_hic(), fold=dict(status="pass", point=900_000, core_start=880_000, core_end=920_000,
                                                       agrees=False)),
    }
    for what, kw in blocked.items():
        d = dec(**kw)
        assert "fam_over_low" not in d["flags"] and d["call"] != (600_000, 800_000), what
        assert "fam_elsewhere" in d["flags"], what
    alt = dict(sat, alt_start=250_000, alt_end=300_000)                                     # Hi-C on the LTR cores at the call,
    d = dec(seq=alt, hic=_failed_hic(te_pair=dict(core=(150_000, 250_000), status="pass", clade="Tekay")))   # beside its rival
    assert "hic_alt_ambiguous" in d["flags"] and "fam_over_low" not in d["flags"] and d["call"] == (150_000, 250_000)


def test_a_moved_call_reads_methylation_at_its_new_place():
    """The moved call's methylation is the provisional call's (at the window), as for a chosen competing cluster: the verdict on
    the abandoned array says nothing about the window. The call stays low: the family line never pairs."""
    meth = {"sequence_call": {"label": "dip", "score": 0.1}, "provisional": {"label": "centromere_like", "score": 0.9}}
    sat = dict(start=150_000, end=250_000, evidence="satellite", confidence="low", flags="no_dominant_satellite", main_family="A")
    d = fuse.decide(fuse.Evidence("c", 1_000_000, sat, None, None, meth, lambda a, b: 0, fam=dict(FAM, start=600_000, end=800_000)))
    assert d["rule"] == "fam_over_low" and d["meth_call"]["label"] == "centromere_like" and d["confidence"] == "low"


def test_the_family_line_speaks_only_where_its_window_holds_a_cluster():
    """0.18: every chromosome gets a best family window, also where no family clusters (the holocentric Leptidea sinapis, whose
    test passes, has 0-2 family K-mers in each window that placed or moved a call): where the window holds fewer than
    famline.MIN_COPIES family K-mers (no cluster beyond a pair), the line places no call (fam_only), moves none (fam_over_low)
    and flags none (fam_elsewhere)."""
    thin = dict(FAM, density=2.0)
    assert _ev(fam=thin)["call"] is None and _ev(fam=dict(FAM, density=3.0))["rule"] == "fam_only"
    nop = "not_one_per_chromosome"
    d = _ev(seq=(150_000, 250_000), fam=dict(thin, start=600_000, end=800_000), flags=nop)
    assert d["call"] == (150_000, 250_000) and not {"fam_over_low", "fam_elsewhere"} & set(d["flags"])
    d = _ev(seq=(150_000, 250_000), fam=dict(thin, start=300_000, end=500_000), flags=nop)
    assert not {"fam_over_low", "fam_elsewhere"} & set(d["flags"])
    assert _ev(seq=(150_000, 250_000), fam=dict(FAM, start=600_000, end=800_000, density=3.0), flags=nop)["rule"] == "fam_over_low"
