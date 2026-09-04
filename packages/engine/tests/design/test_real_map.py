"""G2, against a real Addgene-style lentiviral map instead of a fixture we built.

Everything else in this lane designs into `synthetic_lenti_ef1a.gb` or
`synthetic_mcs_ef1a.gb`. Those are enough for CORRECTNESS -- they exercise the
coordinate model, the origin-spanning and reverse-strand paths, the carried-motif
partition and I9 -- and they are not enough for ANNOTATION QUALITY, because a
hand-built fixture is a fixture the gate was written against. #74 says so, and
recording G2 as passed against one would be #4's "regenerating a fixture is not a
fix" wearing a different hat.

So this file designs into `real_lenti_pFTMGW_EF177827.gb`: GenBank accession
EF177827, the pFTMGW lentiviral transfer vector of Geller et al. 2007 (Mol. Vis.
13:730-739), deposited by its authors and unedited here. 8928 bp, circular, and
messy in the specific ways a synthetic fixture is not:

* no `promoter` feature at all -- the CMV promoter that replaced U3 is mentioned
  only in prose, inside a `/note` on the 5' LTR `repeat_region`;
* an MCS annotated as a bare `misc_feature` whose only marking is the free text
  `/note="MCS; Multiple Cloning Site"`;
* real LTRs, a real Psi, a real WPRE, a real ColE1;
* real repeated DNA the vector's authors shipped and BT5 did not create;
* and an eGFP CDS carrying `/transl_table=11` -- the BACTERIAL table, on a
  mammalian reporter, in a record deposited at GenBank by the people who built
  the plasmid.

That last one is the reason this file exists. Rule 3.1 says the genetic code
table is explicit and never defaulted, because "a wrong table is a silently wrong
protein no assay catches for months". A synthetic fixture cannot test that rule
honestly, because the person writing the fixture also chooses the table. Here
somebody else chose it, and chose wrong, years before BT5 existed.

## A measured caveat on cost, and why the default-speed test is short

`design()` on this map is slow ON `main`, and not because the map is long. The
cost climbs superlinearly in PROTEIN length here and does not on the synthetic
map. Same settings (`sweep_steps=1`, 12-variant nulls), same protein, one
machine:

    protein   synthetic_mcs_ef1a (2000 bp)   real pFTMGW (8928 bp)
     39 aa                --                        1.0 s
     72 aa                --                        4.7 s
    100 aa               0.4 s                     13.4 s
    140 aa               0.5 s                    255.5 s

Shrinking the null does nothing: 140 aa costs 255 s at 12 draws and 256 s at 2.
So the null is not the term, and neither is the length of the construct being
scanned -- 39 aa on the same 8928 bp map is a second.

A profile at 100 aa says where it goes: 19.3 s of 22.8 s is inside
`solver/repair.py::repair`, which calls `solver/catalog.py::find` 758 times --
about 250 repair iterations per sweep pick -- and each of those re-runs
`f1_direct_repeats`, whose k-mer duplicate scan over the whole 9 kb construct is
9.7 s of the total on its own. This map contains real repeated DNA its authors
shipped (see `test_liabilities_in_the_vectors_own_dna_are_reported_not_silently_owned`),
none of it inside the insert, so repair keeps re-attacking breaches no codon
choice can clear. The synthetic fixture has no such repeats and converges in a
handful of iterations, which is why 100 aa -> 140 aa costs it a tenth of a
second.

THAT IS ISSUE #111, AND IT IS ALREADY FIXED. Measured on
`claude/m1-repair-stagnation` (PR #120, "stop re-attacking a breach while the
sequence stands still"), same map, same settings, only `solver/repair.py`
differing from `main`:

    protein     main     with #120
    100 aa     13.4 s       5.8 s
    140 aa    255.5 s       6.5 s

The superlinearity is gone, not merely reduced. So no new issue is filed for
this: the number in the first table is what `main` costs today, and #120 is the
fix. The first thing to do when #120 lands is delete the `slow` marker below and
give the default-speed test the full 140-residue fixture protein -- at 6.5 s it
no longer needs to be a separate slow test at all.

Until then the default-speed test designs a 72-residue protein, and
`TestTheFullFixtureProtein` runs the lane's real 140-residue one under `slow` --
which today means nothing runs it: `ci.yml` and `scripts/gates.sh` both pass
`-m "not slow"` and the repo has no nightly job. That is called out in the PR
rather than hidden here, because a test nobody runs covers nothing.

## What this file does and does not let the plan record

G2 -- annotation quality against a real map -- is what #74 asked for and what the
assertions below cover.

G7 ("500 residues, end to end, in 10 s") it does NOT cover, and nothing here
should be read as evidence for it. `test_timing.py` measures G7 against the
synthetic map alone, and on `main` this map takes four minutes at 140 residues.
With #120 the picture changes, but G7 at 500 residues on a real vector has not
been measured by anything in this repo.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import HostId, Modality
from bt5.core.registry import all_specs, discover
from bt5.core.types import Interval
from bt5.design import design
from bt5.vector import read_genbank
from bt5.vector.annotate import ORIGIN_QUALIFIER
from bt5.vector.backbone import VectorBackbone, VectorError

ROOT = Path(__file__).resolve().parents[4]
REAL_MAP = ROOT / "tests" / "data" / "backbones" / "real_lenti_pFTMGW_EF177827.gb"

CODE = FileTableProvider().genetic_code(1)

#: The MCS the vector annotates, as a half-open interval. Written out rather than
#: read from the feature table on purpose: `TestTheMapIsWhatThisFileThinksItIs`
#: checks these two against each other, so a silent re-download that moved the
#: site fails loudly instead of relocating every test in the file with it.
MCS = Interval(2766, 2835, 1)

#: The eGFP CDS's `/transl_table`, which is wrong for a mammalian reporter.
VECTOR_TABLE = 11

#: 72 residues. Short for the reason the module docstring measures; long enough to
#: fill more than three codons of every synonymous family the DP has to choose in.
PROTEIN_72 = "MKLVTAAFERSKSVQNYVVSTKDSPLYYLRKWVRSGYKFDCEEVGLREHQGPAATYTPTQAIWRLTLPSPLL"


def real_map() -> VectorBackbone:
    return read_genbank(REAL_MAP)


def run(protein: str, **kw: Any) -> Any:
    """`design()` into the annotated MCS, with the sweep and nulls at test size.

    The site is passed explicitly. `choose_site()` would land on the eGFP CDS
    instead -- see `TestTheWrongTableIsRefusedBothWays`, which is about exactly
    that -- and this map annotates no `promoter`, so `cloning_sites()` finds
    nothing to fall back to.
    """
    discover()
    params: dict[str, Any] = {
        "backbone": real_map(),
        "protein": protein,
        "table_id": 1,
        "modality": Modality.LENTIVIRAL,
        "hosts": [HostId.HEK293],
        "site_interval": MCS,
        "site_label": "MCS",
        "sweep_steps": 1,
        "null_sizes": {cls.id: 12 for cls in all_specs()},
    }
    params.update(kw)
    return design(**params)


@pytest.fixture(scope="module")
def designed() -> Any:
    """One `design()` run shared by the round-trip assertions.

    Module-scoped because it is the expensive thing in this file and every
    assertion below reads the same result; re-running it per test would multiply
    the cost in the docstring's table by the number of tests.
    """
    return run(PROTEIN_72)


class TestTheMapIsWhatThisFileThinksItIs:
    """Preconditions. Every other test in this file reads a coordinate or a
    qualifier off this record, so if the file on disk is not the record those
    tests were written against, they should say THAT rather than fail somewhere
    downstream with a coordinate nobody can place."""

    def test_it_is_the_deposited_record_unedited(self) -> None:
        bb = real_map()
        assert bb.length == 8928
        assert bb.is_circular
        assert bb.name == "EF177827"
        assert len(bb.features) == 12

    def test_the_mcs_constant_matches_the_annotation(self) -> None:
        bb = real_map()
        mcs = [f for f in bb.features if "MCS" in "".join(f.qualifiers.get("note", ()))]
        assert len(mcs) == 1, "the record annotates exactly one MCS"
        assert mcs[0].interval == MCS
        assert mcs[0].kind == "misc_feature", (
            "the MCS is a bare misc_feature, which is why cloning_sites() cannot see it"
        )

    def test_the_vector_annotates_the_bacterial_table_on_egfp(self) -> None:
        """The defect this file is really here for. Not ours -- GenBank's copy of
        somebody's real plasmid says table 11 on a mammalian reporter."""
        bb = real_map()
        cds = [f for f in bb.features if f.kind == "CDS"]
        assert len(cds) == 1
        assert cds[0].qualifiers["transl_table"] == (str(VECTOR_TABLE),)
        assert cds[0].qualifiers["gene"] == ("eGFP",)

    def test_the_record_annotates_no_promoter(self) -> None:
        """`cloning_sites()` keys off a `promoter` feature, and this map has none:
        its CMV promoter lives in prose inside the 5' LTR's `/note`. Real maps do
        this; the synthetic fixtures do not, which is why site selection here
        falls all the way through to the annotated-CDS branch."""
        bb = real_map()
        assert not [f for f in bb.features if f.kind.lower() == "promoter"]
        ltr = [f for f in bb.features if f.kind == "repeat_region"][0]
        assert "CMV promoter" in ltr.qualifiers["note"][0]


class TestTheWrongTableIsRefusedBothWays:
    """Rule 3.1, exercised by a record somebody else mis-annotated.

    Both directions refuse, and that is the whole point: there is no table BT5
    will accept here without the user resolving the contradiction themselves. A
    design that silently picked either one would produce a real construct with a
    real wrong protein.
    """

    def test_deferring_to_the_host_is_refused_because_the_vector_disagrees(self) -> None:
        with pytest.raises(VectorError, match=r"annotates /transl_table=11 at 'eGFP'"):
            design(
                backbone=real_map(),
                protein=PROTEIN_72,
                table_id=1,
                modality=Modality.LENTIVIRAL,
                hosts=[HostId.HEK293],
                sweep_steps=1,
                null_sizes={},
            )

    def test_deferring_to_the_vector_is_refused_because_the_host_disagrees(self) -> None:
        with pytest.raises(ValueError, match=r"host hek293 is locked to NCBI translation table 1"):
            design(
                backbone=real_map(),
                protein=PROTEIN_72,
                table_id=VECTOR_TABLE,
                modality=Modality.LENTIVIRAL,
                hosts=[HostId.HEK293],
                sweep_steps=1,
                null_sizes={},
            )

    def test_naming_the_site_explicitly_is_how_the_user_resolves_it(self) -> None:
        """The escape hatch, and it is deliberately not a flag that says 'trust
        the vector'. Pointing at the MCS instead of the CDS makes the disagreement
        moot, because the span being replaced is not the mis-annotated one."""
        res = run(PROTEIN_72[:40])
        assert res.result.candidates


class TestTheDesignItself:
    def test_every_candidate_encodes_the_protein(self, designed: Any) -> None:
        assert designed.result.candidates
        for candidate in designed.result.candidates:
            assert CODE.translate(candidate.cds)[:-1] == PROTEIN_72

    def test_the_insert_lands_in_the_annotated_mcs(self, designed: Any) -> None:
        assert designed.assembly.cds_interval.start == MCS.start

    def test_the_construct_is_not_rotated(self, designed: Any) -> None:
        """The coordinate assertions below compare construct positions against
        backbone positions directly, which is only legitimate at rotation 0. The
        assembler rotates when an insert would straddle the origin; this site is
        mid-plasmid, so it does not. Pinned rather than assumed, because a
        rotation would shift every one of those comparisons by a constant and
        several would still pass."""
        assert designed.assembly.rotation == 0


class TestTheGenBankRoundTrips:
    """G2. The output is read back with the same reader that read the input, and
    every fact the vector's authors recorded is still there afterwards."""

    def test_it_parses_and_the_length_is_the_constructs(self, designed: Any) -> None:
        # Through a handle: read_genbank path-tests a string first and OSErrors on
        # text this long.
        parsed = read_genbank(io.StringIO(designed.genbank))
        assert parsed.length == designed.assembly.construct.length

    def test_topology_and_locus_name_survive(self, designed: Any) -> None:
        parsed = read_genbank(io.StringIO(designed.genbank))
        assert parsed.is_circular, "a lentiviral transfer plasmid that came back linear"
        assert parsed.name == real_map().name

    def test_every_vector_feature_but_the_replaced_one_survives(self, designed: Any) -> None:
        """Feature-for-feature, not a count. A count passes while two features
        swap places, and on a map with overlapping `misc_feature`s that is the
        failure mode to expect."""
        parsed = read_genbank(io.StringIO(designed.genbank))
        kept = {
            (f.kind, tuple(sorted(f.qualifiers.get("note", ()) + f.qualifiers.get("gene", ()))))
            for f in parsed.features
            if f.qualifiers.get(ORIGIN_QUALIFIER) == ("provided",)
        }
        for feature in real_map().features:
            if feature.interval == MCS:
                continue  # the span the insert replaced; it is gone by design
            key = (
                feature.kind,
                tuple(
                    sorted(feature.qualifiers.get("note", ()) + feature.qualifiers.get("gene", ()))
                ),
            )
            assert key in kept, f"{feature.kind} {feature.qualifiers} did not survive"

    def test_the_replaced_mcs_is_the_only_thing_that_went_missing(self, designed: Any) -> None:
        parsed = read_genbank(io.StringIO(designed.genbank))
        provided = [
            f for f in parsed.features if f.qualifiers.get(ORIGIN_QUALIFIER) == ("provided",)
        ]
        assert len(provided) == len(real_map().features) - 1
        assert not [f for f in provided if "MCS" in "".join(f.qualifiers.get("note", ()))]

    def test_the_egfp_cds_keeps_every_qualifier_including_the_wrong_table(
        self, designed: Any
    ) -> None:
        """The richest feature on the record: `/translation`, `/protein_id`,
        `/product`, `/codon_start` and the mis-annotated `/transl_table`. BT5 must
        not quietly CORRECT the table it just refused to design against -- the
        vector said 11 and the output has to keep saying 11, or the annotation
        stops describing the plasmid the user actually holds."""
        parsed = read_genbank(io.StringIO(designed.genbank))
        out = [f for f in parsed.features if f.qualifiers.get("gene") == ("eGFP",)]
        src = [f for f in real_map().features if f.qualifiers.get("gene") == ("eGFP",)]
        assert len(out) == len(src) == 2  # the gene and the CDS
        out_cds = next(f for f in out if f.kind == "CDS")
        src_cds = next(f for f in src if f.kind == "CDS")
        for key, value in src_cds.qualifiers.items():
            assert out_cds.qualifiers[key] == value, f"eGFP lost or changed /{key}"

    def test_features_upstream_of_the_insert_do_not_move(self, designed: Any) -> None:
        parsed = read_genbank(io.StringIO(designed.genbank))
        by_note = {
            "".join(f.qualifiers.get("note", ())): f
            for f in parsed.features
            if f.qualifiers.get(ORIGIN_QUALIFIER) == ("provided",)
        }
        for feature in real_map().features:
            note = "".join(feature.qualifiers.get("note", ()))
            if not note or feature.interval.end > MCS.start or feature.interval == MCS:
                continue
            assert by_note[note].interval == feature.interval, f"{note} moved"

    def test_features_downstream_shift_by_exactly_the_length_delta(self, designed: Any) -> None:
        """The arithmetic a coordinate bug shows up in. Every feature after the
        site moves by the same amount, and that amount is the change in the
        construct's total length -- insert minus the MCS it replaced."""
        before = real_map()
        parsed = read_genbank(io.StringIO(designed.genbank))
        delta = parsed.length - before.length
        assert delta == designed.assembly.cds_interval.length - MCS.length
        by_note = {
            "".join(f.qualifiers.get("note", ())): f
            for f in parsed.features
            if f.qualifiers.get(ORIGIN_QUALIFIER) == ("provided",)
        }
        moved = 0
        for feature in before.features:
            note = "".join(feature.qualifiers.get("note", ()))
            if not note or feature.interval.start < MCS.end or feature.interval == MCS:
                continue
            out = by_note[note]
            assert out.interval.start == feature.interval.start + delta, f"{note} start"
            assert out.interval.end == feature.interval.end + delta, f"{note} end"
            moved += 1
        assert moved, "no downstream feature was checked; the map or MCS changed"


class TestProvenanceOnRealAnnotation:
    """Whose annotation is whose, on a record with 12 features somebody else
    wrote. The three classes have to stay separable after the merge, or a user
    cannot tell what BT5 added from what their vector already said."""

    def test_the_designed_insert_is_marked_designed_and_carries_the_asked_table(
        self, designed: Any
    ) -> None:
        parsed = read_genbank(io.StringIO(designed.genbank))
        made = [f for f in parsed.features if f.qualifiers.get(ORIGIN_QUALIFIER) == ("designed",)]
        assert len(made) == 1
        assert made[0].kind == "CDS"
        assert made[0].interval.start == MCS.start
        assert made[0].qualifiers["transl_table"] == ("1",), (
            "the insert is annotated with the table the DESIGN used, not the "
            "table the vector wrongly claims 150 bp downstream"
        )

    def test_no_provided_feature_is_marked_designed(self, designed: Any) -> None:
        parsed = read_genbank(io.StringIO(designed.genbank))
        for f in parsed.features:
            if f.qualifiers.get(ORIGIN_QUALIFIER) == ("designed",):
                assert f.qualifiers.get("gene") != ("eGFP",)
                assert not f.qualifiers.get("note", ())

    def test_liabilities_in_the_vectors_own_dna_are_reported_not_silently_owned(
        self, designed: Any
    ) -> None:
        """Real vectors contain real repeats. This one ships a 72 bp tandem repeat
        and two 21 bp ones, none of them anywhere near the insert -- so BT5 cannot
        fix them by codon choice and must not imply it did. They come back as
        `noted`, which is neither `provided` (the user did not annotate them) nor
        `designed` (BT5 did not create them).

        A synthetic fixture has none of this: it was written by someone who would
        not have put a tandem repeat in it.
        """
        parsed = read_genbank(io.StringIO(designed.genbank))
        noted = [f for f in parsed.features if f.qualifiers.get(ORIGIN_QUALIFIER) == ("noted",)]
        assert noted, "the real vector's repeated DNA went unreported"
        insert = designed.assembly.cds_interval
        for f in noted:
            assert f.qualifiers["note"], "a liability with no explanation is not a report"
            assert f.interval.end <= insert.start or f.interval.start >= insert.end, (
                f"liability at {f.interval} overlaps the insert; the solver owns that span"
            )

    def test_the_report_still_refuses_to_predict(self, designed: Any) -> None:
        """The one claim BT5 never makes, checked on the output a real map
        produces rather than only on a fixture's."""
        assert "never a predicted expression level" in designed.rendered


@pytest.mark.slow
class TestTheFullFixtureProtein:
    """#74 asked for the end-to-end run at full size, and this is it.

    Marked `slow` because on `main` it takes about four minutes -- see the module
    docstring's table. With PR #120 the same run is 6.5 s, at which point this
    marker should go and this protein should simply be the default-speed one.
    """

    def test_a_140_residue_design_round_trips_on_the_real_map(self) -> None:
        protein = PROTEIN_72 + (
            "NVDVWQNSCKSLQHTASWKKHRFGLFTLVISPLIRLGEVASLCGLCEHTATSEVKVCPIDCLQSPTSF"
        )
        assert len(protein) == 140
        res = run(protein)
        assert res.result.candidates
        for candidate in res.result.candidates:
            assert CODE.translate(candidate.cds)[:-1] == protein
        parsed = read_genbank(io.StringIO(res.genbank))
        assert parsed.length == res.assembly.construct.length
        assert parsed.is_circular
        provided = [
            f for f in parsed.features if f.qualifiers.get(ORIGIN_QUALIFIER) == ("provided",)
        ]
        assert len(provided) == len(real_map().features) - 1
