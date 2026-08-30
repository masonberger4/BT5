"""VendorSelection: the value a rule takes instead of a vendor string (#43 V2/V3).

The load-bearing test is `TestNothingReintroducesTheNLoop`. Issue #43 argued for
evaluating each rule once per selected vendor; with adapters off the table that
would produce N identical findings and N-fold a SOFT weight keyed by spec_id, so
the design evaluates ONCE and attributes per vendor instead. The pin proves the
four synthesis rules return identical findings under a three-vendor selection as
under one -- the moment they stop, someone has put the loop back.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.services import Services
from bt5.rules.vendors import (
    DEFAULT_SELECTION,
    DEFAULT_VENDOR,
    VendorSelection,
    accepting_length,
    orderable_keys,
    require_selection,
)
from conftest import construct, context


class TestOfRefusesEveryNonDesign:
    def test_empty_is_not_a_choice(self) -> None:
        with pytest.raises(ValueError, match="empty selection is not a choice"):
            VendorSelection.of()

    def test_an_unknown_key_reuses_the_registry_message(self) -> None:
        with pytest.raises(ValueError, match="unknown vendor 'acme'"):
            VendorSelection.of("acme")

    def test_none_cannot_be_combined_with_a_real_vendor(self) -> None:
        with pytest.raises(ValueError, match="not a vendor you can add another to"):
            VendorSelection.of("none", "idt_gblocks")

    def test_adapter_on_may_not_be_mixed_with_adapter_free(self) -> None:
        """The refusal that makes the different-molecule hazard unreachable.

        Compared by PAYLOAD, not object identity: two adapter-free products are a
        legal selection though their `Adapters` objects differ.
        """
        with pytest.raises(ValueError, match="physically different molecules"):
            VendorSelection.of("idt_gblocks", "twist_gene_fragment_adapter_on")

    def test_two_adapter_free_products_are_a_legal_selection(self) -> None:
        s = VendorSelection.of("idt_gblocks", "twist_gene_fragment")
        assert s.keys == ("idt_gblocks", "twist_gene_fragment")

    def test_duplicates_collapse_and_order_is_kept(self) -> None:
        assert VendorSelection.of("idt_gblocks", "idt_gblocks").keys == ("idt_gblocks",)
        assert VendorSelection.of("twist_gene_fragment", "idt_gblocks").keys == (
            "twist_gene_fragment",
            "idt_gblocks",
        )


class TestRequireSelection:
    def test_a_bare_string_is_refused_not_iterated(self) -> None:
        """The old call shape. Silently iterating it would spec a design for five
        single-character 'vendors'; TypeError makes the mistake loud."""
        with pytest.raises(TypeError, match="not a bare string"):
            require_selection("idt_gblocks")

    def test_a_selection_passes_through(self) -> None:
        s = VendorSelection.of("idt_gblocks")
        assert require_selection(s) is s


class TestLabelAndAdapters:
    def test_a_single_key_reads_exactly_as_before(self) -> None:
        s = VendorSelection.of("idt_gblocks")
        assert s.label == "idt_gblocks"
        # Byte-identical to IDT_GBLOCKS: empty payload, vendor is the key.
        assert (s.adapters.five, s.adapters.three, s.adapters.vendor) == ("", "", "idt_gblocks")

    def test_the_label_names_the_whole_selection(self) -> None:
        s = VendorSelection.of("idt_gblocks", "twist_gene_fragment")
        assert s.label == "idt_gblocks, twist_gene_fragment"
        assert s.adapters.vendor == s.label

    def test_default_selection_is_the_one_default(self) -> None:
        assert DEFAULT_SELECTION.keys == (DEFAULT_VENDOR,)


class TestOrderableOnly:
    def test_none_is_refused(self) -> None:
        with pytest.raises(ValueError, match="not orderable from anyone"):
            VendorSelection.of("none").orderable_only()

    def test_an_orderable_selection_returns_itself(self) -> None:
        s = VendorSelection.of("idt_gblocks")
        assert s.orderable_only() is s


class TestHomopolymerLimits:
    def test_the_stricter_axis_binds_and_names_its_publisher(self) -> None:
        s = VendorSelection.of("idt_gblocks", "twist_gene_fragment")
        (at_min, at_keys), (gc_min, gc_keys) = s.homopolymer_limits()
        # IDT (9, 5) is stricter on both axes than Twist (13, 13).
        assert (at_min, at_keys) == (9, ("idt_gblocks",))
        assert (gc_min, gc_keys) == (5, ("idt_gblocks",))

    def test_a_run_over_the_strict_limit_may_be_within_a_looser_vendors(self) -> None:
        s = VendorSelection.of("idt_gblocks", "twist_gene_fragment")
        # 11 A: over IDT's 9, under Twist's 13.
        assert s.homopolymer_accepts(11, "A/T") == ("twist_gene_fragment",)
        # 14 A: over both.
        assert s.homopolymer_accepts(14, "A/T") == ()


class TestLengthVerdicts:
    def test_accept_ambiguous_refuse(self) -> None:
        s = VendorSelection.of("idt_gblocks", "twist_gene_fragment")
        # 4002 bp: over gBlocks' 3000, inside Twist's 5000.
        assert dict(s.verdicts_for_length(4002)) == {
            "idt_gblocks": "refuse",
            "twist_gene_fragment": "accept",
        }
        # 200 bp: inside gBlocks (125-3000), under Twist's 300 floor.
        assert dict(s.verdicts_for_length(200)) == {
            "idt_gblocks": "accept",
            "twist_gene_fragment": "refuse",
        }

    def test_ambiguous_only_arises_for_an_adapter_on_order(self) -> None:
        # 280 bp insert + 44 bp adapter straddles Twist's 300 floor.
        s = VendorSelection.of("twist_gene_fragment_adapter_on")
        assert s.verdicts_for_length(280) == (("twist_gene_fragment_adapter_on", "ambiguous"),)


class TestAlternativesCannotDriftFromAcceptingLength:
    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    @pytest.mark.parametrize("bp", [100, 200, 299, 300, 1500, 3000, 4002, 6000])
    def test_single_key_alternatives_equal_accepting_length_exclude(
        self, key: str, bp: int
    ) -> None:
        """`alternatives_for` delegates to `accepting_length`, so the two answers
        to 'where else could this go' cannot diverge."""
        assert VendorSelection.of(key).alternatives_for(bp) == accepting_length(bp, exclude=key)

    def test_a_selected_vendor_is_never_its_own_alternative(self) -> None:
        s = VendorSelection.of("idt_gblocks", "twist_gene_fragment")
        for bp in (100, 200, 4002):
            assert not (set(s.alternatives_for(bp)) & set(s.keys))


def _dna(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


class TestNothingReintroducesTheNLoop:
    """The four synthesis rules read only adapters from a profile, and the owner
    orders none. So a three-vendor selection must return the SAME findings as one
    -- identical apart from `detail["vendor"]`, which is exactly what V3's second
    clause changed. A bare byte-identity claim would be false for E4 and E5, which
    already carried the vendor key; hence the modulo.
    """

    @pytest.fixture
    def svc(self) -> Services:
        from bt5.vector.kmers import ConstructKmerIndex

        return Services(
            fold=None,
            kmer=ConstructKmerIndex,
            tables=None,  # type: ignore[arg-type]
            rng=np.random.default_rng(42),
        )

    def _rules(self, sel: VendorSelection) -> list[object]:
        from bt5.rules.catalog.e4_gc_extent import GCExtent
        from bt5.rules.catalog.e5_synthesis_repeats import SynthesisRepeats
        from bt5.rules.catalog.e6_repeat_density import RepeatDensity
        from bt5.rules.catalog.e7_short_tandem_repeats import ShortTandemRepeats

        return [
            GCExtent(vendors=sel),
            SynthesisRepeats(vendors=sel),
            RepeatDensity(vendors=sel),
            ShortTandemRepeats(vendors=sel),
        ]

    def test_three_adapter_free_vendors_evaluate_like_one(self, svc: Services) -> None:
        one = VendorSelection.of("idt_gblocks")
        many = VendorSelection.of("idt_gblocks", "idt_eblocks", "twist_gene_fragment")
        # A construct rigged to trip every one of the four: a GC-extreme, repetitive insert.
        cds = "ATG" + "GGCGGCGGCGGCAGC" * 40 + _dna(300, 3) + "TAA"
        c = construct(cds, _dna(400, 9))
        for r_one, r_many in zip(self._rules(one), self._rules(many), strict=True):
            e1 = r_one.evaluate(c, context(), svc)  # type: ignore[attr-defined]
            e2 = r_many.evaluate(c, context(), svc)  # type: ignore[attr-defined]
            name = type(r_one).__name__
            assert e1.raw_score == e2.raw_score, name
            assert e1.n_evaluated == e2.n_evaluated, name
            assert e1.windows == e2.windows, name
            assert len(e1.breaches) == len(e2.breaches), name
            for b1, b2 in zip(e1.breaches, e2.breaches, strict=True):
                assert b1.interval == b2.interval, name
                assert b1.magnitude == b2.magnitude, name
                assert b1.fixable_by_codon_choice == b2.fixable_by_codon_choice, name
                # Everything but the vendor key is identical; the vendor key is
                # the whole point of V3's second clause.
                d1 = {k: v for k, v in b1.detail.items() if k != "vendor"}
                d2 = {k: v for k, v in b2.detail.items() if k != "vendor"}
                assert d1 == d2, name
                assert b1.detail["vendor"] == "idt_gblocks"
                assert b2.detail["vendor"] == "idt_gblocks, idt_eblocks, twist_gene_fragment"
