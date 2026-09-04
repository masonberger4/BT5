"""The vendor registry: one namespace, one default, every number sourced.

The load-bearing test here is `TestTheGCBandsMatchTheLadder`. `global_gc` is the
only field in a profile that no rule reads yet -- E2 carries the universal band
and #43 V3 wires the per-vendor one -- so it is the field most able to drift
wrong unnoticed. It is pinned directly to the 18-probe ladder recorded in
docs/design/vendor-gc-calibration.md, which is the measurement it exists to
encode. If a band moves without the ladder moving, this fails.
"""

from __future__ import annotations

import datetime as dt

import pytest
from bt5.rules.fragment import TWIST_ADAPTER_ON, TWIST_GENE_FRAGMENT
from bt5.rules.vendors import (
    DEFAULT_SELECTION,
    DEFAULT_VENDOR,
    PROFILES,
    VendorProfile,
    accepting_length,
    all_keys,
    orderable,
    orderable_keys,
    profile,
)

#: The global-GC half of the 18-probe ladder, verdicts as the vendors returned
#: them. `True` is "they will build it" -- IDT's "Accepted, Moderate Complexity"
#: is an acceptance with a surcharge, not a refusal, so 70% and 75% are True.
#: docs/design/vendor-gc-calibration.md
LADDER: dict[str, tuple[tuple[int, bool], ...]] = {
    "twist_gene_fragment": (
        (20, False),  # Not Accepted
        (25, False),  # Not Accepted
        (30, True),
        (40, True),
        (50, True),
        (60, True),
        (65, True),
        (70, True),
        (75, True),
        (80, True),  # Standard -- not "complex", Standard
    ),
    "idt_gblocks": (
        (20, False),  # Denied, score 64.6
        (25, False),  # Denied, score 27.7
        (30, True),
        (40, True),
        (50, True),
        (60, True),
        (65, True),
        (70, True),  # Accepted, Moderate (14.2)
        (75, True),  # Accepted, Moderate (21.2)
        (80, False),  # Denied, score 28.2
    ),
}


class TestTheGCBandsMatchTheLadder:
    @pytest.mark.parametrize("key", sorted(LADDER))
    def test_the_band_admits_exactly_what_the_vendor_built(self, key: str) -> None:
        band = PROFILES[key].global_gc
        assert band is not None
        lo, hi = band
        for pct, accepted in LADDER[key]:
            gc = pct / 100
            assert (lo <= gc <= hi) is accepted, (
                f"{key} at {pct}% GC: the ladder says accepted={accepted}, "
                f"the profile band {band} says otherwise"
            )

    def test_the_two_vendors_disagree_at_the_top_and_that_is_the_point(self) -> None:
        """80% GC is Standard at Twist and Denied at IDT.

        A single shipped GC ceiling cannot serve both, which is the whole reason
        the band is a per-vendor field rather than one constant. If these ever
        become equal, either a vendor changed or somebody collapsed the finding.
        """
        twist = PROFILES["twist_gene_fragment"].global_gc
        idt = PROFILES["idt_gblocks"].global_gc
        assert twist is not None
        assert idt is not None
        assert twist[1] > idt[1]
        assert twist[1] >= 0.80
        assert idt[1] < 0.80

    def test_the_floor_is_shared_because_the_measurement_was(self) -> None:
        """Both vendors refuse <=25% and accept >=30%; one measurement, one floor."""
        floors = {p.global_gc[0] for p in PROFILES.values() if p.global_gc is not None}
        assert len(floors) == 1, f"the ladder found one floor, the registry has {floors}"


class TestInheritanceStaysInsideAVendor:
    """A value may be carried between products of one vendor, never across two.

    The ladder's finding IS the vendor difference, so a number copied from IDT
    onto a Twist profile would be transferring away the thing that was measured.
    These two tests state that as a property of the whole registry rather than as
    a comment on the two profiles that currently inherit.
    """

    # `global_gc` is deliberately absent: its FLOOR is shared across vendors by
    # measurement, so the whole-band version of this property is the test below.
    @pytest.mark.parametrize("field", ["homopolymer_at", "homopolymer_gc"])
    def test_a_shared_value_implies_a_shared_vendor(self, field: str) -> None:
        by_value: dict[object, set[str]] = {}
        for p in PROFILES.values():
            value = getattr(p, field)
            if value is None:
                continue
            by_value.setdefault(value, set()).add(p.vendor)
        for value, vendors in by_value.items():
            assert len(vendors) == 1, f"{field}={value} appears under vendors {vendors}"

    def test_no_gc_band_is_shared_across_vendors(self) -> None:
        """The floor is shared by measurement; the BAND must not be."""
        by_band: dict[object, set[str]] = {}
        for p in PROFILES.values():
            if p.global_gc is None:
                continue
            by_band.setdefault(p.global_gc, set()).add(p.vendor)
        for band, vendors in by_band.items():
            assert len(vendors) == 1, f"band {band} appears under vendors {vendors}"

    def test_an_inheriting_profile_says_so(self) -> None:
        """eBlocks carries gBlocks' numbers; the notes have to admit it, because
        `last_verified` would otherwise date a measurement that never happened
        for that product."""
        assert "INHERITED" in PROFILES["idt_eblocks"].notes


class TestTheOneDefault:
    def test_the_default_is_real_and_orderable(self) -> None:
        assert DEFAULT_VENDOR in PROFILES
        assert PROFILES[DEFAULT_VENDOR].is_orderable

    def test_every_rule_that_takes_a_selection_gets_the_same_one(self) -> None:
        """The bug: E1 defaulted to IDT while E4-E7 and E9 defaulted to Twist, so
        with nothing chosen BT5 answered with IDT's run limits and Twist's
        lengths. One import, one answer -- and now that answer is a value object,
        so the property is stronger: every rule holds the SAME selection, not
        merely the same string."""
        from bt5.rules.catalog.e1_homopolymers import Homopolymers
        from bt5.rules.catalog.e2_gc_band import GCBand
        from bt5.rules.catalog.e4_gc_extent import GCExtent
        from bt5.rules.catalog.e5_synthesis_repeats import SynthesisRepeats
        from bt5.rules.catalog.e6_repeat_density import RepeatDensity
        from bt5.rules.catalog.e7_short_tandem_repeats import ShortTandemRepeats
        from bt5.rules.catalog.e9_length_tiers import LengthTiers

        for rule in (
            Homopolymers(),
            GCBand(),
            GCExtent(),
            SynthesisRepeats(),
            RepeatDensity(),
            ShortTandemRepeats(),
            LengthTiers(),
        ):
            assert rule.vendors == DEFAULT_SELECTION, type(rule).__name__
        assert DEFAULT_SELECTION.keys == (DEFAULT_VENDOR,)

    def test_the_schema_enum_offers_only_configurations_that_exist(self) -> None:
        """The `vendors` array's enum, and its default, must name real products.

        Guarded by the key name: renaming `vendor` -> `vendors` had to move this
        check with it, or it would pass vacuously on the `if "vendors" not in
        props: continue` line and stop protecting anything -- the exact way a
        schema check quietly dies.
        """
        from bt5.core.registry import all_specs, discover

        discover()
        checked = 0
        for spec in all_specs():
            props = spec.param_schema.get("properties", {})
            assert isinstance(props, dict)
            if "vendors" not in props:
                continue
            checked += 1
            enum = props["vendors"]["items"]["enum"]
            for key in enum:
                assert key in PROFILES, f"{spec.id} offers {key!r}, which is not a profile"
            for key in props["vendors"]["default"]:
                assert key in enum, f"{spec.id} defaults to {key!r}, not in its own enum"
        assert checked == 7, f"expected seven vendor-taking rules, checked {checked}"


class TestLookup:
    def test_an_unknown_key_lists_the_real_ones(self) -> None:
        with pytest.raises(ValueError, match="unknown vendor 'acme'"):
            profile("acme")

    def test_none_resolves_but_is_not_orderable(self) -> None:
        assert profile("none").is_orderable is False
        with pytest.raises(ValueError, match="not orderable from anyone"):
            orderable("none")

    def test_orderable_keys_excludes_none_and_all_keys_does_not(self) -> None:
        assert "none" not in orderable_keys()
        assert "none" in all_keys()
        assert set(orderable_keys()) | {"none"} == set(all_keys())

    def test_accepting_length_finds_the_configurations_that_fit(self) -> None:
        # 200 bp: only gBlocks reach below 300.
        assert accepting_length(200) == ("idt_gblocks",)
        # 4000 bp: over the gBlocks and eBlocks ceilings, inside Twist's.
        assert accepting_length(4000) == (
            "twist_gene_fragment",
            "twist_gene_fragment_adapter_on",
        )
        assert accepting_length(50) == ()

    def test_accepting_length_can_exclude_the_configuration_being_reported_on(self) -> None:
        assert "idt_gblocks" not in accepting_length(1000, exclude="idt_gblocks")


class TestTheProfileRefusesToBeHalfTrue:
    def test_a_length_without_the_rest_is_refused(self) -> None:
        with pytest.raises(ValueError, match="half-specified"):
            VendorProfile(
                key="twist_gene_fragment",
                vendor="Twist Bioscience",
                product="Gene Fragment",
                adapters=TWIST_GENE_FRAGMENT,
                homopolymer_at=13,
            )

    def test_adapters_must_name_their_own_configuration(self) -> None:
        with pytest.raises(ValueError, match="one configuration, one name"):
            VendorProfile(
                key="twist_gene_fragment",
                vendor="Twist Bioscience",
                product="Gene Fragment",
                adapters=TWIST_ADAPTER_ON,
            )

    def test_an_orderable_profile_must_carry_its_provenance(self) -> None:
        """`notes` and `last_verified` are not decoration. Vendor numbers drift --
        Twist's own homopolymer limit moved from 14 to 30 bp between 2023 and
        2026 -- and an undated number is one that goes wrong silently."""
        common = {
            "key": "twist_gene_fragment",
            "vendor": "Twist Bioscience",
            "product": "Gene Fragment",
            "adapters": TWIST_GENE_FRAGMENT,
            "length_bp": (300, 5000),
            "homopolymer_at": 13,
            "homopolymer_gc": 13,
            "global_gc": (0.28, 0.80),
        }
        with pytest.raises(ValueError, match="measured, published or inherited"):
            VendorProfile(**common, last_verified="2026-08-28")  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="need a verification date"):
            VendorProfile(**common, notes="measured")  # type: ignore[arg-type]

    def test_an_inverted_range_is_refused(self) -> None:
        with pytest.raises(ValueError, match="length range .* is empty"):
            VendorProfile(
                key="twist_gene_fragment",
                vendor="Twist Bioscience",
                product="Gene Fragment",
                adapters=TWIST_GENE_FRAGMENT,
                length_bp=(5000, 300),
                homopolymer_at=13,
                homopolymer_gc=13,
                global_gc=(0.28, 0.80),
                last_verified="2026-08-28",
                notes="measured",
            )


class TestEveryProfileIsWellFormed:
    @pytest.mark.parametrize("key", sorted(PROFILES))
    def test_the_key_matches_the_profile(self, key: str) -> None:
        assert PROFILES[key].key == key

    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_last_verified_is_a_real_date(self, key: str) -> None:
        dt.date.fromisoformat(PROFILES[key].last_verified)

    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_it_names_a_vendor_and_a_product_a_person_can_order(self, key: str) -> None:
        """A finding says "orderable instead as X", so X has to be a product and
        not a tier. `twist_standard`, the old homopolymer key, was a complexity
        VERDICT the checker returns -- there is nothing to add to a cart."""
        assert PROFILES[key].vendor.strip()
        assert PROFILES[key].product.strip()

    def test_only_the_adapter_on_option_carries_adapters(self) -> None:
        with_adapters = {k for k in all_keys() if PROFILES[k].adapters.total}
        assert with_adapters == {"twist_gene_fragment_adapter_on"}


class TestTheSidecarAndTheCodeAgree:
    """`_provenance.json` is only worth having if it cannot drift from the code.

    A provenance file that quietly disagrees with the values it documents is
    worse than none: it reads as a source for a number that came from somewhere
    else. So every number in a profile is checked against the sidecar entry that
    claims to source it, and every profile has to have one.
    """

    @staticmethod
    def sidecar() -> dict:
        import json
        from pathlib import Path

        import bt5.rules.vendors as mod

        path = Path(mod.__file__).with_name("_provenance.json")
        return json.loads(path.read_text(encoding="utf-8"))

    def test_every_orderable_profile_is_documented(self) -> None:
        doc = self.sidecar()
        documented = {k for k in doc if not k.startswith("_")}
        assert documented == set(orderable_keys())

    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_the_numbers_match(self, key: str) -> None:
        entry = self.sidecar()[key]
        p = PROFILES[key]
        assert entry["last_verified"] == p.last_verified
        assert tuple(entry["length_bp"]["value"]) == p.length_bp
        assert tuple(entry["homopolymer"]["value"]) == (p.homopolymer_at, p.homopolymer_gc)
        assert tuple(entry["global_gc"]["value"]) == p.global_gc
        assert entry["vendor"] == p.vendor
        assert entry["product"] == p.product

    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_every_value_declares_how_it_was_established(self, key: str) -> None:
        entry = self.sidecar()[key]
        for field in ("adapters", "length_bp", "homopolymer", "global_gc"):
            status = entry[field]["status"]
            assert status in {"MEASURED", "PUBLISHED", "INHERITED"}, f"{key}.{field}={status}"
            if status == "INHERITED":
                # Inside a vendor only. A number carried across vendors would be
                # discarding the one difference the ladder measured.
                parent = entry[field]["from"]
                assert PROFILES[parent].vendor == PROFILES[key].vendor, (
                    f"{key}.{field} inherits from {parent}, a different vendor"
                )
            else:
                assert entry[field]["source"].strip()

    def test_a_measured_value_says_when_it_was_measured(self) -> None:
        doc = self.sidecar()
        for key in orderable_keys():
            for entry in doc[key].values():
                if isinstance(entry, dict) and entry.get("status") == "MEASURED":
                    dt.date.fromisoformat(entry["measured_on"])

    def test_the_dropped_keys_say_why_they_went(self) -> None:
        """`genscript` and `twist_standard` were removed, not renamed. A registry
        that silently loses a key is a config that silently changes meaning."""
        dropped = self.sidecar()["_dropped"]
        assert set(dropped) == {"genscript", "twist_standard"}
        for key, entry in dropped.items():
            assert entry["why"].strip(), key
            assert key not in PROFILES
