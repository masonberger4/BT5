"""Issue #56: the error-free length, and the one that is deliberately absent.

`score/report.py` already behaves correctly -- `screening_burden` returns None for
`idt_gblocks` and `build_report` turns that into a stated degradation. What #56 names is
that the *provenance* was missing: `_provenance.json` documented every other vendor
number's source and said nothing at all about fidelity, so the absence lived only in a
code comment in another lane.

These tests bind the two together. `ERROR_FREE_BP` lives in `bt5/score/`, which is M3's
lane and read-only here; the sidecar is M4's. This file is the seam, and it is written so
that a change on either side that does not update the other fails loudly rather than
leaving the report authoritative about a number nobody sourced.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from bt5.rules.vendors import PROFILES, orderable_keys
from bt5.score.report import ERROR_FREE_BP

FIELD = "error_free_bp"


def sidecar() -> dict:
    import bt5.rules.vendors as mod

    return json.loads(Path(mod.__file__).with_name("_provenance.json").read_text(encoding="utf-8"))


class TestEveryProfileStatesItsFidelityProvenance:
    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_the_field_exists(self, key: str) -> None:
        """Including the absent one. A missing key is indistinguishable from one
        nobody considered, which is the whole complaint in #56."""
        assert FIELD in sidecar()[key], f"{key} does not say where its fidelity came from"

    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_the_status_is_one_of_the_four(self, key: str) -> None:
        assert sidecar()[key][FIELD]["status"] in {
            "MEASURED",
            "PUBLISHED",
            "INHERITED",
            "ABSENT",
        }

    @pytest.mark.parametrize("key", sorted(orderable_keys()))
    def test_a_sourced_value_carries_its_source(self, key: str) -> None:
        entry = sidecar()[key][FIELD]
        if entry["status"] in {"MEASURED", "PUBLISHED"}:
            assert entry["source"].strip()
        elif entry["status"] == "INHERITED":
            assert entry["from"] in PROFILES


class TestTheSidecarAndTheCodeAgree:
    """The seam. `ERROR_FREE_BP` is M3's and the sidecar is M4's, so nothing but a
    test stops them drifting into a report that cites a source for a different number.
    """

    def test_every_value_in_the_code_is_sourced_in_the_sidecar(self) -> None:
        doc = sidecar()
        for key, value in ERROR_FREE_BP.items():
            entry = doc[key][FIELD]
            assert entry["status"] != "ABSENT", (
                f"{key} has a value in ERROR_FREE_BP but is recorded ABSENT"
            )
            assert entry["value"] == value, (
                f"{key}: sidecar says {entry['value']}, code says {value}"
            )

    def test_every_absent_entry_is_absent_from_the_code_too(self) -> None:
        doc = sidecar()
        for key in orderable_keys():
            if doc[key][FIELD]["status"] == "ABSENT":
                assert key not in ERROR_FREE_BP, (
                    f"{key} is recorded ABSENT but ERROR_FREE_BP has a value for it. "
                    "If a figure was found, update the sidecar in the same change."
                )


class TestTheGblocksFigure:
    """The hole #56 named, now filled -- and pinned so it cannot be filled WRONGLY.

    This class replaces `TestTheGblocksGap`, which pinned the absence "so it cannot
    be filled quietly". It was filled loudly instead: a published IDT figure, its
    URL, its retrieval date and its verbatim wording, all in the sidecar. What
    those tests were really defending was never the emptiness -- it was that
    gBlocks must not silently wear eBlocks' number. Every one of those assertions
    survives below, aimed at the figure rather than at its absence.
    """

    def test_gblocks_is_recorded_published_with_its_source(self) -> None:
        entry = sidecar()["idt_gblocks"][FIELD]
        assert entry["status"] == "PUBLISHED"
        assert entry["value"] == 5000
        assert entry["source"].strip()
        assert entry["retrieved"].strip(), "an undated vendor figure is unmaintainable"

    def test_it_cites_the_gblocks_page_and_not_the_eblocks_one(self) -> None:
        """THE assertion this class exists for, and the one #56 called the subtle
        way to get it wrong. The two products publish the same 1:5000, so the
        value alone cannot tell a correct entry from an inherited one -- only the
        citation can. eBlocks and gBlocks are different manufacturing processes at
        one company, and a fidelity figure is a property of the process."""
        entry = sidecar()["idt_gblocks"][FIELD]
        assert "gblocks-gene-fragments" in entry["source"]
        assert "eblocks" not in entry["source"].lower()
        assert entry["status"] != "INHERITED"

    def test_it_is_the_standard_product_not_hifi(self) -> None:
        """gBlocks HiFi is a different product on the same page with a different
        figure (1:12,000). Applying it here would understate the screening burden
        by more than a factor of two -- the unsafe direction."""
        entry = sidecar()["idt_gblocks"][FIELD]
        assert entry["value"] != 12000
        assert "hifi" in entry["notes"].lower(), "the confusable product must be named"
        assert PROFILES["idt_gblocks"].length_bp == (125, 3000), (
            "the standard product's range; HiFi is 1000-3000 and would need its own key"
        )

    def test_the_default_vendor_now_carries_a_figure(self) -> None:
        """Why closing this mattered: with no figure, BT5's OWN default produced a
        report with no 'pick N colonies' line."""
        from bt5.rules.vendors import DEFAULT_VENDOR

        assert DEFAULT_VENDOR == "idt_gblocks"
        assert ERROR_FREE_BP[DEFAULT_VENDOR] == 5000

    def test_no_fidelity_figure_is_inherited_across_products(self) -> None:
        """Inheritance elsewhere in this file is within a VENDOR. For a fidelity
        figure it must also be within a PRODUCT, because the number is a property of
        the synthesis process rather than of the company."""
        doc = sidecar()
        for key in orderable_keys():
            entry = doc[key][FIELD]
            if entry["status"] != "INHERITED":
                continue
            parent = entry["from"]
            assert PROFILES[parent].vendor == PROFILES[key].vendor
            assert PROFILES[parent].product.split(",")[0] == PROFILES[key].product.split(",")[0], (
                f"{key} inherits a fidelity figure from {parent}, a different product"
            )


class TestThePolicyDocumentsTheNewStatus:
    def test_absent_is_named_in_the_about_block(self) -> None:
        """A status the policy does not name reads as a typo to the next reader."""
        assert "ABSENT" in sidecar()["_about"]["policy"]
