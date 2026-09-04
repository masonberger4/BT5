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


class TestTheGblocksGap:
    """The specific hole #56 names, pinned so it cannot be filled quietly."""

    def test_gblocks_is_recorded_absent(self) -> None:
        entry = sidecar()["idt_gblocks"][FIELD]
        assert entry["status"] == "ABSENT"
        assert entry["value"] is None

    def test_it_names_the_issue_and_what_it_blocks(self) -> None:
        entry = sidecar()["idt_gblocks"][FIELD]
        assert "56" in entry["issue"]
        assert entry["blocks"].strip()
        assert entry["to_close"].strip(), "an absence with no route to closing it is a shrug"

    def test_it_records_the_rejected_ways_of_filling_it(self) -> None:
        """The subtle one is inheriting eBlocks' number across IDT's two product
        lines. A fidelity figure comes out of a manufacturing process, not a design
        guideline, and the two products are different processes."""
        rejected = sidecar()["idt_gblocks"][FIELD]["rejected"]
        assert "eblocks" in rejected.lower()
        assert "5000" in rejected

    def test_the_gap_is_on_the_default_vendor_which_is_why_it_matters(self) -> None:
        from bt5.rules.vendors import DEFAULT_VENDOR

        assert DEFAULT_VENDOR == "idt_gblocks"
        assert DEFAULT_VENDOR not in ERROR_FREE_BP

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
