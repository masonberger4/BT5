"""The three simultaneous contexts, and the genetic-code guard."""

from __future__ import annotations

import pytest
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)


def slot(
    role: str = "propagation", host: HostId = HostId.E_COLI_K12, table: int = 11
) -> ContextSlot:
    return ContextSlot(role, host, Modality.PLASMID_TRANSIENT, table)  # type: ignore[arg-type]


class TestGeneticCodeGuard:
    def test_accepts_the_locked_table(self) -> None:
        assert slot().table_id == 11

    def test_rejects_a_mismatched_table(self) -> None:
        """A wrong table is a silently wrong protein no assay catches for months."""
        with pytest.raises(ValueError, match="locked to NCBI translation table"):
            ContextSlot("target", HostId.HEK293, Modality.LENTIVIRAL, 11)

    def test_mammalian_hosts_are_locked_to_table_1(self) -> None:
        ContextSlot("target", HostId.HEK293, Modality.LENTIVIRAL, 1)
        ContextSlot("target", HostId.CHO, Modality.PLASMID_STABLE, 1)


class TestDesignContext:
    def ctx(self, *slots: ContextSlot) -> DesignContext:
        return DesignContext(
            slots=slots or (slot(),),
            cassette_orientation=1,
            seed=42,
            screen=BiosecurityVerdict("not_run"),
        )

    def test_requires_at_least_one_slot(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            DesignContext(
                slots=(), cassette_orientation=1, seed=1, screen=BiosecurityVerdict("not_run")
            )

    def test_rejects_more_than_three_slots(self) -> None:
        slots = (
            slot("propagation"),
            ContextSlot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1),
            ContextSlot("target", HostId.HUMAN, Modality.LENTIVIRAL, 1),
        )
        extra = ContextSlot("target", HostId.CHO, Modality.LENTIVIRAL, 1)
        with pytest.raises(ValueError, match="at most three"):
            DesignContext(
                slots=(*slots, extra),
                cassette_orientation=1,
                seed=1,
                screen=BiosecurityVerdict("not_run"),
            )

    def test_rejects_duplicate_roles(self) -> None:
        with pytest.raises(ValueError, match="duplicate"):
            self.ctx(slot("propagation"), slot("propagation", HostId.E_COLI_BL21))

    def test_the_compound_case_is_three_simultaneous_slots(self) -> None:
        """plasmid -> virus -> transduce: E. coli propagation, the packaging cell,
        and the target cell all constrain ONE sequence at the same time."""
        ctx = self.ctx(
            slot("propagation"),
            ContextSlot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1),
            ContextSlot("target", HostId.HUMAN, Modality.LENTIVIRAL, 1),
        )
        assert len(ctx.active_slots) == 3
        assert ctx.slot("producer") is not None
        assert ctx.slot("target").modality is Modality.LENTIVIRAL

    def test_disabled_slots_drop_out_of_active(self) -> None:
        disabled = ContextSlot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1, enabled=False)
        ctx = self.ctx(slot("propagation"), disabled)
        assert len(ctx.active_slots) == 1


class TestBiosecurityVerdict:
    def test_not_run_is_never_reported_as_clear(self) -> None:
        v = BiosecurityVerdict("not_run")
        assert v.status == "not_run"
        assert v.status != "clear", "a screen that did not run must not read as clean"

    def test_block_stops_the_pipeline(self) -> None:
        assert not BiosecurityVerdict("block").may_proceed
        assert BiosecurityVerdict("flag").may_proceed
        assert BiosecurityVerdict("clear").may_proceed
