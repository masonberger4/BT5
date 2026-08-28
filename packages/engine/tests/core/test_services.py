"""The fold contract: what an energy has to carry, and what an absent engine is.

Every literature threshold in BT5 is a kcal/mol number calibrated against one
engine and one parameter set, so an energy that travels without its ruler is not
comparable to the threshold it is being checked against. These tests are about
what has to ride along with the number.
"""

from __future__ import annotations

from bt5.core.services import FoldEnergy


def energy(**kw: object) -> FoldEnergy:
    base: dict[str, object] = {
        "dg_kcal_mol": -10.9,
        "engine": "viennarna",
        "engine_version": "2.7.2",
        "param_set": "rna_turner2004",
    }
    base.update(kw)
    return FoldEnergy(**base)  # type: ignore[arg-type]


class TestFoldEnergyProvenance:
    def test_the_calibration_key_identifies_engine_and_parameters(self) -> None:
        """`Spec.engine_calibration` is declared in this form, so a rule can
        compare its calibration against the engine that produced a result."""
        assert energy().calibration_key == "viennarna:rna_turner2004"

    def test_conditions_default_to_viennarnas_own_defaults(self) -> None:
        """37 C and dangles=2 are what ViennaRNA uses; a mismatch here would put
        every default-constructed energy on a different ruler from the engine."""
        e = energy()
        assert e.temperature_c == 37.0
        assert e.dangles == 2


class TestStructureTravelsWithTheEnergy:
    def test_structure_is_empty_when_not_recorded(self) -> None:
        assert energy().structure == ""

    def test_a_recorded_structure_survives(self) -> None:
        """Rules stated in terms of an individual hairpin's position -- the
        cap-proximal ladder, the AAV hairpin flag -- cannot be answered by a
        scalar energy, and re-folding to recover the structure would be a second
        O(n^3) call for information the first one already had."""
        assert energy(structure="..((...))..").structure == "..((...)).."


class TestDuplexIsDistinguishable:
    def test_an_ordinary_fold_is_not_a_duplex(self) -> None:
        e = energy(structure="..((...))..")
        assert not e.is_duplex
        assert e.duplex_split is None

    def test_a_duplex_records_where_the_second_molecule_starts(self) -> None:
        """ViennaRNA's dimer structure carries no separator, so without the split
        index the two halves are indistinguishable in one flat string."""
        e = energy(structure="......(((((((....)))))))..", duplex_split=36)
        assert e.is_duplex
        assert e.duplex_split == 36

    def test_a_duplex_split_of_zero_still_reads_as_a_duplex(self) -> None:
        """Guards the obvious falsy-integer bug: `if duplex_split:` would call
        this intramolecular."""
        assert energy(duplex_split=0).is_duplex
