"""The default weight vectors, and why each number is what it is.

docs/PLAN.md is blunt about this: "The default weight vector is the product.
90% of users never move a slider." A preset is therefore a scientific artefact,
not a convenience -- it is the claim BT5 makes about what matters, for most
users the ONLY claim it makes, and it has to be defensible line by line. So
`rationale` is required and non-empty, and so is the `note` on any entry that
departs from the rule's own `default_weight`.

Presets key on `brief_ref`, not on spec id. That is deliberate. `brief_ref` is
the row in docs/research/brief.md the rule implements ("2.B1"), fixed by the
brief before any rule exists; the spec id contains a slug the rule's author
picks later. Keying on ids would mean this file guessing forty slugs that M4 has
not chosen yet, and going quietly stale on every one it guessed wrong. Keying on
brief_ref means a preset written today resolves correctly against a rule written
next month, and a preset naming a row nobody has implemented reports that as an
unimplemented objective instead of silently weighting nothing.

The weights are percentile weights. Everything in the sum has been normalised
against the same null, so they are commensurable by construction -- which is the
only reason a number like 1.0 next to a number like 0.2 means anything.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import get_args

from bt5.core.context import LOCKED_TRANSLATION_TABLE, ContextSlot, Modality, SlotRole
from bt5.core.registry import all_specs, discover
from bt5.core.spec import Enforcement, Spec


class PresetError(ValueError):
    """A preset is internally inconsistent, or asks for something incoherent."""


@dataclass(frozen=True, slots=True)
class WeightEntry:
    """One objective's weight in one preset.

    `note` is required whenever this preset weights an objective differently
    from the rule's own `default_weight`. The rule's `weight_provenance` says
    why the rule's default is what it is; it cannot also explain why THIS
    context departs from it, and a preset that silently overrides a documented
    default is exactly the unreviewable magic number the contract test exists
    to prevent.
    """

    brief_ref: str
    weight: float
    note: str = ""

    def __post_init__(self) -> None:
        if not self.brief_ref.strip():
            raise PresetError("a weight entry must name a brief_ref")
        if self.weight < 0.0:
            raise PresetError(
                f"{self.brief_ref}: weight must be >= 0, got {self.weight}. "
                f"A negative weight inverts an objective's direction, which is a "
                f"different rule, not a different weighting."
            )


@dataclass(frozen=True, slots=True)
class Preset:
    """A named starting point, with the argument for it attached."""

    id: str
    title: str
    #: Why this vector, in prose, for the report and the UI. CI-enforced.
    rationale: str
    modality: Modality
    entries: tuple[WeightEntry, ...] = ()
    #: Wet-lab advice that travels with the design. Kept as data rather than
    #: prose in a template so the report cannot claim more than the preset does.
    strain_protocol: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.rationale.strip():
            raise PresetError(
                f"{self.id}: a preset must carry a rationale. For most users the "
                f"default vector is the only claim BT5 makes about what matters."
            )
        refs = [e.brief_ref for e in self.entries]
        duplicated = sorted({r for r in refs if refs.count(r) > 1})
        if duplicated:
            raise PresetError(f"{self.id}: brief_ref weighted more than once: {duplicated}")

    @property
    def by_ref(self) -> Mapping[str, WeightEntry]:
        return {e.brief_ref: e for e in self.entries}


@dataclass(frozen=True, slots=True)
class ResolvedPreset:
    """A preset bound to the rules that actually exist in this build."""

    preset: Preset
    #: spec_id -> weight, ready for DesignContext.weights.
    weights: Mapping[str, float] = field(default_factory=dict)
    #: brief_refs this preset weights that no rule implements yet. Reported, not
    #: swallowed: an objective the user believes is being optimised and that is
    #: silently absent is the difference between a ranking that means something
    #: and one that does not.
    unimplemented: tuple[str, ...] = ()

    @property
    def degradations(self) -> tuple[str, ...]:
        """Strings for Provenance.degradations."""
        if not self.unimplemented:
            return ()
        return (
            f"preset {self.preset.id} weights {len(self.unimplemented)} objective(s) "
            f"with no rule in this build: {', '.join(self.unimplemented)}",
        )


def _index_by_ref(specs: Iterable[type[Spec]]) -> Mapping[str, type[Spec]]:
    index: dict[str, type[Spec]] = {}
    for spec in specs:
        ref = spec.brief_ref
        if ref in index:
            raise PresetError(
                f"two rules claim brief_ref {ref!r}: {index[ref].id} and {spec.id}. "
                f"Presets key on brief_ref, so a duplicate makes a preset ambiguous."
            )
        index[ref] = spec
    return index


#: The slot roles `ContextSlot` admits, read from the contract rather than
#: retyped, so a role added to `core` is probed without editing this file.
_PROBE_ROLES: tuple[SlotRole, ...] = get_args(SlotRole)


def _slots_admitted_by(modality: Modality) -> tuple[ContextSlot, ...]:
    """Every slot this modality could occur in, for asking `enforcement_for`.

    A preset pins a MODALITY, not a slot: `Preset.modality` is the whole of
    what it knows about the context its weights will run in. `ContextSlot`
    cannot be built from that alone -- `role` and `host` have no defaults, and
    `__post_init__` locks `table_id` to the host, because a guessed host is a
    guessed genetic code (CLAUDE.md 3.1).

    So this enumerates rather than guesses. Naming one representative host per
    modality would answer correctly today -- no rule in the catalog keys
    enforcement on `host` or `role` -- and would go silently wrong the first
    time one does, which is the same failure this whole guard exists to fix:
    asking a narrower question than the one that decides the answer.
    Enumerating makes the guard's claim exactly "not scored in ANY slot this
    modality admits", and errs toward refusing a weight. That direction is
    deliberate: a refused weight costs a preset an objective and says so
    loudly, while a permitted one puts a constraint the validator already
    enforces into the weighted sum, silently.

    Pairs come from `LOCKED_TRANSLATION_TABLE`, which IS the declared set of
    host/table combinations, so no table is defaulted here either.
    """
    return tuple(
        ContextSlot(role=role, host=host, modality=modality, table_id=table_id)
        for host, table_id in LOCKED_TRANSLATION_TABLE.items()
        for role in _PROBE_ROLES
    )


def _unscored_enforcement(spec: type[Spec], modality: Modality) -> Enforcement | None:
    """The enforcement keeping `spec` out of the weighted sum, or None if scored.

    Asks BOTH questions and returns the first answer that bars a weight. The
    class-level `enforcement` is only the FLOOR; `enforcement_for(slot)` is
    what the solver actually routes on (`solver/catalog.py:_enforcement_of`),
    and the two disagree exactly where it matters -- `d4_internal_polya`'s
    ClassVar is SOFT while its `enforcement_for` returns HARD_REPAIR on every
    packaged modality. Reading the ClassVar alone is how two shipped presets
    came to weight a rule that is hard in the modality they are pinned to.

    A rule that `gate`s off in a slot is skipped: it contributes nothing to
    that slot's sum, so it is not being weighted there either.
    """
    if not spec.enforcement.is_scored:
        return spec.enforcement
    try:
        rule = spec()
    except TypeError:
        # A spec needing constructor arguments cannot be probed per slot. The
        # ClassVar floor above is then all this guard can honestly assert.
        return None
    enforcement_for = getattr(rule, "enforcement_for", None)
    if enforcement_for is None:
        return None
    gate = getattr(rule, "gate", None)
    for slot in _slots_admitted_by(modality):
        if gate is not None and not gate(slot):
            continue
        # Annotated because `enforcement_for` came through getattr and is
        # therefore Any: a rule returning something that is not an Enforcement
        # would otherwise slip through this guard untyped.
        level: Enforcement = enforcement_for(slot)
        if not level.is_scored:
            return level
    return None


def resolve(preset: Preset, specs: Iterable[type[Spec]] | None = None) -> ResolvedPreset:
    """Bind a preset's brief_refs to the spec ids present in this build.

    Raises when the preset weights a rule that is not SOFT **in the modality
    the preset is pinned to**. That is not a preference: `Enforcement` exists
    so that a hard constraint can never be enforced by a penalty weight, and a
    preset assigning weight to a HARD_LATTICE or HARD_REPAIR rule is precisely
    the move it forbids -- it would put a term in the weighted sum for
    something the solver has already guaranteed, and imply to the user that
    the guarantee is a trade-off.

    Enforcement is per SLOT, not per class, so this asks per slot. Reading the
    class-level `enforcement` was a guard that could not see the case it was
    written for: a rule whose floor is SOFT and whose `enforcement_for`
    escalates. `preset.modality` is what makes the question answerable here
    without a `DesignContext` -- which could not be required in any case,
    since `ResolvedPreset.weights` is an INPUT to `DesignContext.weights`.
    """
    if specs is None:
        # Autodiscovery is idempotent, and skipping it here would make an
        # un-imported catalog look identical to an unimplemented objective.
        discover()
        specs = all_specs()
    index = _index_by_ref(specs)
    weights: dict[str, float] = {}
    unimplemented: list[str] = []

    for entry in preset.entries:
        spec = index.get(entry.brief_ref)
        if spec is None:
            unimplemented.append(entry.brief_ref)
            continue
        unscored = _unscored_enforcement(spec, preset.modality)
        if entry.weight and unscored is not None:
            raise PresetError(
                f"{preset.id} gives {spec.id} ({entry.brief_ref}) weight "
                f"{entry.weight}, but it is {unscored.value} in the "
                f"{preset.modality.value} modality this preset is pinned to. "
                f"Hard constraints are guaranteed by construction or by repair "
                f"plus the independent validator, never by a penalty weight. "
                f"Use steering_weight on the rule if the DP needs a nudge."
            )
        if entry.weight != spec.default_weight and not entry.note.strip():
            raise PresetError(
                f"{preset.id} overrides {spec.id}'s default weight "
                f"({spec.default_weight} -> {entry.weight}) with no note saying "
                f"why. The rule's weight_provenance explains its own default; it "
                f"cannot explain this preset's departure from it."
            )
        weights[spec.id] = entry.weight

    return ResolvedPreset(
        preset=preset,
        weights=weights,
        unimplemented=tuple(unimplemented),
    )


# --- the shipped presets -------------------------------------------------
#
# Weights are deliberately sparse. An objective absent from a preset is not
# weighted zero by accident -- it is not part of that preset's claim, and
# `unimplemented` will say so if a rule for it exists and the preset ignores it.

_REPEAT_NOTE = (
    "recA- strains (Stbl3, NEB Stable) are standard for this workflow and "
    "suppress only the RecA-DEPENDENT pathway, which needs >~200-300 bp of "
    "homology -- the LTRs and ITRs themselves. The repeats codon choice creates "
    "or removes are 15-100 bp, squarely inside the RecA-INDEPENDENT regime the "
    "strain does not touch: a 28 bp repeat pair still recombined at 7.8e-7 to "
    "3.1e-5 in four different recA- strains (Oliveira 2008). So short-repeat "
    "avoidance is BT5's job and cannot be delegated to the strain, and it "
    "carries more weight here than in a preset with no packaging step."
)

# Why neither packaging preset weights 2.D4 (internal polyA), which is the
# hazard that most defines them.
#
# Because in these modalities it is not a trade-off. An internal polyA signal
# raised expression 3-6.5x and cut FUNCTIONAL TITER 8-9x with CMV or EF1a
# (PubMed 18627247): the genome is truncated before packaging, so the assay a
# user runs says the construct got BETTER while the thing they need collapsed.
# `d4_internal_polya.enforcement_for` therefore returns HARD_REPAIR on every
# packaged modality -- removed by Tier-B repair and proven by the independent
# validator, which refuses to emit rather than pricing a hexamer against
# expression. A weight on top of that would add a term for something already
# guaranteed and tell the user it is a slider (CLAUDE.md 3.5).
#
# These presets DID carry `WeightEntry("2.D4", 1.0, ...)`. It survived because
# `resolve()` read d4's class-level `enforcement`, which is SOFT -- the floor
# for the plasmid case, where the same hexamer costs a little expression and
# nothing else -- instead of `enforcement_for(slot)`, which is what the solver
# routes on. The guard now asks per slot, so this cannot be reintroduced.

_NATIVE_NOTE = (
    "Weighted below the mechanical objectives on purpose. Nine benchmarked "
    "optimizers were a coin flip against native sequence (Ranaghan 2021), and an "
    "18-glycoprotein Expi293F benchmark found codon optimization of human "
    "proteins in a human line did not increase yields at all. Where the evidence "
    "is that weak, the honest default is to spend the sequence's freedom on the "
    "constraints that are mechanically real."
)

LENTIVIRAL: Preset = Preset(
    id="lentiviral_hek293",
    title="Lentiviral vector, HEK293 packaging",
    rationale=(
        "The compound case BT5 exists for: propagated in E. coli, packaged in a "
        "producer line, expressed in a target cell. Weight goes to the short "
        "repeats the recA- strain does not cover, which are mechanically real "
        "and which codon choice actually controls, rather than to expression "
        "proxies whose published ceiling is ~14% of protein-level variance. "
        "The hazard that most defines this modality -- an internal polyA "
        "signal on the packaged strand -- carries no weight here precisely "
        "because it is not a trade-off: here it is hard, removed by repair and "
        "proven by the independent validator, which refuses to emit rather "
        "than pricing a truncated genome against expression."
    ),
    modality=Modality.LENTIVIRAL,
    entries=(
        WeightEntry("2.F2", 1.0, _REPEAT_NOTE),
        WeightEntry("2.C1", 0.2, _NATIVE_NOTE),
        WeightEntry("2.C3", 0.3, _NATIVE_NOTE),
    ),
    strain_protocol=(
        "Propagate in a recA- strain (Stbl3 or NEB Stable) at 30 C. This "
        "protects the LTRs, which are long perfect direct repeats -- Stbl2 to "
        "Stbl3 alone rescued an HIV vector lost entirely in 0.5 L Stbl2 "
        "cultures. It does NOT cover the 15-100 bp repeats codon choice "
        "controls; those are handled in the design, not the strain.",
    ),
)

AAV: Preset = Preset(
    id="aav_hek293",
    title="AAV vector, HEK293 packaging",
    rationale=(
        "As lentiviral, plus a hard size budget: the ITR-to-ITR payload is the "
        "binding constraint and no codon choice moves it. The one controlled "
        "dataset on AAV stuffer composition found a designed stuffer at GC "
        "43.5-44.8% cut yield up to 68% while a LOWER-GC natural stuffer of the "
        "same length cost nothing, so repetitiveness is weighted and GC is not "
        "steered toward any target on titer grounds."
    ),
    modality=Modality.AAV,
    entries=(
        WeightEntry("2.F2", 1.0, _REPEAT_NOTE),
        WeightEntry("2.C1", 0.2, _NATIVE_NOTE),
        WeightEntry("2.C3", 0.3, _NATIVE_NOTE),
    ),
    strain_protocol=(
        "Propagate in a recA- strain; for ITR palindromes prefer a sbcC-deficient "
        "strain at 42 C. As above, this covers the long repeats only.",
    ),
)

BACTERIAL: Preset = Preset(
    id="ecoli_expression",
    title="E. coli expression",
    rationale=(
        "The one context where a computable objective earns real weight. Kudla "
        "2009 measured 5' folding free energy over the -4..+37 window explaining "
        "44% of expression variance across 154 synonymous GFP variants (59% in a "
        "second promoter system), against CAI's r = 0.14, not significant. So B1 "
        "carries the highest weight in any BT5 preset -- and it is the objective "
        "that REQUIRES the vector's real 5'UTR, since the window spans the "
        "UTR/CDS junction and cannot be computed from the CDS alone."
    ),
    modality=Modality.BACTERIAL_EXPRESSION,
    entries=(
        WeightEntry(
            "2.B1",
            1.0,
            "The highest-weighted objective in BT5, and the only one whose "
            "published effect size justifies it: r = 0.66 over 154 measured "
            "variants. This MATCHES the rule's own default rather than "
            "departing from it: b1_five_prime gates to bacterial expression "
            "and nothing else, so its default weight and its weight here are "
            "necessarily the same judgement, and there is no second context "
            "for the default to have been set for.",
        ),
        WeightEntry("2.C1", 0.2, _NATIVE_NOTE),
        WeightEntry("2.C3", 0.3, _NATIVE_NOTE),
        WeightEntry("2.F2", 0.5, "Propagation still happens; no packaging step competes."),
    ),
)

PRESETS: tuple[Preset, ...] = (LENTIVIRAL, AAV, BACTERIAL)

_BY_MODALITY: Mapping[Modality, Preset] = {p.modality: p for p in PRESETS}


def preset_for(modality: Modality) -> Preset | None:
    """The shipped preset for a modality, or None when there is not one yet.

    None rather than a fallback: handing a lentiviral weight vector to an IVT
    mRNA design because it was the nearest thing on the shelf is worse than
    saying there is no curated default for that modality.
    """
    return _BY_MODALITY.get(modality)


def get(preset_id: str) -> Preset:
    for p in PRESETS:
        if p.id == preset_id:
            return p
    raise PresetError(f"no preset {preset_id!r}; have {[p.id for p in PRESETS]}")
