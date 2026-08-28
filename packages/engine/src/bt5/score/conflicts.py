"""Two rules that cannot both be satisfied here.

docs/PLAN.md makes the conflict panel a first-class output rather than an error
path, and the data model is what makes it nearly free: every Breach already
carries an interval and a slot_role, so two rules demanding opposite things over
the same bases produce two Breaches over the same Interval. Finding them is a
grouping problem, not an inference problem.

There are two kinds and they need different machinery, which is why
`Conflict.kind` distinguishes them:

- POSITIONAL conflicts are DISCOVERED, by overlap. NcoI at the start codon
  against a Kozak rule is visible because both fired over the same bases.
- STRUCTURAL conflicts are DECLARED, via `Spec.conflicts_with`, because overlap
  alone misses them. CpG depletion against a vendor GC floor is a real conflict
  that may never produce two breaches in the same window -- it is a conflict
  about the whole sequence's composition, and waiting for a positional collision
  to surface it means never surfacing it.

Overlap is computed wrap-aware on the assembled construct. A conflict at the
origin is still a conflict, and the linear predicate would report the plasmid as
having no conflicts precisely where it has one.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from bt5.core.result import Conflict
from bt5.core.spec import Breach, Enforcement, Evaluation, Spec
from bt5.core.types import Interval


@dataclass(frozen=True, slots=True)
class _Group:
    interval: Interval
    breaches: tuple[Breach, ...]


def _merge(a: Interval, b: Interval, *, length: int, circular: bool) -> Interval:
    """The smallest interval covering both, in the canonical representation.

    On a circular construct the two may be expressed a turn apart, so `b` is
    first brought into `a`'s frame. Taking a plain min/max of the raw endpoints
    would produce an interval spanning almost the whole plasmid, and the
    conflict panel's whole value is that its intervals are small enough to point
    at something.
    """
    if not circular:
        return Interval(min(a.start, b.start), max(a.end, b.end))
    best: Interval | None = None
    for shift in (0, length, -length):
        start, end = b.start + shift, b.end + shift
        if not (a.start < end and start < a.end):
            continue
        candidate_start = min(a.start, start)
        candidate_end = max(a.end, end)
        if candidate_start < 0:
            candidate_start += length
            candidate_end += length
        candidate = Interval(candidate_start, candidate_end)
        if best is None or candidate.length < best.length:
            best = candidate
    return best if best is not None else Interval(min(a.start, b.start), max(a.end, b.end))


def _positional_groups(
    breaches: Sequence[Breach], *, length: int, circular: bool
) -> tuple[_Group, ...]:
    """Cluster breaches into overlapping runs, regardless of which rule made them."""
    groups: list[_Group] = []
    for breach in sorted(breaches, key=lambda b: (b.interval.start, b.interval.end, b.spec_id)):
        placed = False
        for i, group in enumerate(groups):
            if group.interval.overlaps(breach.interval, length, circular):
                groups[i] = _Group(
                    interval=_merge(
                        group.interval, breach.interval, length=length, circular=circular
                    ),
                    breaches=(*group.breaches, breach),
                )
                placed = True
                break
        if not placed:
            groups.append(_Group(interval=breach.interval, breaches=(breach,)))

    # A group can grow past a group formed later -- most easily around the
    # origin, where sorting by start says nothing about adjacency. Two groups
    # reported where the bases are one contiguous problem would split one
    # conflict into two, each naming half the rules involved, so settle it.
    merged = True
    while merged:
        merged = False
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                if groups[i].interval.overlaps(groups[j].interval, length, circular):
                    groups[i] = _Group(
                        interval=_merge(
                            groups[i].interval, groups[j].interval, length=length, circular=circular
                        ),
                        breaches=(*groups[i].breaches, *groups[j].breaches),
                    )
                    del groups[j]
                    merged = True
                    break
            if merged:
                break
    return tuple(groups)


def _kind(breaches: Sequence[Breach]) -> str:
    """Which sort of positional conflict this is.

    `immutable_region` wins whenever any member is unfixable by codon choice:
    that changes what the user can DO about it, from "relax one of these rules"
    to "no codon anywhere resolves this", and it is the more useful thing to say.
    """
    if any(not b.fixable_by_codon_choice for b in breaches):
        return "immutable_region"
    sides = {str(b.detail.get("binding_side", "")) for b in breaches}
    sides.discard("")
    if len(sides) > 1:
        return "opposing_gradient"
    return "mutually_exclusive"


def _binding(breaches: Sequence[Breach], specs: Mapping[str, type[Spec]]) -> str:
    """Which rule is actually forcing the issue here.

    A hard rule beats a soft one outright -- the soft one is a preference and the
    hard one is not, so naming the soft rule as binding would point the user at
    the slider that cannot help. Among equals, the largest magnitude.
    """

    def rank(b: Breach) -> tuple[int, float]:
        spec = specs.get(b.spec_id)
        hard = 1 if spec is not None and spec.enforcement.is_hard else 0
        return (hard, b.magnitude)

    return max(breaches, key=rank).spec_id


def detect_conflicts(
    evaluations: Iterable[Evaluation],
    specs: Mapping[str, type[Spec]],
    *,
    length: int,
    circular: bool,
) -> tuple[Conflict, ...]:
    """Every conflict the evidence supports, positional and declared.

    `specs` maps spec_id to the rule class, and is passed in rather than read
    from the registry so this is a pure function of its inputs -- the conflict
    panel is a reported artefact and it must not change because an unrelated
    import registered another rule.
    """
    breaches = [b for ev in evaluations for b in ev.breaches]
    out: list[Conflict] = []

    for group in _positional_groups(breaches, length=length, circular=circular):
        ids = sorted({b.spec_id for b in group.breaches})
        if len(ids) < 2:
            continue  # one rule unhappy in one place is a breach, not a conflict
        out.append(
            Conflict(
                interval=group.interval,
                spec_ids=tuple(ids),
                kind=_kind(group.breaches),  # type: ignore[arg-type]
                binding_spec_id=_binding(group.breaches, specs),
            )
        )

    out.extend(_declared(breaches, specs, length=length))
    return tuple(out)


def _declared(
    breaches: Sequence[Breach],
    specs: Mapping[str, type[Spec]],
    *,
    length: int,
) -> tuple[Conflict, ...]:
    """Structural conflicts named by `Spec.conflicts_with`.

    Reported when BOTH sides are present in this run and at least one of them
    actually fired. Requiring both to fire would hide the common case -- one
    rule is satisfied precisely because it won, and the losing rule's breach is
    the only visible trace of the disagreement. Requiring neither to fire would
    report a declared conflict on every construct, which is noise.
    """
    fired = {b.spec_id for b in breaches}
    present = set(specs)
    seen: set[tuple[str, str]] = set()
    out: list[Conflict] = []

    for spec_id, spec in sorted(specs.items()):
        for other in spec.conflicts_with:
            if other not in present:
                continue
            pair = (min(spec_id, other), max(spec_id, other))
            if pair in seen or not (fired & set(pair)):
                continue
            seen.add(pair)
            spans = [b.interval for b in breaches if b.spec_id in pair]
            out.append(
                Conflict(
                    interval=spans[0] if spans else Interval(0, max(1, length)),
                    spec_ids=pair,
                    kind="declared",
                    binding_spec_id=_binding([b for b in breaches if b.spec_id in pair], specs),
                )
            )
    return tuple(out)


def hard_versus_soft(conflict: Conflict, specs: Mapping[str, type[Spec]]) -> bool:
    """Is this a hard rule overriding a soft one?

    Worth separating in the UI: it is not a trade-off the user can slide their
    way out of, and presenting it beside genuine trade-offs invites them to try.
    """
    levels = {specs[s].enforcement for s in conflict.spec_ids if s in specs}
    return any(e.is_hard for e in levels) and Enforcement.SOFT in levels
