"""One registry of orderable vendor configurations.

**Why this file exists is a bug that already shipped.** Vendor facts lived in two
places that each validated their own keys: `fragment.VENDOR_ADAPTERS` /
`VENDOR_LENGTHS` keyed on `twist_gene_fragment`, and
`e1_homopolymers.VENDOR_LIMITS` keyed on `twist_standard` and `idt_gblocks`.
Both lookups succeeded, so a single run could be spec'd for IDT's homopolymer
limits and Twist's length range at the same time and nothing anywhere noticed --
the two dicts did not even agree on what a key *was*, since `twist_standard` is
a Twist complexity TIER label and not a product at all. E1 defaulted to IDT
while E4-E7 and E9 defaulted to Twist, which is the same incoherence wearing a
different hat: with no vendor chosen, BT5's answer was "IDT's runs, Twist's
lengths".

So: one dataclass, one dict, one default. A `VendorProfile` is everything BT5
knows about ONE orderable configuration, and a rule that wants a vendor fact
takes the profile rather than reaching into a table keyed on a string it hopes
matches.

**Every value carries where it came from.** Vendor numbers drift -- Twist moved
its published homopolymer limit from 14 to 30 bp between 2023 and 2026 -- so a
profile without provenance is a number that goes quietly wrong. `last_verified`
is per profile, `notes` says which values are measured and which are published,
and `_provenance.json` beside this file records the source of each.

**Inheritance is allowed within a vendor and never across one.** The GC bands
and the IDT homopolymer limits were established for ONE product line each
(gBlocks plate entry; the Twist Gene Fragment order form). Where a sibling
product from the same vendor has no measurement of its own it inherits, and the
`notes` field says so out loud. Nothing is ever carried from one vendor to
another: the whole finding of the 18-probe ladder is that Twist ships 80% GC as
Standard where IDT denies it, so a value transferred across vendors would be
transferring the one thing that was measured to differ.
"""

from __future__ import annotations

from dataclasses import dataclass

from bt5.rules.fragment import (
    IDT_EBLOCKS,
    IDT_GBLOCKS,
    NO_ADAPTERS,
    TWIST_ADAPTER_ON,
    TWIST_GENE_FRAGMENT,
    Adapters,
)


@dataclass(frozen=True, slots=True)
class VendorProfile:
    """One configuration a user can actually place an order for.

    The orderable fields are all `| None` together, and `__post_init__` refuses a
    profile that fills in some and not others. That is what makes `none` -- "no
    vendor chosen" -- expressible as a member of this registry rather than as a
    key that happens to be absent from one dict and present in another, which is
    precisely how the two-namespace bug hid.
    """

    #: Registry key. Also the string a rule reports in a finding, so it has to
    #: name a product a person can order, not a tier or a vendor.
    key: str
    #: Human-facing vendor name, for report prose. Empty for `none`.
    vendor: str
    #: Human-facing product name. Empty for `none`.
    product: str
    #: Fixed sequence the vendor synthesises onto the fragment. Empty for every
    #: current product except the Twist adapter-on option.
    adapters: Adapters

    #: (min, max) bp of ORDERED DNA. See `e9_length_tiers` for why "ordered DNA"
    #: rather than "ordered DNA plus adapters" is an assumption and not a fact.
    length_bp: tuple[int, int] | None = None
    #: Longest accepted A/T run, then longest accepted G/C run. Two numbers and
    #: not one: the asymmetry is the vendor evidence that G/C runs are chemically
    #: worse rather than merely repetitive. See `e1_homopolymers`.
    homopolymer_at: int | None = None
    homopolymer_gc: int | None = None
    #: (min, max) fraction GLOBAL GC the vendor will accept. Measured, not
    #: published -- see docs/design/vendor-gc-calibration.md. No rule reads this
    #: yet; E2 carries the universal band and #43 V3 wires the per-vendor one.
    global_gc: tuple[float, float] | None = None

    last_verified: str = ""
    #: Which values here are measured, which are published, and which are
    #: inherited from a sibling product. Required for an orderable profile.
    notes: str = ""

    def __post_init__(self) -> None:
        # The adapters carry a vendor key of their own, and a profile whose
        # adapters name a different configuration is the two-namespace bug
        # reconstituted inside one object.
        if self.adapters.vendor != self.key:
            raise ValueError(
                f"profile {self.key!r} carries adapters labelled "
                f"{self.adapters.vendor!r}; one configuration, one name"
            )

        orderable = [
            self.length_bp,
            self.homopolymer_at,
            self.homopolymer_gc,
            self.global_gc,
        ]
        present = [v is not None for v in orderable]
        if any(present) and not all(present):
            raise ValueError(
                f"profile {self.key!r} is half-specified: an orderable configuration "
                f"needs a length range, both homopolymer limits and a GC band, and a "
                f"non-orderable one needs none of them. Partial is how a rule ends up "
                f"silently unconstrained"
            )

        if self.length_bp is not None:
            lo, hi = self.length_bp
            if not 0 < lo < hi:
                raise ValueError(f"profile {self.key!r}: length range {self.length_bp} is empty")
            assert self.global_gc is not None  # narrowed by the all-or-nothing check
            gc_lo, gc_hi = self.global_gc
            if not 0.0 <= gc_lo < gc_hi <= 1.0:
                raise ValueError(f"profile {self.key!r}: GC band {self.global_gc} is not a band")
            if not self.notes.strip():
                raise ValueError(
                    f"profile {self.key!r}: an orderable profile must say which of its "
                    f"numbers are measured, published or inherited"
                )
            if not self.last_verified.strip():
                raise ValueError(f"profile {self.key!r}: vendor numbers need a verification date")

    @property
    def is_orderable(self) -> bool:
        """Whether a person can place an order against this configuration.

        Derived rather than stored, so it cannot come to disagree with the data
        it summarises.
        """
        return self.length_bp is not None

    def accepts_length(self, bp: int) -> bool:
        """Whether `bp` of ordered DNA is inside this configuration's range."""
        if self.length_bp is None:
            return False
        lo, hi = self.length_bp
        return lo <= bp <= hi


#: Twist's published homopolymer bound, read conservatively and deliberately so.
#:
#: The brief recorded "under 14 bp", which is the 2023 number; Twist's current
#: FAQ says 30 bp. We keep the tighter reading because it is the safe direction
#: for a value we cannot pin to the Gene Fragment product specifically, and
#: because loosening a HARD_LATTICE limit on an unverified number is the one
#: change here that would be invisible if it were wrong -- the automaton would
#: simply stop forbidding runs and nothing would report it.
_TWIST_HOMOPOLYMER = (13, 13)

#: Measured against the IDT gBlocks plate-entry and Twist Gene Fragment order
#: forms, 18-sequence ladder, 2026-08-28. Both refuse global GC <=25% and accept
#: >=30%, so 0.28 is the split; at the top Twist ships 80% as Standard while
#: IDT's solved penalty line (score = 1.40*GC% - 83.8) puts Denied at 77.0%.
#: docs/design/vendor-gc-calibration.md
_TWIST_GC = (0.28, 0.80)
_IDT_GC = (0.28, 0.77)

_TWIST_NOTES = (
    "GC band MEASURED (18-probe ladder through the Twist order form, 2026-08-28): "
    "<=25% Not Accepted, 30-80% all Standard. Length range PUBLISHED. Homopolymer "
    "limits PUBLISHED but stale-leaning: the brief's 'under 14 bp' is the 2023 "
    "number and the current FAQ says 30 bp, kept tight because loosening a "
    "HARD_LATTICE bound on an unverified number fails silently."
)

_IDT_GBLOCKS_NOTES = (
    "GC band MEASURED (18-probe ladder through the gBlocks plate-entry complexity "
    "checker, 2026-08-28): <=25% and >=80% Denied, 70-75% Accepted-Moderate, "
    "30-65% clean; the ceiling is the solved Denied threshold, 77.0%. Length range "
    "and homopolymer limits PUBLISHED."
)

_IDT_EBLOCKS_NOTES = (
    "Length range PUBLISHED for eBlocks. GC band and homopolymer limits INHERITED "
    "from idt_gblocks -- same vendor, different product line, and the calibration "
    "ladder was run through gBlocks plate entry only. Inheriting within a vendor "
    "is the policy; nothing is ever carried across one, since Twist shipping 80% "
    "as Standard where IDT denies it is the finding that ladder exists to record."
)

#: Every configuration BT5 can spec a fragment for, plus the degenerate one.
#:
#: Gene fragments only, Twist and IDT only. GenScript's (15, 15) homopolymer
#: entry is deliberately dropped rather than carried: it had no adapters and no
#: length range, so it could never become a complete profile, and a half-filled
#: profile is what `__post_init__` now refuses on principle.
PROFILES: dict[str, VendorProfile] = {
    #: "No vendor chosen." Orderable fields absent, and absent TOGETHER. Rules
    #: that only need the synthesis SCOPE (E4-E7: the fragment is linear, the
    #: backbone is not synthesised) accept it; rules that ask what a vendor will
    #: take (E1, E9) refuse it, because there is no vendor to ask.
    "none": VendorProfile(
        key="none",
        vendor="",
        product="",
        adapters=NO_ADAPTERS,
    ),
    "twist_gene_fragment": VendorProfile(
        key="twist_gene_fragment",
        vendor="Twist Bioscience",
        product="Gene Fragment",
        adapters=TWIST_GENE_FRAGMENT,
        length_bp=(300, 5000),
        homopolymer_at=_TWIST_HOMOPOLYMER[0],
        homopolymer_gc=_TWIST_HOMOPOLYMER[1],
        global_gc=_TWIST_GC,
        last_verified="2026-08-28",
        notes=_TWIST_NOTES,
    ),
    "twist_gene_fragment_adapter_on": VendorProfile(
        key="twist_gene_fragment_adapter_on",
        vendor="Twist Bioscience",
        product="Gene Fragment, adapters appended",
        adapters=TWIST_ADAPTER_ON,
        length_bp=(300, 5000),
        homopolymer_at=_TWIST_HOMOPOLYMER[0],
        homopolymer_gc=_TWIST_HOMOPOLYMER[1],
        global_gc=_TWIST_GC,
        last_verified="2026-08-28",
        notes=(
            _TWIST_NOTES + " Adapter-on is a checkout option on the same product, so "
            "every number here is the plain Gene Fragment's, unchanged."
        ),
    ),
    "idt_gblocks": VendorProfile(
        key="idt_gblocks",
        vendor="IDT",
        product="gBlocks Gene Fragment",
        adapters=IDT_GBLOCKS,
        length_bp=(125, 3000),
        homopolymer_at=9,
        homopolymer_gc=5,
        global_gc=_IDT_GC,
        last_verified="2026-08-28",
        notes=_IDT_GBLOCKS_NOTES,
    ),
    "idt_eblocks": VendorProfile(
        key="idt_eblocks",
        vendor="IDT",
        product="eBlocks Gene Fragment",
        adapters=IDT_EBLOCKS,
        length_bp=(300, 1500),
        homopolymer_at=9,
        homopolymer_gc=5,
        global_gc=_IDT_GC,
        last_verified="2026-08-28",
        notes=_IDT_EBLOCKS_NOTES,
    ),
}

#: The one default, and it is IDT's gBlocks rather than Twist's Gene Fragment
#: for a reason that is not a coin flip.
#:
#: Unifying the split default has to move one rule or the other, and the two
#: directions are not symmetric. Defaulting to Twist would take E1's run limits
#: from (9, 5) to (13, 13) -- a HARD_LATTICE loosening, which means a 13 nt G-run
#: simply stops being unreachable and NOTHING REPORTS IT, because a lattice
#: constraint's whole output is the absence of a finding. Defaulting to IDT
#: instead moves E9's length range from 300-5000 to 125-3000. That is a real
#: change, but it is a visible one -- E9 is HARD_CHECK and every finding it emits
#: already names the other configurations that would take the fragment -- and it
#: is the CORRECT answer to the question actually being asked, since 200 bp is
#: genuinely orderable as a gBlock.
#:
#: A synthetic "strictest of every vendor" profile was considered and rejected:
#: it would name a product nobody can order, which is exactly what E9's findings
#: must not do.
#:
#: This default is what BT5 assumes when the user has not chosen. #43 V3 replaces
#: the assumption with the user's actual selection, at which point the choice
#: here stops mattering.
DEFAULT_VENDOR = "idt_gblocks"


def profile(key: str) -> VendorProfile:
    """Look up a configuration, with an error that lists the real ones."""
    try:
        return PROFILES[key]
    except KeyError:
        raise ValueError(f"unknown vendor {key!r}; have {sorted(PROFILES)}") from None


def orderable(key: str) -> VendorProfile:
    """Look up a configuration that can actually be ordered from somebody.

    For rules whose whole question is "will a vendor take this" -- E1's run
    limits, E9's length range. `none` reaches them as a request to answer a
    question with no answer, and guessing a vendor to answer it with is how the
    split default happened in the first place.
    """
    p = profile(key)
    if not p.is_orderable:
        raise ValueError(
            f"{key!r} is not orderable from anyone, so it has no vendor limits to "
            f"check against; choose one of {sorted(orderable_keys())}"
        )
    return p


def orderable_keys() -> tuple[str, ...]:
    """Every configuration a fragment can actually be ordered as, sorted."""
    return tuple(sorted(k for k, p in PROFILES.items() if p.is_orderable))


def all_keys() -> tuple[str, ...]:
    """Every configuration including `none`, sorted. The `param_schema` enum."""
    return tuple(sorted(PROFILES))


def accepting_length(bp: int, exclude: str = "") -> tuple[str, ...]:
    """Every orderable configuration whose range contains `bp`, sorted.

    The advisory half of E9: a finding no codon can act on is only useful if it
    names somewhere the order could go instead.
    """
    return tuple(
        sorted(k for k in orderable_keys() if k != exclude and PROFILES[k].accepts_length(bp))
    )


@dataclass(frozen=True, slots=True)
class VendorSelection:
    """One or more configurations a fragment is being spec'd for at once.

    #43 V2/V3. A rule used to take a single `vendor: str`; it now takes a
    selection, because a user may want a design that every one of several vendors
    can build. The selection carries the answer to the three questions that
    genuinely differ per vendor -- run limits (E1), length range (E9), GC band
    (E2, wired in #43 V3b) -- and nothing else, because with adapters off the
    table (the owner never orders them) the ordered molecule is identical under
    every configuration and the four synthesis rules see no difference at all.

    `of()` is the ONLY public constructor, and it is varargs on purpose:
    `VendorSelection.of("idt_gblocks")` is one key and can never be misread as
    five characters. `require_selection` refuses a bare string outright, so the
    old `vendor="idt_gblocks"` call site fails loudly rather than silently
    iterating a string.

    Multi-vendor is N verdicts, not one merged threshold. But it is still ONE
    evaluation and ONE breach per physical fact -- the rules attribute a finding
    to the configurations it is a finding against, they do not re-run themselves
    per vendor. Where two vendors' limits differ, enforcement is the stricter
    (there is one automaton, one lattice) and only the ATTRIBUTION is per vendor.
    """

    #: Deduplicated, order-preserving. `of()` is the only thing that builds it.
    keys: tuple[str, ...]

    @staticmethod
    def of(*keys: str) -> VendorSelection:
        """Build a selection, refusing every combination that is not a design.

        Four refusals, each reusing the registry's own message where one exists:
        an empty selection is not a choice; an unknown key is `profile()`'s
        error; `none` cannot be combined with a real vendor; and adapter-on may
        not be mixed with adapter-free, because those are physically different
        molecules and BT5 designs one. That last refusal is what makes the
        "different molecule per vendor" hazard structurally unreachable rather
        than something each rule has to remember to handle.
        """
        if not keys:
            raise ValueError(
                "an empty selection is not a choice; omit the argument for the "
                "default or name a configuration"
            )
        seen: dict[str, None] = {}
        for k in keys:
            profile(k)  # raises "unknown vendor 'x'; have [...]" -- one message, reused
            seen[k] = None
        unique = tuple(seen)

        if len(unique) > 1 and "none" in unique:
            raise ValueError(
                "'no vendor chosen' is not a vendor you can add another to; pick "
                "configurations or pick none, not both"
            )
        # Compare the PAYLOAD, not the Adapters object: NO_ADAPTERS,
        # TWIST_GENE_FRAGMENT, IDT_GBLOCKS and IDT_EBLOCKS are four distinct
        # objects that all carry the same empty payload, so object identity would
        # refuse a selection of two adapter-free products it should accept.
        payloads = {(profile(k).adapters.five, profile(k).adapters.three) for k in unique}
        if len(payloads) > 1:
            raise ValueError(
                "a selection may not mix adapter-on and adapter-free configurations: "
                "they are physically different molecules and BT5 designs one. Order "
                "them as separate designs"
            )
        return VendorSelection(unique)

    @property
    def profiles(self) -> tuple[VendorProfile, ...]:
        return tuple(profile(k) for k in self.keys)

    @property
    def label(self) -> str:
        """The configurations, comma-joined -- what a finding names. For a single
        key this is the key itself, so a single-vendor finding reads exactly as it
        did before the selection existed."""
        return ", ".join(self.keys)

    def orderable_only(self) -> VendorSelection:
        """Assert every member is a product a vendor will take, and return self.

        For rules whose whole question is "will a vendor build this" -- E1's run
        limits, E9's length range. `none` reaches them as a request to answer a
        question with no answer; `orderable()` raises for it with the message the
        registry already owns. Returns self so it composes at a call site.
        """
        for k in self.keys:
            orderable(k)
        return self

    @property
    def adapters(self) -> Adapters:
        """The one adapter payload the selection shares, labelled with the whole
        selection. `of()` proved the payload is shared, so any member's is the
        selection's. The label -- not a member's key -- is what `fragments()`
        stamps onto `Fragment.vendor`, so a finding names the molecule's whole
        selection rather than one arbitrary member of it. Single key: identical
        to that product's own adapters today."""
        first = self.profiles[0].adapters
        return Adapters(first.five, first.three, vendor=self.label)

    def homopolymer_limits(
        self,
    ) -> tuple[tuple[int, tuple[str, ...]], tuple[int, tuple[str, ...]]]:
        """Per axis, the strictest run limit across the selection and the real
        configurations that publish it.

        Enforcement is min-merged and cannot be otherwise: there is one automaton
        and `lattice_terms` returns one forbidden set, so the union of forbidden
        runs IS the minimum. What the binders buy is attribution -- naming which
        selected product a run is too long for. Requires an orderable selection
        (call `orderable_only()` first); an orderable profile always carries both
        axes.
        """

        def axis(values: list[tuple[int, str]]) -> tuple[int, tuple[str, ...]]:
            m = min(v for v, _ in values)
            return m, tuple(k for v, k in values if v == m)

        at: list[tuple[int, str]] = []
        gc: list[tuple[int, str]] = []
        for k, p in zip(self.keys, self.profiles, strict=True):
            assert p.homopolymer_at is not None  # orderable: both axes present
            assert p.homopolymer_gc is not None
            at.append((p.homopolymer_at, k))
            gc.append((p.homopolymer_gc, k))
        return axis(at), axis(gc)

    def homopolymer_accepts(self, length: int, base_class: str) -> tuple[str, ...]:
        """Selected configurations that accept a run of `length` on this axis.

        A run over the strictest limit may still be within a looser vendor's, and
        E1 says so when the run is recodeable. `base_class` is "A/T" or "G/C".
        """
        axis = "homopolymer_at" if base_class == "A/T" else "homopolymer_gc"
        out = []
        for k, p in zip(self.keys, self.profiles, strict=True):
            limit = getattr(p, axis)
            if limit is not None and length <= limit:
                out.append(k)
        return tuple(out)

    def verdicts_for_length(self, bp: int) -> tuple[tuple[str, str], ...]:
        """Per orderable member, whether it accepts `bp` of ordered DNA.

        `"accept"` / `"ambiguous"` / `"refuse"`. Ambiguous is the insert-versus-
        insert+adapters question E9 already reported; it can only arise for an
        adapter-on selection, which `of()` guarantees is single-key. Non-orderable
        members (`none`) are skipped; E9 calls `orderable_only()` so there are
        none.
        """
        out: list[tuple[str, str]] = []
        for k, p in zip(self.keys, self.profiles, strict=True):
            if p.length_bp is None:
                continue
            lo, hi = p.length_bp
            total = bp + p.adapters.total
            ordered_ok = lo <= bp <= hi
            total_ok = lo <= total <= hi
            if ordered_ok and total_ok:
                out.append((k, "accept"))
            elif ordered_ok != total_ok:
                out.append((k, "ambiguous"))
            else:
                out.append((k, "refuse"))
        return tuple(out)

    def alternatives_for(self, bp: int) -> tuple[str, ...]:
        """Orderable configurations OUTSIDE the selection whose range contains
        `bp`. Delegates to `accepting_length` so the two cannot drift; for a
        single-key selection it equals `accepting_length(bp, exclude=key)`,
        pinned by test."""
        chosen = set(self.keys)
        return tuple(k for k in accepting_length(bp) if k not in chosen)


#: What every rule assumes when the user has not chosen. One key, the one default.
DEFAULT_SELECTION = VendorSelection.of(DEFAULT_VENDOR)


def require_selection(value: object) -> VendorSelection:
    """The one gate every rule constructor runs its `vendors` argument through.

    A `VendorSelection` passes; anything else -- above all a bare `str`, the old
    call shape -- raises `TypeError` rather than being silently iterated into a
    per-character selection. This is the structural half of what `of()`'s varargs
    started: the selection cannot be built wrong, and it cannot be passed wrong.
    """
    if not isinstance(value, VendorSelection):
        raise TypeError(
            f"vendors must be a VendorSelection, not {type(value).__name__}; "
            f"pass VendorSelection.of('idt_gblocks'), not a bare string"
        )
    return value


# Import-time proof, not a test: a registry whose keys and members disagree is
# the exact failure this module was written to end, and it should be impossible
# to import a broken one rather than merely impossible to ship it.
for _key, _profile in PROFILES.items():
    if _key != _profile.key:
        raise ValueError(f"PROFILES[{_key!r}] holds a profile named {_profile.key!r}")

__all__ = [
    "DEFAULT_SELECTION",
    "DEFAULT_VENDOR",
    "PROFILES",
    "VendorProfile",
    "VendorSelection",
    "accepting_length",
    "all_keys",
    "orderable",
    "orderable_keys",
    "profile",
    "require_selection",
]
