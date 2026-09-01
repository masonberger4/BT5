## 2026-09-01 — The §3.5 weighting guard asks per slot, and 2.D4 leaves two presets (#72)

**Decided:** `score/presets.py` `resolve()` guards on `enforcement_for(slot)` for every
slot the preset's modality admits, not on the `enforcement` ClassVar. The ClassVar is a
FLOOR — `d4_internal_polya` declares `SOFT` (right for the plasmid case) while its
`enforcement_for` returns `HARD_REPAIR` on every packaged modality — so the old guard
asked whether a rule is hard *everywhere* instead of whether it is hard *here*, and
passed. `LENTIVIRAL` and `AAV` both shipped `WeightEntry("2.D4", 1.0)` through it. Those
entries are removed; that is the actual fix, and the resolver change is what stops it
being reintroduced.

**No signature change.** `resolve()` still takes `specs`, not a `DesignContext`.
`Preset.modality` is the pin, and every `enforcement_for` in the catalog keys on
`slot.modality` alone. Requiring a `DesignContext` would also be circular:
`ResolvedPreset.weights` is an *input* to `DesignContext.weights`, so the context does
not exist yet at the moment `resolve()` runs.

**Rejected:**
- *Naming one representative host per modality to build the probe slot.* `ContextSlot`
  cannot be built from a modality alone — `role` and `host` have no defaults and
  `__post_init__` locks `table_id` to the host — so a single probe means guessing a
  host. It answers correctly today (no catalog rule keys `gate` or `enforcement_for` on
  `host` or `role`) and goes silently wrong the first time one does, which is the same
  shape as the bug being fixed: asking a narrower question than the one that decides the
  answer. `_slots_admitted_by` enumerates `LOCKED_TRANSLATION_TABLE` × `SlotRole`
  instead, so nothing is guessed and nothing is defaulted.
- *Guarding on `is_hard`.* `is_scored` is `is SOFT`, so `is_hard` would have started
  admitting `REPORT_ONLY` into the weighted sum — a weakening smuggled in as a fix.
  The guard refuses anything not scored, as before.
- *Keeping the ClassVar guard and documenting the per-slot check as the consumer's job*
  (option 2 in #72). `bt5/design/catalog.py` already does that, but a guard that has to
  be re-implemented by every consumer is not a guarantee. Left untouched — it is PR #71's
  file.
- *Deleting `_POLYA_NOTE`'s science with its entry.* The 8–9× functional titer loss is
  why d4 is HARD in these modalities; it is now a comment above the presets explaining
  why 2.D4 is deliberately absent, plus the reason the old guard let it through.

**Also corrected:** `LENTIVIRAL.rationale` claimed weight went to "internal polyA on the
packaged strand, cryptic splice donors". Neither was weighted after this change, and the
splice claim was never true in any form — there is no splice-donor rule in the catalog at
all (15 rules, none for splicing), so the prose asserted an objective that does not
exist. It now names internal polyA only, and says it carries no weight *because* it is
hard.

**Scientific impact — not "none".** d4 leaves the weighted sum for lentiviral and AAV
designs, so candidate ranking moves for both. That is the correct direction: the
constraint is enforced by Tier-B repair plus the independent validator, which refuses to
emit, and weighting it as well both double-counts it and tells the user a guarantee is a
trade-off. Owner merges under §7b.

**Where:** branch `claude/hard-rule-weighting-guard-nitkbe`; lane M3, `score/` only.
