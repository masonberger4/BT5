# BT5 — Architecture & Buildout Proposal
## Framing C: risk-first vertical slice

---

## 0. The thesis in one paragraph

The dangerous thing about this project is not that any single algorithm is hard. It is that **six of the eight hardest problems all push on the same interface** — the shape of "a thing that evaluates a rule against a sequence." Circularity, immutable backbone regions, three simultaneous contexts, strand-of-interest, per-window conflict attribution, minimal-conflicting-set reporting, and objective normalization are *all* properties of that one type. If it is wrong, eight agents' worth of rule code is wrong with it, and the error surfaces only when the vector lane lands in week four. So the plan is: **spend week one building a working, ugly, one-host, three-rule, one-plasmid product that exercises every one of those pressures at once**, run seven numeric go/no-go gates against it, then freeze. Wave 1 is three agents on the three lanes that can still force a core change. Only after wave 1 lands does the interface freeze become real and the remaining seven lanes start.

---

## 1. Risk register, ranked

Ranked by *expected damage × probability × how late you find out*. The last factor is what Framing C exists to attack.

### Kill-class — these end the project or reduce it to a DnaChisel reskin

| # | Risk | Why fatal | Detected by |
|---|---|---|---|
| **R1** | **The `Spec`/`Evaluation` interface cannot carry localized, attributable, context-gated, circular-coordinate results.** | Every rule in lanes 3, 4, 5, 6, 7 implements it. A change in week 4 invalidates ~60% of engine work. Nothing else in the system has this blast radius, and it is the risk *not* on the user's list. | Gate **G0** (below), and by making wave 1 = the three lanes that stress it hardest |
| **R2** | **Circular assembled construct with immutable backbone regions is not expressible in the mutation space.** | This is *the* competitive differentiator (brief §6.1). It is not bolt-on-able: it changes what a "position" means, what the Aho–Corasick seed state is, what a window is, and what infeasibility means. Lose it and BT5 is a free-floating-CDS optimizer, i.e. everything else on the market. | G2 |
| **R3** | **Tier-A automaton DP with windowed-GC state is intractable or does not actually guarantee.** | Hard constraints stop being guaranteed by construction; you fall back to reject-and-repair, which is DnaChisel. Survivable but it removes the correctness story. | G1 |
| **R4** | **The candidate gallery degenerates.** A λ-sweep over a weighted sum recovers only *supported* (convex-hull) Pareto points. With hard constraints doing most of the work, 24 λ-points can collapse to 3 near-identical sequences. | The entire trade-off UI premise (locked decision #3) and the honesty story (brief §4, §8.11) rest on showing *genuinely different* designs. If they are all 98% identical, the product's central claim is false and you find out after the UI is built. **Nobody named this risk.** | G4 |

### Severe — schedule and trust damage, recoverable

| # | Risk | Note |
|---|---|---|
| **R5** | **Infeasibility with a minimal conflicting set.** critique §Tier-5-9: over-constrained failure is the *dominant* real-world experience. "No solution" is a useless product. | G5 |
| **R6** | **Objective normalization against a random-synonymous null.** critique Tier-1 #3 — sliders dead over most of their range is the single most likely reason the app feels broken. Also a latency risk: 200–500 extra evaluations per job. | G3 |
| **R7** | **Folding inside a 10 s interactive loop**, and the seqfold-vs-ViennaRNA ΔG offset (brief §7-tech-iii names this the most likely correctness bug in the whole feature). | G6 |
| **R8** | **GenBank feature round-tripping** — 1-based inclusive vs 0-based half-open, origin-spanning `join()`, `SeqRecord` slicing silently dropping features and all `.annotations`. | G2 |

### Contained — real, but discoverable late without much damage

R9 scale ceiling (FVIII 2332 aa, dystrophin 3685 aa); R10 biosecurity DB install friction; R11 cross-platform determinism; R12 vendor-rule staleness (Twist changed homopolymer 14→30); R13 CI concurrency starvation at 8 agents (githubPlan §8.1 — you are 4–5× oversubscribed).

**The single most important sentence in this document:** R1 is an *interface* risk, not an algorithm risk, and interface risks are the only kind that parallel agents cannot absorb. Everything in §5 exists to retire R1 before agent #4 starts.

---

## 2. The walking skeleton — what exists at the end of week one

The foundation is **not scaffolding**. It is a working product with the dial turned to 1 on every axis:

- **one host** (*E. coli* K-12, NCBI table 11) — not a host framework
- **three rules** — `AvoidPatterns` (Tier-A capable), `GCBand` (window-reporting, Tier-B only), `Homopolymer` (trivially both)
- **one vector** (a real Addgene lentiviral transfer plasmid GenBank, with an origin-spanning feature)
- **one folding backend** (seqfold) behind the `FoldEngine` protocol, with one real structure term (−4…+37 accessibility)
- **one preset**, real λ-sweep, real 200-variant null model, real percentile normalization
- **one report** (JSON QC + annotated GenBank export) and one UI page

and it goes **protein → validated → biosecurity-screened → cassette planned → spliced into the circular backbone → mutation space over editable intervals only → Tier-A DP → Tier-B repair → independent `verify_construct` on the assembled circular product → normalized scorecard → 5-candidate gallery → annotated GenBank + QC report** in under 10 seconds.

Every one of R1–R8 is touched. That is the point. Everything after week one **thickens** this path; nothing after week one *creates* a stage of it.

### Foundation lands as two merges, not one

| Merge | Contents | Rulesets |
|---|---|---|
| **PR #0 — Contract** | `bt5/core/**` (frozen types, `Spec` protocol, registry, services, result types), `verify.py`, CI workflows, `path_guard.py`, `CLAUDE.md`, labels, `data/**`, `benchmarks/{panel,tolerances}` with **every metric pre-registered**, `tests/invariants/**`, `pyproject.toml` with **every third-party dependency pre-declared**, `uv.lock` | pushed to `main` *before* rulesets exist (githubPlan §7) |
| **PR #1 — Skeleton** | the working thin path above, `tests/fixtures/{constructs,api}/**`, first goldens | **this is your ruleset rehearsal** — githubPlan §7 step 3 wants a throwaway PR; make it this one instead |

This is a deliberate refinement of locked decision #7. It is still one hand-written foundation authored by you; splitting it means the contract is immutable before the skeleton's implementation pressure can quietly bend it.

---

## 3. The frozen interfaces

Everything in `packages/engine/src/bt5/core/` is owner-only and label-gated. These are the actual signatures.

### 3.1 Geometry — `core/types.py`

```python
Strand = Literal[1, -1]


class Topology(StrEnum):
    LINEAR = "linear"
    CIRCULAR = "circular"


@dataclass(frozen=True, slots=True)
class Interval:
    """Half-open [start, end), 0-based, in CONSTRUCT coordinates.

    On a circular construct, `end > construct.length` means the interval wraps
    the origin. There is exactly ONE representation of a wrapping interval;
    GenBank join() is normalised into it on import and re-split on export.
    """

    start: int
    end: int
    strand: Strand = 1


@dataclass(frozen=True, slots=True)
class Feature:
    interval: Interval
    kind: str  # GenBank feature key
    qualifiers: Mapping[str, tuple[str, ...]]  # ordered, multi-valued, byte-preserving
    uid: str


@dataclass(frozen=True)
class Construct:
    """The ONLY thing a rule is ever evaluated against.

    There is no API anywhere in BT5 that evaluates a rule against a bare `str`.
    That single decision is what makes 'evaluate on the assembled circular
    plasmid' the default rather than a feature someone has to remember.
    """

    sequence: str  # ACGT only, linearised at the GenBank origin
    topology: Topology
    features: tuple[Feature, ...]
    editable: tuple[Interval, ...]  # THE immutable-backbone model: the complement
    # of this is immutable, by construction
    codon_map: tuple[Interval, ...]  # one interval per codon, in translation order;
    # may wrap; may be discontiguous across a 2A
    annotations: Mapping[str, str]  # survives round-trip; SeqRecord drops these
    provenance: Provenance

    @property
    def length(self) -> int: ...
    def slice(self, iv: Interval) -> str: ...  # wrap- and strand-aware
    def tripled(self) -> tuple[str, int]: ...  # DnaChisel circular trick, offset L
    def is_editable(self, iv: Interval) -> bool: ...
```

`editable` and `codon_map` are the two fields that make R2 tractable, and they cannot be retrofitted. They are in the contract PR.

### 3.2 The rule protocol — `core/spec.py`

```python
class Severity(StrEnum):  HARD = "hard"; SOFT = "soft"; INFO = "info"
class Evidence(StrEnum):  A = "evidence_backed"; B = "contested"
                          C = "vendor_asserted"; D = "folklore"

@dataclass(frozen=True, slots=True)
class Breach:
    """One localized, attributable problem. The unit of everything downstream:
    the conflict panel, the infeasibility certificate, the 'which side is
    binding' per-window display, and the report all consume ONLY Breaches."""
    spec_id: str
    interval: Interval                 # construct coordinates
    magnitude: float                   # rule-native, > 0 == worse
    message: str                       # must name the exact offending substring
    slot_role: str | None              # "propagation" | "producer" | "target"
    detail: Mapping[str, float | str] = field(default_factory=dict)

@dataclass(frozen=True, slots=True)
class Evaluation:
    spec_id: str
    passes: bool                       # meaningful only when severity is HARD
    raw_score: float                   # SIGN CONVENTION: lower is always better
    breaches: tuple[Breach, ...]
    windows: tuple[tuple[Interval, float], ...] = ()   # per-window values
    n_evaluated: int = 0               # denominator for honest reporting

@dataclass(frozen=True)
class LatticeTerms:
    """A rule's OPT-IN into Tier-A. Returning None means Tier-B only.

    This is the interface that decouples lane 1 from lanes 3/4/5/6/7 completely:
    a rule can ship Tier-B-only on day 3 and upgrade itself to Tier-A on day 20
    WITHOUT the solver agent touching a single file."""
    forbidden: tuple[str, ...] = ()          # IUPAC; automaton-consumed; see §4.2
    codon_weights: Mapping[str, float] | None = None
    codon_pair_weights: Mapping[tuple[str, str], float] | None = None
    positional: Callable[[int, str], float] | None = None   # h(codon_index, codon)

class Spec(Protocol):
    id: ClassVar[str]
    version: ClassVar[str]
    evidence: ClassVar[Evidence]
    citation: ClassVar[str]                  # URL; CI asserts it starts with https://
    last_verified: ClassVar[str]             # ISO date; vendor rules go stale
    requires_construct: ClassVar[bool]       # False == cheap, evaluable on the bare ORF
    default_enabled: ClassVar[bool]          # folklore rules ship False
    severity: Severity
    param_schema: ClassVar[Mapping[str, object]]   # JSON Schema; the UI renders from this

    def gate(self, slot: ContextSlot) -> bool: ...
    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation: ...
    def localize(self, iv: Interval) -> Interval: ...
    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None: ...
```

Three properties of this protocol are load-bearing and each retires a specific risk:

- `Breach.interval` + `Breach.slot_role` → the conflict panel and the "which of the three contexts fired this" display are *free*, because every rule is forced to produce them. Retrofitting attribution later would touch every rule. **(R1)**
- `lattice_terms()` → solver ⊥ rules. **(R1, and the whole no-waiting property)**
- `requires_construct` → the null model evaluates 500 variants of cheap specs on the bare ORF and only the expensive specs on the full plasmid. Without this flag, normalization costs 500 circular plasmid folds and R6 becomes fatal. **(R6)**

### 3.3 Context and the three simultaneous constraint sets — `core/context.py`

```python
@dataclass(frozen=True)
class ContextSlot:
    role: Literal["propagation", "producer", "target"]
    host: HostId
    modality: Modality
    table_id: int  # NCBI genetic code. EXPLICIT. Never defaulted silently.
    # (critique Tier-1 #2: CTG=Ser in table 12 is a silent
    #  protein-changing bug)
    strand_of_interest: Strand  # which strand is the packaged genome / the mRNA
    enabled: bool = True


@dataclass(frozen=True)
class DesignContext:
    slots: tuple[ContextSlot, ...]  # 1..3, evaluated SIMULTANEOUSLY, never merged
    cassette_orientation: Strand  # required input; without it polyA/splice is backwards
    seed: int
    strict_biosecurity: bool
    engine_versions: Mapping[str, str]  # goes verbatim into Provenance
```

The three contexts are never collapsed into one rule set. Each `Spec` is instantiated once **per gating slot**, so a `Breach` always knows which context produced it, and two slots demanding opposite things produce two `Breach`es over the same `Interval` — which is exactly the input the conflict detector needs. Surfacing conflicts (brief §3.3) becomes a *consequence of the data model* rather than a feature somebody has to build.

### 3.4 Injected services — `core/services.py`

```python
class FoldEngine(Protocol):
    name: ClassVar[str]
    version: ClassVar[str]
    params: ClassVar[str]  # "rna_turner2004"; temperature; dangles

    def mfe(self, seq: str) -> float: ...
    def mfe_window(self, seq: str, iv: Interval) -> float: ...
    def accessibility(self, seq: str, iv: Interval, u: int) -> float | None: ...


class KmerIndex(Protocol):
    """Constructed ONLY from a Construct. There is deliberately no constructor
    that accepts an external sequence database — see §7 (biosecurity)."""

    @classmethod
    def of(cls, c: Construct, k: int) -> KmerIndex: ...
    def duplicates(self, min_len: int) -> Iterator[tuple[Interval, Interval]]: ...
    def revcomp_pairs(
        self, min_stem: int, max_loop: int
    ) -> Iterator[tuple[Interval, Interval]]: ...


class TableProvider(Protocol):
    def genetic_code(self, table_id: int) -> GeneticCode: ...
    def usage(self, host: HostId) -> CodonUsage: ...
    def weights(
        self, host: HostId, kind: Literal["cai", "tai", "stai", "csc"]
    ) -> Mapping[str, float]: ...


@dataclass(frozen=True)
class Services:
    fold: FoldEngine
    kmer: Callable[[Construct, int], KmerIndex]
    tables: TableProvider
    rng: np.random.Generator  # np.random.default_rng(seed), threaded everywhere
```

Rules receive `Services`; they never `import` lane 2 or lane 4. Injection is what makes lane 3 buildable on day 1 against the skeleton's seqfold and E. coli table, and *automatically* correct on day 30 against ViennaRNA and thirteen hosts, **with zero diff in lane 3**.

### 3.5 Results and failure — `core/result.py`

```python
@dataclass(frozen=True)
class ObjectiveScore:
    spec_id: str
    raw: float
    percentile: float  # [0,1] vs the empirical null. 1.0 == better than all nulls.
    null_n: int
    null_mean: float
    null_sd: float
    null_kind: Literal["host_frequency", "uniform_synonymous"]


@dataclass(frozen=True)
class Conflict:
    interval: Interval
    spec_ids: tuple[str, ...]
    kind: Literal["mutually_exclusive", "opposing_gradient", "immutable_region"]
    binding_spec_id: str  # which side is binding, per window
    relaxations: tuple[Relaxation, ...]


@dataclass(frozen=True)
class Relaxation:
    spec_id: str
    change: str  # "raise gc_max 0.60 -> 0.64"
    predicted_cost: Mapping[str, float]  # delta percentile on EVERY other objective


@dataclass(frozen=True)
class InfeasibilityCertificate:
    """A zero-variant merged segment in the mutation space is a PROOF, not a guess."""

    interval: Interval
    protein_span: tuple[int, int]
    minimal_conflicting_specs: tuple[str, ...]
    proof: Literal["empty_mutation_space", "automaton_dead_state", "immutable_region"]
    relaxations: tuple[Relaxation, ...]


class InfeasibleConstraints(Exception):
    certificate: InfeasibilityCertificate


@dataclass(frozen=True)
class Candidate:
    label: str  # "balanced", "min-structure", "native", "host-frequency"
    construct: Construct
    orf: str
    scorecard: ScoreCard
    design_hash: str  # content hash; goes on the tube label
    codon_distance_to: Mapping[str, float]


@dataclass(frozen=True)
class DesignResult:
    candidates: tuple[Candidate, ...]
    native_baseline: Candidate | None  # "don't optimize" as a FIRST-CLASS output
    conflicts: tuple[Conflict, ...]
    biosecurity: BiosecurityReport
    provenance: Provenance
```

`native_baseline` being a field of the result type — not a UI afterthought — is how critique Tier-3 #9 ("use the native sequence" as an honest recommendation) becomes structurally impossible to forget.

### 3.6 The engine entry point — `core/api.py`

```python
def design(req: DesignRequest) -> DesignResult: ...
def design_stream(req: DesignRequest) -> Iterator[ProgressEvent | DesignResult]: ...
```

### 3.7 The oracle — a deliberate deviation from githubPlan §6

githubPlan's `verify_solution(protein, cons, dna, *, table_id)` takes a bare string. That signature cannot express the thing this product exists to guarantee. PR #0 ships:

```python
def verify_construct(plan: CassettePlan, ctx: DesignContext, c: Construct) -> None:
    """Independently re-derives everything. Raises VerificationError.

    I1  alphabet          I5  stops (interior forbidden, terminal per declaration)
    I2  frame             I6  forbidden motifs on the CIRCULAR construct, pattern set
    I3  round trip            closed under revcomp, including junction- and
    I4  initiator             origin-spanning hits, including inside immutable backbone
    I7  GC band, global + windowed, windows wrapping the origin
    I8  homopolymer / repeat ceiling across the whole construct
    I9  every immutable interval is byte-identical to the input backbone   <-- NEW
    I10 cassette frame invariant (brief A4) across the assembled cassette  <-- NEW
    """


def verify_solution(protein: str, cons: Constraints, dna: str, *, table_id: int = 1) -> None:
    """Thin linear wrapper. Kept so githubPlan's invariant suite and the
    `python-quality` grep for `_verify: bool = True` are unchanged."""
```

**I9 is the highest-value new invariant in the system.** It makes "the optimizer silently edited the user's backbone" — the worst possible bug in a vector-context tool — a raised exception rather than a shipped plasmid.

---

## 4. Module decomposition

### 4.1 Directory tree and ownership

```
BT5/
├── CLAUDE.md                              OWNER
├── pyproject.toml, uv.lock                OWNER  — every dep PRE-DECLARED in PR #0
├── pnpm-workspace.yaml, pnpm-lock.yaml    OWNER
├── data/                                  OWNER, label-gated approved:data-change
│   ├── genetic_codes/                       NCBI tables — can silently change the PROTEIN
│   ├── codon_usage/                         per-host usage + PaxDb HEG reference sets
│   ├── plasmids/                            real Addgene GenBank goldens
│   └── MANIFEST.sha256
├── benchmarks/                            OWNER, label-gated — EVERY metric pre-registered
├── tests/
│   ├── invariants/  goldens/  data_integrity/     OWNER, label-gated
│   ├── fixtures/constructs/*.json         OWNER writes v1 — unblocks lanes 6, 8
│   ├── fixtures/api/*.json                OWNER writes v1 — unblocks lane 10
│   └── e2e/                               LANE 10
└── packages/
    ├── engine/src/bt5/
    │   ├── core/          OWNER, FROZEN     types, spec, context, services, registry, result
    │   ├── verify.py      OWNER, gated      the oracle
    │   ├── pipeline.py    OWNER             fixed stage order; delegates via protocols only
    │   ├── solver/        LANE 1  ─ mutation space, Aho–Corasick, Tier-A DP, β-sampler,
    │   │                             Tier-B localized repair, infeasibility certificates
    │   ├── codon/         LANE 2  ─ genetic code tables, host usage, CAI/tAI/stAI/%MinMax/
    │   │                             CFD/ENc/CSC, harmonization, user-supplied tables
    │   ├── rules/         LANE 3  ─ patterns/ : restriction, recombinase, polyA hexamers,
    │   │                             non-B DNA, ARE, Pol III, Chi, CpG×3, out-of-frame stops
    │   │                    models/ : MaxEntScan splice, Salis TIR/RBS, promoter calculator,
    │   │                             Kozak/uAUG, polyA downstream-element escalation
    │   ├── structure/     LANE 4  ─ FoldEngine impls (seqfold, ViennaRNA), windowed ΔG with
    │   │                             incremental invalidation, per-engine threshold table,
    │   │                             5′ terms, IVT terms, DegScore
    │   ├── vector/        LANE 5  ─ GenBank/SnapGene/FASTA I/O, IntervalRemapper, insertion-
    │   │                             site detection, Golden-Gate cassette auto-detect,
    │   │                             construct assembly, editable-region computation, KmerIndex
    │   ├── assembly/      LANE 6  ─ vendor profiles + complexity precheck, repeats/homopolymer/
    │   │                             GC-extent manufacturability, overhang design (tatapov),
    │   │                             domestication, primers, fragment splitting, order files
    │   ├── cassette/      LANE 7  ─ protein validator, genetic-code selection, tag/linker/2A/
    │   │                             protease library, repetitive-protein analysis, protein
    │   │                             liability scan, signal peptide, BIOSECURITY screening
    │   └── score/         LANE 8  ─ null model, percentile normalization, weight vectors,
    │                                 presets + provenance, λ-sweep, diversity-filtered gallery,
    │                                 conflict detection, QC report, design hash, provenance
    ├── engine/tests/{solver,codon,rules,structure,vector,assembly,cassette,score}/   lane-owned
    ├── server/            LANE 9  ─ FastAPI on 127.0.0.1, ProcessPoolExecutor jobs, SSE
    └── ../apps/web/       LANE 10 ─ React/TS, sequence viewer, sliders, gallery, conflict panel
```

Ten lanes. Six run concurrently in steady state (see §6 on why not ten).

### 4.2 Two decisions inside lane 1 worth stating now

**Never scan two strands.** The forbidden-pattern set is closed under reverse complement at construction time, and only the forward strand is scanned. A minus-strand hit *is* a plus-strand hit of the revcomp. This deletes githubPlan's #1 named agent bug class (`minus_strand_bamhi`) at the type level rather than testing for it.

**Windowed GC in the DP — the honest answer.** Exact sliding-window GC needs ~17 codons of history; naive augmentation is exponential and there is no literature precedent (brief §7-tech-i). Three options were evaluated:

| Option | Guarantee | Cost |
|---|---|---|
| Lagrangian multiplier on GC deviation as a per-codon linear term | none (steering only) | free |
| Block-state DP: state = (automaton, GC count in current 17-codon block ∈ 0..51) | *block-wise* GC, a strictly weaker surrogate for sliding-window | 600 × 52 = 31.2 k states → ~92 M relaxations → **0.5–1.5 s** vectorized |
| Full 17-codon history | exact | intractable |

**Recommendation: Lagrangian term in Tier A (steering) + Tier-B repair for the hard bound + independent validator, with the block-state DP available behind an opt-in "guaranteed GC window" flag.** This means **windowed GC is the one hard constraint not guaranteed by construction** — permitted by the brief's "or by reject-and-repair plus an INDEPENDENT final validator," but it is a decision the owner should ratify in week one rather than discover in week five. G1 measures the block-state cost so the decision is made on a number.

---

## 5. Why nobody waits — the six mechanisms

The hard constraint is: no agent blocked on another, no agent editing another's files. Six mechanisms, each removing one class of coupling.

**M1 — Registry autodiscovery removes the last shared file.** `core/registry.py` walks `bt5.rules.*`, `bt5.assembly.*`, `bt5.structure.*` with `pkgutil.walk_packages` and collects `@register`-decorated `Spec` classes. Adding a rule edits **zero** shared files — not even an `__init__.py`. This is the single highest-leverage anti-conflict decision in the design; the obvious alternative (a hand-maintained `SPECS = [...]` list) would put every rule lane into the same file on every PR.

**M2 — `lattice_terms()` decouples the solver from every rule.** Lane 1 never imports a rule. Lane 3 never imports the solver. A rule upgrading from Tier-B to Tier-A is a lane-3-only diff.

**M3 — Services injection decouples rules from tables, folding and k-mer indexing.** Lane 3 codes against `svc.fold` on day 1, backed by the skeleton's seqfold; lane 4 replaces the implementation on day 25; lane 3's diff is empty.

**M4 — The skeleton already implements every layer.** A lane that conceptually consumes another's output consumes the skeleton's degenerate version today and the real one after that lane merges — *without an import change*, because the interface is identical. Lane 6 (assembly) needs assembled circular constructs from lane 5; it gets six of them from `tests/fixtures/constructs/*.json` on day 1. Lane 10 (web) needs API responses; it gets them from `tests/fixtures/api/*.json` on day 1.

**M5 — The UI is manifest-driven, not rule-aware.** `GET /capabilities` returns the rule catalogue built from the registry:

```ts
export interface RuleDescriptor {
  id: string; title: string;
  severity: 'hard' | 'soft' | 'info';
  evidence: 'evidence_backed' | 'contested' | 'vendor_asserted' | 'folklore';
  defaultEnabled: boolean;
  unit: string;
  appliesTo: { hosts: string[]; modalities: string[]; roles: string[] };
  params: JSONSchema7;      // the slider/toggle panel renders FROM this
  citation: string;         // rendered as the evidence badge's link
  lastVerified: string;     // stale vendor rules render with a warning chip
}

export interface Capabilities {
  engineVersion: string;
  hosts: HostDescriptor[];
  rules: RuleDescriptor[];
  presets: Preset[];
  foldEngine: { name: string; version: string; params: string };
  biosecurity: { available: boolean; databaseVersion: string | null };
}
```

Lane 10 hard-codes no rule name. Lanes 3–7 add rules without touching `apps/web/`. Evidence badges (brief §6.8) arrive for free because `evidence` and `citation` are `ClassVar`s the protocol forces every rule to declare.

**M6 — Everything shared is pre-registered in PR #0, not added incrementally.** This is the mechanism that would be most easily overlooked, and githubPlan's label gates make its absence fatal:

| Shared thing | Why lanes would otherwise serialize on the owner | Pre-registration |
|---|---|---|
| `pyproject.toml` / `uv.lock` | every dep addition → conflict + the §8.2 lockfile storm | **all deps declared in PR #0**: biopython, numpy, pyahocorasick, seqfold, primer3-py, pydantic, fastapi, tatapov, pyyaml, `viennarna` and `commec` as optional extras |
| `benchmarks/tolerances.yaml` | label-gated; a new metric needs an owner PR | **every planned metric name pre-registered**, `blocking: false` until its lane lands; flipping to blocking is an owner action |
| `benchmarks/panel.json` | label-gated | full panel up front: GFP, mCherry, Cas9, trastuzumab HC, p53, GC-rich human, FVIII, scFv, plus the adversarial set |
| `data/**` | label-gated | **split the line**: data that can silently change the *protein* (genetic codes, codon usage, HEG sets) stays in root `data/` and stays gated. Data that changes *which liabilities are flagged* (motif tables, vendor profiles, MaxEnt matrices) lives **in-lane** at `bt5/<lane>/data/` with a mandatory `_provenance.json` sidecar (source URL, license, retrieval date, sha256) and a lane-local integrity test. Same discipline, decentralized. |
| `tests/invariants/**` | label-gated | universal properties pre-registered; lanes own `packages/engine/tests/<lane>/` including their own Hypothesis properties |
| generated TS API types | server lane would have to write into `apps/web/` | types are **generated** into `apps/web/src/api/generated/` by a script the web lane owns; a CI job fails if stale. Lane 9 never touches `apps/web/`. |

Lane-local snapshots are safe to be lane-owned because `goldens-not-hand-edited` deletes and regenerates *all* snapshots and diffs — it catches hand-editing regardless of who owns the file. Extend `path_guard.RULES` to `**/__snapshots__/**` only for `tests/goldens/`.

### 5.1 Dependency graph

```
                      ┌──────────────────────────────────────┐
                      │  core/  (OWNER, frozen after wave 1) │
                      │  types · spec · context · services   │
                      │  registry · result · api · schema    │
                      └──────────────────────────────────────┘
                          ▲   ▲   ▲   ▲   ▲   ▲   ▲   ▲
     ┌────────────────────┘   │   │   │   │   │   │   └────────────────┐
     │        ┌───────────────┘   │   │   │   │   └───────┐            │
  ┌──┴───┐ ┌──┴────┐ ┌────┴──┐ ┌──┴───┴──┐ ┌─┴──────┐ ┌──┴──────┐ ┌───┴────┐
  │ L1   │ │ L2    │ │ L3    │ │ L4      │ │ L5     │ │ L6      │ │ L7     │
  │solver│ │codon  │ │rules  │ │structure│ │vector  │ │assembly │ │cassette│
  └──────┘ └───────┘ └───────┘ └─────────┘ └────────┘ └─────────┘ └────────┘
       ▲        ▲         ▲          ▲          ▲          ▲          ▲
       └────────┴─────────┴──────────┴──────────┴──────────┴──────────┘
                                     │  (registry + Services only — NO imports)
                              ┌──────┴──────┐
                              │ L8  score   │   consumes Evaluation/ObjectiveScore
                              └──────┬──────┘   generically; names no spec
                                     │
                              ┌──────┴──────┐
                              │ L9  server  │   consumes core/api + core/schema
                              └──────┬──────┘
                                     │  OpenAPI → generated TS
                              ┌──────┴──────┐
                              │ L10 web     │   consumes /capabilities + fixtures
                              └─────────────┘
```

Every arrow points at `core/` or at a generated artifact. **There is not one arrow between two lanes.** L8, L9 and L10 look sequential but are not: L8 works against the skeleton's three registered rules, L9 against `core/api` which exists in PR #0, L10 against recorded fixtures. All three can start on day 8.

---

## 6. Build sequence

### Wave 0 — days 1–7, owner only

PR #0 (Contract) then PR #1 (Skeleton). Then the seven gates.

### Week-one go/no-go gates

These are the point of Framing C. Each is a number, measured against the skeleton, before any agent starts.

| Gate | Measurement | Pass | Fail → |
|---|---|---|---|
| **G0** (R1) | Implement the three skeleton rules against `Spec` and confirm each needed **zero** additions to `Evaluation`/`Breach` | 0 additions | the protocol is wrong; iterate before wave 1 |
| **G1** (R3) | Tier-A DP, 1000 codons, ~600 automaton states, forward + revcomp closure, vectorized | ≤ 50 ms; block-state variant ≤ 2 s; `verify_construct` finds 0 motifs | > 500 ms plain → drop Tier-A, re-plan around repair-only and say so in the README |
| **G2** (R2, R8) | Real Addgene lentiviral GenBank round-trips byte-identically incl. an origin-spanning `join()`; a BsaI site planted **across the insert/backbone junction** and one **across the origin** are both caught on the circular construct; I9 catches a deliberate backbone edit | all three | the coordinate representation is wrong — this is the most expensive thing to discover late |
| **G3** (R6) | 200-variant null, 500 aa protein, cheap specs on bare ORF | < 2 s; optimized design > 0.90 percentile on its own objective, ≈ 0.50 on an ignored one | > 10 s → normalization leaves the interactive path and becomes a report-time-only feature |
| **G4** (R4) | λ-sweep of 24 → greedy max-min pick of 5 | pairwise codon distance ≥ 15% | < 5% → the gallery premise is dead; switch to ε-constraint enumeration + explicit diversity constraints **before** lane 10 builds a UI around it |
| **G5** (R5) | `junction_trap` adversarial fixture | raises `InfeasibleConstraints` with ≤ 3 spec ids and an interval ≤ 10 codons | "no solution" only → redesign the mutation-space merge before rule lanes multiply the constraint count |
| **G6** (R7) | seqfold vs ViennaRNA ΔG on −4…+37 over 100 sequences | offset is ~constant, and the per-engine threshold table demonstrably switches | non-constant → every structure threshold needs independent per-engine calibration; scope lane 4 up accordingly |
| **G7** | End-to-end skeleton, 500 aa | ≤ 10 s wall | budget re-allocation before rules multiply the cost |

**G2 and G4 are the two that should scare you most.** G2 because a wrong coordinate model is the most expensive possible late discovery; G4 because it is the one gate whose failure invalidates a *product* decision rather than a technical one, and the only one nobody had on their list.

### Wave 1 — days 8–20, three agents

Only the three lanes that can still force a `core/` change:

- **L1 solver** — R3 (windowed GC decision), R5 (infeasibility certificates)
- **L5 vector** — R2 (circular + immutable + insertion detection), R8 (round-trip)
- **L8 score** — R6 (normalization at scale), R4 (gallery diversity), default weight vector with written provenance

Three agents, not eight, because a `core/` change during wave 1 costs three rebases. Wave 1 has a documented escalation path: a `core/` change is an owner-authored PR, batched at most once per week, announced to all three sessions.

**`core/` freezes for real when wave 1 lands.** This is the honest version of "interfaces frozen in the foundation PR": frozen from the moment agent #4 starts, deliberately provisional for the three agents whose job is to break it. That is precisely "discover fatal design errors before 8 agents have built on top of them."

### Wave 2 — days 18–50, remaining seven lanes

L2 codon · L3 rules · L4 structure · L6 assembly · L7 cassette · L9 server · L10 web — all pure consumers of a frozen core.

### Concurrency reality

githubPlan §8.1: a free personal account gets 20 concurrent jobs; a Python PR consumes ~12. Eight open PRs demand ~96 slots. **Run 6 concurrent open PRs maximum**, ten lanes rotating, with githubPlan's mitigations applied in order: `cancel-in-progress` (already in), in-job path filters (already in), draft PRs by default, and — first thing to cut — move CodeQL to `push: [main] + schedule` and drop `codeql-passed` from required checks. This is a CI-capacity limit, not an architecture limit: the decomposition supports ten simultaneously.

### How each lane thickens the same skeleton path

| Lane | The skeleton stage it thickens | The skeleton keeps working because… |
|---|---|---|
| L1 solver | Tier-A DP → adds codon-pair state, β-sampler, block-state GC, real Tier-B repair | the `design()` contract is unchanged; only the quality of the sequence improves |
| L2 codon | one host → thirteen hosts, all genetic code tables, harmonization | `TableProvider` protocol unchanged |
| L3 rules | three rules → ~40 registered rules | autodiscovery; the pipeline never enumerates rules |
| L4 structure | seqfold → ViennaRNA + windowed incremental + IVT terms | `FoldEngine` protocol unchanged |
| L5 vector | one hand-specified insertion point → annotation-driven detection, Golden Gate auto-detect, SnapGene import | `Construct` shape unchanged |
| L6 assembly | GenBank export → overhangs, domestication, primers, fragments, vendor order files | additive output fields |
| L7 cassette | bare protein → tags, linkers, 2A, signal peptides, repetitive-protein mode, biosecurity | `CassettePlan` shape unchanged |
| L8 score | one preset → presets, ε-constraint mode, conflict panel, protocol output, calibration view | `ScoreCard` shape unchanged |
| L9 server | three endpoints → jobs, cancellation, streaming, file handling | OpenAPI generated from frozen `core/schema` |
| L10 web | one page → viewer, sliders, gallery, conflict panel, order-file download | manifest-driven |

At no point does the end-to-end path stop working. That is the whole plan.

---

## 7. The four Tier-1 critique items, structurally

These are mandatory, so each gets a *structural* home rather than a task:

1. **Biosecurity.** `bt5/cassette/biosecurity/` (lane 7). `pipeline.py` calls `screen_protein()` **before** the mutation space is built — protein-level screening is the one layer BT5's own output cannot defeat. `commec` is an optional extra with a clear "screening unavailable" banner, never a silent skip. Crucially: **`KmerIndex.of()` takes a `Construct` and nothing else.** There is no constructor accepting an external database and no "minimize identity to reference" objective anywhere in the type system, with a CI grep asserting it. The policy is enforced by the interface, not by a rule an agent might not read.
2. **Genetic code table.** `ContextSlot.table_id` is required with no default. `verify_construct` re-derives translation from Biopython's NCBI tables independently of the optimizer's own table. Property test `translate(back_translate(p, t), t) == p` for every shipped table, in `tests/invariants/` (owner-gated). `/transl_table` written to output GenBank.
3. **Objective normalization.** `ObjectiveScore.percentile` is the *only* quantity the weighted sum operates on, so sliders are linear and unit-free by construction. `requires_construct` keeps the cost affordable. Null cached in local SQLite keyed by `(protein sha256, host, table_id, spec-set version, N, seed)` so slider moves never recompute. Default weight vectors are checked-in preset JSON with a mandatory `rationale` field; CI asserts every nonzero weight carries a citation.
4. **Repetitive proteins / antibodies.** Repeat detection runs on the **input protein**, in lane 7, *before any codon is chosen*, and emits `editable`-region metadata that forces divergent synonymous assignment across copies — a constraint that overrides codon optimality, not an objective competing with it. The scFv/CAR and trastuzumab heavy chain are in `benchmarks/panel.json` from PR #0, so the property is gated from day one rather than added as a mode later.

---

## 8. What I would watch for, and the honest weak points

- **The block-state GC decision (G1) is the one place I am recommending a weaker guarantee than the brief's ideal.** Windowed GC will be enforced by repair + independent validator, not by construction. Ratify it in week one.
- **G4 has no fallback that is cheap.** If the gallery collapses, ε-constraint enumeration with explicit diversity constraints is a substantially larger lane-8 scope. Measure it in week one, not week five.
- **Lane 3 is oversized** (~40 rules across `patterns/` and `models/`). It is the designated split point: `bt5/rules/patterns/` and `bt5/rules/models/` are already separate directories precisely so a second agent can take `models/` with no restructuring.
- **`core/` is provisionally frozen for three agents for two weeks.** That is a real cost, deliberately paid, and it is the entire point of the framing.
- **The foundation PR is large** — realistically 4–6k hand-written lines. That is the price of the walking skeleton, and it is cheaper than eight agents building on a wrong `Breach`.