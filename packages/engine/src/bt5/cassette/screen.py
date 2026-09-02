"""Protein-level biosecurity screening behind the `Screen` protocol.

BT5's core function is producing a functionally identical sequence with
maximally different nucleotides -- the textbook method for evading
nucleotide-homology screening. Protein-level screening is therefore the one
layer BT5's own output cannot defeat, and it is the only place BT5 flags a
hazard before an order goes out. `core/context.BiosecurityVerdict` states that
threat model; this module makes the screen real without letting it lie.

**The detection is not ours.** `commec` (the IBBIS Common Mechanism) does the
biology: HMM biorisk profiles, BLASTX/BLASTN best-match against regulated
pathogens, and a low-concern clear step. This module is the plumbing around it
-- a protocol, a concrete `commec`-backed implementation, and an honest
degradation path -- exactly the way `structure/vienna.py` wraps ViennaRNA
behind `FoldEngine`.

**The one failure this module exists to prevent is a screen that reports
"clear" when it never ran.** Every other bug here is recoverable; that one
hands a user a false assurance about a hazard, and it is the bug that passes
every mechanical check -- the types are right, the field says "clear". So the
degradation is fail-safe by construction:

  * `commec` absent, no database configured, a timeout, a non-zero exit, an
    unparseable result, or an outcome word we do not recognise -> `not_run`.
  * Only an explicit `commec` "Clear" outcome becomes `status="clear"`, and
    only when a real screen actually ran against a named database.

`not_run` is never upgraded to `clear` by inference. `_status_for` defaults
every unrecognised outcome to `not_run`, so even a future `commec` that adds a
new outcome word fails safe rather than reading as clean.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Protocol

from bt5.core.context import BiosecurityVerdict

#: The console script the `screen` extra installs (`commec>=0.2`). Availability
#: is "can we run the CLI", checked with `shutil.which`, because the CLI -- not
#: an importable Python API -- is commec's documented, version-stable surface.
COMMEC_BINARY = "commec"

#: A screen that overruns this is `not_run`, never `clear`: an unfinished screen
#: is an unknown result, and an unknown result is exactly what must not read as
#: clean. Generous because commec runs BLAST against large databases; the point
#: is to bound a hang, not to race a real screen.
DEFAULT_TIMEOUT_S = 900.0

#: commec's per-query outcome vocabulary is Clear / Warning / Flag. The mapping
#: to BT5's status is deliberately conservative and is the one decision here
#: that changes what the app REFUSES to build, so it is stated in one place:
#:
#:   commec "Clear"   -> "clear"  (screen ran, nothing matched)
#:   commec "Warning" -> "flag"   (surfaced for review; still emits)
#:   commec "Flag"    -> "block"  (a regulated-pathogen match; refuses to emit)
#:
#: Anything else -- an empty string, an error marker, a word a newer commec
#: introduces -- is absent from this table and `_status_for` returns "not_run".
#: That default is the safety property: a result we cannot interpret is not a
#: clean one.
_STATUS_BY_OUTCOME: dict[str, str] = {
    "clear": "clear",
    "warning": "flag",
    "warn": "flag",
    "flag": "block",
}


class ScreenUnavailableError(RuntimeError):
    """A biosecurity screen could not be run.

    Distinct from a screen that ran and returned a verdict: this is the absence
    of a screen. Callers do not let it become a "clear" -- `load_screen`
    returns a `NullScreen` reporting `not_run` rather than raising, so the
    report cannot imply a clean result that was never obtained.
    """


class BiosecurityBlockedError(RuntimeError):
    """Raised to refuse emission of a blocked design.

    `BiosecurityVerdict.may_proceed` is `status != "block"`. A "block" verdict
    means commec matched a regulated pathogen; the design must not be emitted.
    `guard_emission` raises this so the refusal is an exception a caller cannot
    forget to check, not a boolean it can drop on the floor.
    """

    def __init__(self, verdict: BiosecurityVerdict) -> None:
        self.verdict = verdict
        super().__init__(
            f"biosecurity screen returned status={verdict.status!r}; this design "
            f"must not be emitted. {verdict.detail}".rstrip()
        )


@dataclass(frozen=True, slots=True)
class CommecOutcome:
    """The two facts BT5 needs out of a commec run, normalised.

    The seam between "invoke commec" and "decide a status". Isolating it keeps
    the status mapping -- the part with the correctness risk -- pure and
    testable without commec installed, and quarantines the commec-specific
    argv-and-JSON handling that CI never exercises (the `screen` extra is not in
    the bootstrap install) into one clearly-marked function.
    """

    #: commec's per-query outcome word, e.g. "Clear" / "Warning" / "Flag".
    outcome: str
    #: The reference-database version the screen ran against, when commec
    #: reports one. Recorded on the verdict so a result traces to what produced
    #: it; None when commec did not name it.
    database_version: str | None = None


class Screen(Protocol):
    """Protein-level biosecurity screen. `commec` is the bundled default.

    Mirrors `FoldEngine`: a protocol with a real implementation and a
    degradation. `screen` takes the INPUT protein, before any codon is chosen,
    and returns a `BiosecurityVerdict` whose `status` a caller renders verbatim
    and never upgrades.
    """

    name: ClassVar[str]

    def screen(self, protein: str) -> BiosecurityVerdict:
        """Screen `protein` and return a verdict. Never raises for a hazard.

        A detected hazard is a `status="block"` verdict, not an exception; a
        screen that could not run is `not_run`. The only inputs that raise are
        programmer errors (an empty protein), never a screening outcome.
        """
        ...


def _status_for(outcome: str) -> str:
    """Map a commec outcome word to a BT5 status, failing safe.

    Unrecognised -> "not_run". This is the single line that stops an
    uninterpretable result from reading as clean.
    """
    return _STATUS_BY_OUTCOME.get(outcome.strip().lower(), "not_run")


@dataclass(frozen=True, slots=True)
class NullScreen:
    """The screen that admits it did not run.

    Returned by `load_screen` whenever a real screen cannot be performed --
    commec absent, or no database configured. Its verdict is always `not_run`
    with `database_version=None`, so the report says exactly that and no reader
    can mistake it for a clean result.
    """

    reason: str = "protein-level biosecurity screening did not run"

    name: ClassVar[str] = "null"

    def screen(self, protein: str) -> BiosecurityVerdict:
        if not protein:
            raise ValueError("cannot screen an empty protein")
        return BiosecurityVerdict("not_run", None, self.reason)


@dataclass(frozen=True, slots=True)
class CommecScreen:
    """`Screen` over the commec CLI.

    The database directory is required: commec cannot screen without its
    reference databases, and a screen without them is not a screen. The commec
    invocation is a seam (`runner`) so the status-mapping and degradation logic
    is testable without commec installed; the default runner shells out to the
    `commec` binary and is the part CI does not reach.
    """

    database_dir: Path
    timeout_s: float = DEFAULT_TIMEOUT_S
    #: Injectable for tests; defaults to the real commec CLI binding.
    runner: Callable[[str, Path, float], CommecOutcome] = field(
        default_factory=lambda: _run_commec_cli
    )

    name: ClassVar[str] = "commec"

    def screen(self, protein: str) -> BiosecurityVerdict:
        if not protein:
            raise ValueError("cannot screen an empty protein")
        try:
            outcome = self.runner(protein, self.database_dir, self.timeout_s)
        except subprocess.TimeoutExpired:
            return BiosecurityVerdict(
                "not_run",
                None,
                f"commec screen exceeded {self.timeout_s:g}s and was stopped; an "
                f"unfinished screen is an unknown result, not a clean one",
            )
        except ScreenUnavailableError as exc:
            return BiosecurityVerdict("not_run", None, str(exc))
        except (OSError, subprocess.SubprocessError, ValueError) as exc:
            return BiosecurityVerdict(
                "not_run",
                None,
                f"commec screen did not produce a usable result ({exc}); reported as "
                f"not run rather than clear",
            )

        status = _status_for(outcome.outcome)
        if status == "not_run":
            return BiosecurityVerdict(
                "not_run",
                outcome.database_version,
                f"commec returned an outcome BT5 does not recognise ({outcome.outcome!r}); "
                f"reported as not run rather than clear",
            )
        return BiosecurityVerdict(
            status,  # type: ignore[arg-type]  # narrowed to a literal member above
            outcome.database_version,
            f"commec screen outcome: {outcome.outcome}",
        )


def commec_available() -> bool:
    """True when the commec CLI can actually be invoked."""
    return shutil.which(COMMEC_BINARY) is not None


def load_screen(
    *, database_dir: str | Path | None = None, timeout_s: float = DEFAULT_TIMEOUT_S
) -> Screen:
    """The best screen available, always a usable `Screen`.

    Unlike `load_fold_engine`, this never returns None: a screen must always
    produce a verdict, and the honest verdict when no real screen can run is
    `not_run`. So a `CommecScreen` when commec is installed AND a database is
    configured; otherwise a `NullScreen` carrying the specific reason.
    """
    reason = _degradation_reason(database_dir)
    if reason is not None:
        return NullScreen(reason=reason)
    assert database_dir is not None  # _degradation_reason guarantees this
    return CommecScreen(database_dir=Path(database_dir), timeout_s=timeout_s)


def screen_degradation_reason(*, database_dir: str | Path | None = None) -> str | None:
    """What to record in `Provenance.degradations`, or None when a real screen ran.

    Mirrors `structure.vienna.degradation_reason`: the caller records this so a
    missing screen degrades visibly instead of vanishing.
    """
    return _degradation_reason(database_dir)


def _degradation_reason(database_dir: str | Path | None) -> str | None:
    if not commec_available():
        return (
            "commec is not installed, so protein-level biosecurity screening did not "
            'run; install the screen extra with `uv pip install -e ".[screen]"`. The '
            "verdict is reported as not_run, never clear"
        )
    if database_dir is None:
        return (
            "commec is installed but no reference database is configured, so "
            "protein-level biosecurity screening did not run; the verdict is reported "
            "as not_run, never clear"
        )
    return None


def guard_emission(verdict: BiosecurityVerdict) -> None:
    """Refuse to emit a blocked design.

    `raise`s `BiosecurityBlockedError` when `not verdict.may_proceed` (i.e.
    `status == "block"`). A "clear", "flag" or "not_run" verdict passes: this
    guard enforces only the one refusal the frozen type defines, and never
    weakens it -- there is no argument that turns a block into a pass. Whether a
    "not_run" should also fail-closed under `DesignContext.strict_biosecurity`
    is an emit-policy decision owned by the design runner, not this guard.
    """
    if not verdict.may_proceed:
        raise BiosecurityBlockedError(verdict)


def _run_commec_cli(protein: str, database_dir: Path, timeout_s: float) -> CommecOutcome:
    """Invoke the commec CLI on a single protein and read its JSON result.

    NOT EXERCISED BY CI: the bootstrap install is `dev,fold,export`, so the
    `screen` extra -- and therefore commec -- is absent, and `load_screen`
    returns a `NullScreen` before this is ever constructed. It is written to be
    correct and defensive rather than clever: anything unexpected raises, and
    `CommecScreen.screen` turns every raise into `not_run`. commec's exact
    JSON schema is a version-dependent detail this deliberately does not hard-
    code beyond the two fields BT5 needs; an outcome it cannot find is reported
    as not run, not clear.
    """
    if shutil.which(COMMEC_BINARY) is None:
        raise ScreenUnavailableError("commec binary is not on PATH")

    with tempfile.TemporaryDirectory(prefix="bt5-commec-") as tmp:
        workdir = Path(tmp)
        fasta = workdir / "query.faa"
        fasta.write_text(f">bt5_query\n{protein}\n")
        out_prefix = workdir / "query"
        subprocess.run(
            [
                COMMEC_BINARY,
                "screen",
                "--database-dir",
                str(database_dir),
                "--output",
                str(out_prefix),
                str(fasta),
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        return _parse_commec_output(out_prefix, database_dir)


def _parse_commec_output(out_prefix: Path, database_dir: Path) -> CommecOutcome:
    """Extract (outcome, database_version) from commec's JSON, defensively.

    Reads the first per-query recommendation it finds under a small set of the
    keys commec is known to use, and the database version if commec records one.
    A shape it does not recognise raises `ValueError`, which becomes `not_run`.
    """
    candidates = [out_prefix.with_suffix(".json"), Path(f"{out_prefix}.output.json")]
    payload_path = next((p for p in candidates if p.exists()), None)
    if payload_path is None:
        raise ValueError("commec produced no JSON output to read a recommendation from")

    data = json.loads(payload_path.read_text())
    outcome = _first_value(data, ("recommendation", "outcome", "screen_status", "status"))
    if outcome is None:
        raise ValueError("commec JSON carried no recognisable recommendation field")
    version = _first_value(data, ("database_version", "db_version", "version"))
    if version is None:
        version = _read_database_version(database_dir)
    return CommecOutcome(outcome=str(outcome), database_version=version)


def _first_value(data: object, keys: tuple[str, ...]) -> str | None:
    """Depth-first search for the first of `keys` holding a scalar, or None."""
    stack: list[object] = [data]
    while stack:
        node = stack.pop()
        if isinstance(node, dict):
            for key in keys:
                value = node.get(key)
                if isinstance(value, (str, int, float)) and not isinstance(value, bool):
                    return str(value)
            stack.extend(node.values())
        elif isinstance(node, list):
            stack.extend(node)
    return None


def _read_database_version(database_dir: Path) -> str | None:
    """A version stamp beside the reference databases, if one is written there."""
    for name in ("VERSION", "version.txt", "database_version.txt"):
        stamp = database_dir / name
        if stamp.exists():
            text = stamp.read_text().strip()
            if text:
                return text
    return None
