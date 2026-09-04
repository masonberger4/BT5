"""Every text read, write and subprocess decode in this repo must name its encoding.

Python resolves a missing `encoding=` through `locale.getencoding()`. On Linux CI
that is UTF-8, so an omission is invisible here; on a Windows checkout it is cp1252
and the same command against the same commit behaves differently. Three ways:

  * `docs/` holds characters cp1252 cannot decode at all (U+221D, U+207B, U+2510,
    U+03C1, U+FE0F), so a read raises `UnicodeDecodeError` -- that is how #102 was
    found, via `test_every_weighted_ref_is_a_real_row_in_the_brief`;
  * far worse, cp1252 *maps* U+2014, U+2265 and U+00D7 to different characters and
    does not raise. Several gates read engine source as text and pattern-match it
    (`test_no_expression_claims`, `test_oracle_independence`, `surface.py`), and a
    provenance assertion evaluated against mojibake is not reading the file the
    reviewer read;
  * the exports -- the GenBank at `cli.py` and the vendor order CSV at
    `score/order.py` -- would be emitted in the platform encoding. For a tool whose
    product is an auditable design record, platform-dependent output bytes are a
    defect whether or not anyone has hit one.

WHY THIS IS A TEST AND NOT A LINT RULE. Ruff ships the file half of this as `PLW1514`
(`unspecified-encoding`), and enabling it is one entry in `[tool.ruff.lint] select`.
That lives in `pyproject.toml`, which CLAUDE.md section 2 protects and for which no
`approved:*` label exists, so the file half stands in for the rule until the owner
enables it and mirrors its call set deliberately.

The subprocess half has no ruff equivalent and is not redundant with it. `text=True`
without `encoding=` decodes the child's stdout through the same locale default --
`Popen.__init__` calls `_text_encoding()` for exactly the reason `open()` does. It is
listed second here because the enumeration in #102 counted only file I/O and missed
it; `python -X warn_default_encoding` found seven live sites, including git output
read by `.claude/hooks/`.

The whole file is a static stand-in for `-X warn_default_encoding`, which is the
authority. That flag cannot be the gate because it only reports lines that actually
EXECUTE, and most of these are error paths.
"""

from __future__ import annotations

import ast
import subprocess
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Ruff's PLW1514 call set. `open` also matches as an attribute so that
# `Path(p).open(...)` and `io.open(...)` are covered.
BARE_FILE_CALLS = frozenset({"open"})
ATTR_FILE_CALLS = frozenset(
    {
        "open",
        "read_text",
        "write_text",
        "NamedTemporaryFile",
        "TemporaryFile",
        "SpooledTemporaryFile",
    }
)

# Modules whose `.open()` is not a text handle: it returns a file descriptor, opens
# a browser, or defaults to binary. Ruff separates these by type inference, which an
# AST walk cannot do, so they are listed. Only a bare `module.open(...)` is matched --
# `os.path`-style attribute chains fall through and are still checked.
NON_TEXT_OPEN_OWNERS = frozenset(
    {"os", "webbrowser", "zipfile", "tarfile", "gzip", "bz2", "lzma", "socket", "dbm", "shelve"}
)

# The subprocess entry points that can hand back decoded text. `call` and
# `check_call` are absent on purpose: they return only a status code, so no decode
# happens and an `encoding=` on them would be cargo.
SUBPROCESS_CALLS = frozenset({"run", "Popen", "check_output"})

# Any one of these switches a subprocess call into text mode; `Popen.__init__` then
# resolves `encoding=None` through `locale.getencoding()`. `errors=` counts, which is
# the non-obvious one -- it turns on text mode all by itself.
TEXT_MODE_KEYWORDS = frozenset({"text", "universal_newlines", "errors"})

# A floor, not a count. It exists so that a broken file list -- a `git ls-files` that
# returns nothing, a checkout with no history -- fails loudly instead of reporting a
# clean scan of zero files, which is the shape of a gate that silently stopped gating.
MIN_FILES_SCANNED = 100


def _tracked_python_files() -> list[Path]:
    proc = subprocess.run(
        ["git", "ls-files", "-z", "--", "*.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=True,
    )
    return [ROOT / name for name in proc.stdout.split("\0") if name]


def _names_an_encoding(call: ast.Call) -> bool:
    return any(keyword.arg == "encoding" for keyword in call.keywords)


def _forwards_kwargs(call: ast.Call) -> bool:
    """`**kwargs` may carry an encoding; the AST cannot see whether it does, and a
    false positive on a passthrough wrapper is not worth the certainty."""
    return any(keyword.arg is None for keyword in call.keywords)


def _opens_in_binary_mode(call: ast.Call) -> bool:
    """True when a literal mode argument carries `b`. An unknown mode is treated as
    text: guessing the other way would let `open(path, mode)` opt itself out."""
    for keyword in call.keywords:
        if keyword.arg == "mode":
            value = keyword.value
            return (
                isinstance(value, ast.Constant)
                and isinstance(value.value, str)
                and "b" in value.value
            )
    # Positional mode: second argument of the builtin `open(file, mode)`, first of
    # every attribute form -- `Path.open(mode)`, `NamedTemporaryFile(mode)`.
    index = 0 if isinstance(call.func, ast.Attribute) else 1
    if len(call.args) > index:
        arg = call.args[index]
        return isinstance(arg, ast.Constant) and isinstance(arg.value, str) and "b" in arg.value
    return False


def _is_unencoded_file_call(call: ast.Call) -> bool:
    func = call.func
    if isinstance(func, ast.Attribute):
        if func.attr not in ATTR_FILE_CALLS:
            return False
        if (
            func.attr == "open"
            and isinstance(func.value, ast.Name)
            and func.value.id in NON_TEXT_OPEN_OWNERS
        ):
            return False
    elif isinstance(func, ast.Name):
        if func.id not in BARE_FILE_CALLS:
            return False
    else:
        return False
    if _names_an_encoding(call) or _forwards_kwargs(call):
        return False
    return not _opens_in_binary_mode(call)


def _is_unencoded_subprocess_call(call: ast.Call) -> bool:
    func = call.func
    name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
    if name not in SUBPROCESS_CALLS:
        return False
    if _names_an_encoding(call) or _forwards_kwargs(call):
        return False
    # A literal False is not text mode; a non-literal is unknowable, so it counts.
    return any(
        keyword.arg in TEXT_MODE_KEYWORDS
        and not (isinstance(keyword.value, ast.Constant) and keyword.value.value in (False, None))
        for keyword in call.keywords
    )


def _offenders(path: Path, predicate: Callable[[ast.Call], bool]) -> list[str]:
    source = path.read_text(encoding="utf-8")
    found: list[str] = []
    for node in ast.walk(ast.parse(source, filename=str(path))):
        if isinstance(node, ast.Call) and predicate(node):
            rendered = ast.get_source_segment(source, node) or ""
            found.append(f"{path.relative_to(ROOT)}:{node.lineno}: {' '.join(rendered.split())}")
    return found


def _scan(predicate: Callable[[ast.Call], bool]) -> list[str]:
    files = _tracked_python_files()
    assert len(files) >= MIN_FILES_SCANNED, (
        f"only {len(files)} tracked .py files found under {ROOT}; this scan is not "
        f"looking at the repository it thinks it is"
    )
    offenders: list[str] = []
    for path in files:
        offenders.extend(_offenders(path, predicate))
    return offenders


def test_every_text_file_call_site_names_its_encoding() -> None:
    offenders = _scan(_is_unencoded_file_call)
    assert not offenders, (
        "text file I/O without an explicit encoding resolves through the platform "
        "locale, so these read cp1252 on a Windows checkout (see #102). Add "
        'encoding="utf-8"; pass mode="rb"/"wb" if the handle is really binary:\n'
        + "\n".join(offenders)
    )


def test_every_subprocess_text_decode_names_its_encoding() -> None:
    offenders = _scan(_is_unencoded_subprocess_call)
    assert not offenders, (
        "these subprocess calls decode the child's output through the platform "
        "locale, so a filename or a git message that is not cp1252-clean is "
        'mojibake on Windows (see #102). Add encoding="utf-8", or drop text mode '
        "and decode the bytes yourself:\n" + "\n".join(offenders)
    )
