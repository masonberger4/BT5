"""The content hash that goes on the tube label.

Two runs producing two different sequences under one name is how a lab ends up
with two tubes and an irreproducible result, so the hash is over the CONTENT and
travels onto the report, the GenBank note and the order file.

Deliberately not Python's `hash()`: that is salted per process, and CI sets
PYTHONHASHSEED=0 precisely because a value that changes between runs cannot
identify anything. sha256 over a canonical encoding is stable across processes,
machines and Python versions.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping

#: Long enough that a collision within one project is not a practical concern,
#: short enough to write on a tube by hand.
HASH_LENGTH = 12


def design_hash(cds: str, *, context: Mapping[str, object] | None = None) -> str:
    """Stable short hash of a design and the parameters that produced it.

    `context` is sorted before hashing, so two callers that build the same
    dictionary in different orders get the same hash -- otherwise the hash would
    identify the caller's iteration order as much as the design.
    """
    payload = {
        "cds": cds.upper(),
        "context": json.loads(json.dumps(context or {}, sort_keys=True, default=str)),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:HASH_LENGTH]
