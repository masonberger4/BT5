"""The one exception the design lane raises before it hands off to the solver.

A design can fail for reasons that are neither an infeasible constraint search
(that is `InfeasibleConstraints`) nor a refusal to emit (that is
`VerificationError`): no insertion site could be found, a protein that does not
start with the initiator the assembler assumes, a backbone whose own sequence
makes the design impossible to start. Those are `DesignError`, raised with a
message that names what the caller can do about it.
"""

from __future__ import annotations


class DesignError(Exception):
    """The design could not be set up. The message names the remedy."""
