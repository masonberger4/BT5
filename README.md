# BT5

Codon optimization and protein back-translation, evaluated on the **assembled
construct** rather than a free-floating coding sequence.

Given a protein and a vector backbone you already have, BT5 designs a coding
sequence that works *in that backbone* — checking forbidden motifs, GC windows and
repeats across the circular product including the insert/backbone junctions and
the origin, against three simultaneous contexts (E. coli propagation, the
packaging cell, and the target cell).

## What it will not tell you

Codon optimization does not reliably increase expression. Across nine benchmarked
commercial optimizers there was "a roughly equivalent chance that an
algorithm-optimized CDS will increase or diminish recombinant yields", and all
computable design features together explain only ~14% of protein-level variance.

So BT5 reports ranks and percentiles against a random-synonymous null, never a
predicted expression number, and treats "use the native sequence" as a first-class
answer. Every rule carries an evidence badge and a citation; rules resting on
folklore ship disabled.

Where BT5 is genuinely useful is the mechanical half: guaranteeing hard
constraints, doing it in the context of your actual plasmid, and quantifying what
each constraint costs you.

## Status

Foundation (PR #0). See `docs/PLAN.md` for the full design and build sequence.

## Development

```bash
uv venv --python 3.11 .venv && . .venv/bin/activate
uv pip install -e ".[dev,fold,export]"
pytest -q
```

Contributors — including AI sessions — must read `CLAUDE.md` first.
