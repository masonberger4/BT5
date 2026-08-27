"""Hypothesis strategies shared by the invariant suite."""

from __future__ import annotations

from hypothesis import strategies as st

# Selenocysteine (U) and pyrrolysine (O) are deliberately excluded: they are not
# encodable by standard back-translation and must be rejected by the input
# validator, not silently mistranslated.
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

proteins = st.text(alphabet=STANDARD_AA, min_size=1, max_size=120)


def proteins_with_met() -> st.SearchStrategy[str]:
    """Proteins beginning with the initiator methionine."""
    return st.builds(lambda rest: "M" + rest, st.text(alphabet=STANDARD_AA, max_size=120))


#: Proteins whose naive back-translation is guaranteed to produce perfect
#: nucleotide repeats -- antibodies, linkers, tags, tandem 2A peptides.
repetitive_proteins = st.one_of(
    st.builds(lambda n: "M" + "GGGGS" * n, st.integers(min_value=2, max_value=8)),
    st.builds(lambda n: "M" + "H" * n, st.integers(min_value=6, max_value=12)),
    st.builds(lambda n: "M" + ("DYKDDDDK" * n), st.integers(min_value=2, max_value=4)),
)

dna_bases = st.sampled_from("ACGT")
dna = st.text(alphabet="ACGT", min_size=3, max_size=300)


@st.composite
def in_frame_dna(draw: st.DrawFn) -> str:
    """DNA whose length is a multiple of three."""
    n = draw(st.integers(min_value=1, max_value=100))
    return "".join(draw(st.text(alphabet="ACGT", min_size=3, max_size=3)) for _ in range(n))
