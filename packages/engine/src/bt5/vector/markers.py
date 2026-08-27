"""Recognising the parts of a vector that are not the user's gene.

Two families, both of which mislead a naive detector:

  selection markers   AmpR, PuroR and friends are real, long, well-formed ORFs
                      driven by their own promoters. Nothing else in a plasmid
                      looks so much like a gene worth optimising, and nobody
                      wants them back-translated.
  cloning scars       attB/loxP/FRT sites are recombination substrates whose
                      NUCLEOTIDE sequence is what matters. When one sits inside
                      a CDS, redesigning that CDS destroys it -- so it is a
                      liability to report, not an annotation to drop quietly.
"""

from __future__ import annotations

#: Matched case-folded against a feature's label, note, product and gene.
MARKER_TOKENS: tuple[str, ...] = (
    "ampr",
    "kanr",
    "puror",
    "neor",
    "hygr",
    "zeor",
    "bler",
    "blastr",
    "bsdr",
    "cmr",
    "specr",
    "tetr",
    "confers resistance",
    "resistance gene",
    "beta-lactamase",
    "b-lactamase",
    "lactamase",
    "acetyltransferase",
    "aminoglycoside",
    "phosphotransferase",
    "nourseothricin",
)

#: Site-specific recombination substrates. Nucleotide-level by definition.
RECOMBINATION_TOKENS: tuple[str, ...] = (
    "attb",
    "attp",
    "attl",
    "attr",
    "loxp",
    "lox2272",
    "frt",
    "rox",
    "gateway",
)

#: A eukaryotic expression promoter is hundreds of bases; T7, T3 and SP6 are
#: 17-20 bp and drive in-vitro transcription or sequencing, not expression. The
#: length is the cleanest discriminator that does not need a name list.
MIN_EXPRESSION_PROMOTER_BP = 50


def _matches(text: str, tokens: tuple[str, ...]) -> bool:
    lowered = text.casefold()
    return any(token in lowered for token in tokens)


def is_marker(text: str) -> bool:
    """True for a selection marker, from its label or description."""
    return _matches(text, MARKER_TOKENS)


def is_recombination_site(text: str) -> bool:
    """True for an att/lox/FRT-style recombination site."""
    return _matches(text, RECOMBINATION_TOKENS)
