"""M3 score lane: the null, the percentile, and the candidate gallery.

Nothing here reports a predicted expression level. It reports where a design
sits against a distribution of random synonymous variants of the same protein in
the same construct -- a claim the evidence supports -- and which candidates are
genuinely different from one another.
"""

from bt5.score.conflicts import detect_conflicts, hard_versus_soft
from bt5.score.distance import (
    codon_distance,
    distance_matrix,
    nucleotide_distance,
    pairwise_minimum,
)
from bt5.score.gallery import (
    G4_MIN_PAIRWISE_DISTANCE,
    MAX_GALLERY,
    MIN_GALLERY,
    Gallery,
    SweepPoint,
    build_gallery,
    greedy_max_min,
    simplex_weights,
    sweep,
)
from bt5.score.hashing import HASH_LENGTH, design_hash
from bt5.score.null import (
    DEFAULT_NULL_N,
    NullDistribution,
    NullKind,
    band_deviation,
    normalise,
    null_distribution,
    percentile_of,
    synonymous_variant,
)
from bt5.score.order import (
    DEFAULT_PLATE_NAME,
    IDT_HEADERS,
    PLATE_SIZE,
    OrderEntry,
    OrderError,
    entry_name,
    order_entries,
    plates,
    wells,
    write_csv,
    write_idt_plate,
)
from bt5.score.presets import (
    AAV,
    BACTERIAL,
    LENTIVIRAL,
    PRESETS,
    Preset,
    PresetError,
    ResolvedPreset,
    WeightEntry,
    preset_for,
    resolve,
)
from bt5.score.report import (
    DEFAULT_CONFIDENCE,
    ERROR_FREE_BP,
    QcReport,
    ScreeningBurden,
    build_report,
    render,
    screening_burden,
)
from bt5.score.steering import (
    REPEAT_STEERING_PENALTY,
    SWEEP_AXES,
    blended_scorer,
    gc_fraction,
    live_axes,
)

__all__ = [
    "write_idt_plate",
    "write_csv",
    "wells",
    "plates",
    "order_entries",
    "entry_name",
    "OrderError",
    "OrderEntry",
    "PLATE_SIZE",
    "IDT_HEADERS",
    "DEFAULT_PLATE_NAME",
    "screening_burden",
    "render",
    "build_report",
    "ScreeningBurden",
    "QcReport",
    "ERROR_FREE_BP",
    "DEFAULT_CONFIDENCE",
    "AAV",
    "BACTERIAL",
    "DEFAULT_NULL_N",
    "G4_MIN_PAIRWISE_DISTANCE",
    "HASH_LENGTH",
    "LENTIVIRAL",
    "MAX_GALLERY",
    "MIN_GALLERY",
    "PRESETS",
    "REPEAT_STEERING_PENALTY",
    "SWEEP_AXES",
    "Gallery",
    "NullDistribution",
    "NullKind",
    "Preset",
    "PresetError",
    "ResolvedPreset",
    "SweepPoint",
    "WeightEntry",
    "band_deviation",
    "blended_scorer",
    "build_gallery",
    "codon_distance",
    "design_hash",
    "detect_conflicts",
    "distance_matrix",
    "gc_fraction",
    "greedy_max_min",
    "hard_versus_soft",
    "live_axes",
    "normalise",
    "null_distribution",
    "nucleotide_distance",
    "pairwise_minimum",
    "percentile_of",
    "preset_for",
    "resolve",
    "simplex_weights",
    "sweep",
    "synonymous_variant",
]
