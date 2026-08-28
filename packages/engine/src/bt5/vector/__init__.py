"""M2 vector lane: vector I/O, CDS/backbone partitioning, and construct assembly.

The lane's job is to turn "a plasmid map plus a designed CDS" into the single
`Construct` every rule is evaluated against, without ever letting a backbone base
change and without silently guessing anything the input file did not say --
topology, genetic code, and the 5'UTR are all required-or-reported, never
defaulted.
"""

from bt5.vector.annotate import (
    COMMENT_WIDTH,
    ORIGIN_QUALIFIER,
    DesignReport,
    annotate,
)
from bt5.vector.assemble import Assembly, assemble
from bt5.vector.backbone import (
    DEFAULT_EXEMPT_KINDS,
    DEFAULT_EXEMPT_LABELS,
    InsertionSite,
    UtrContext,
    VectorBackbone,
    VectorError,
    insertion_site_from_interval,
    rotate_interval,
)
from bt5.vector.candidates import (
    SiteCandidate,
    cloning_sites,
    find_orfs,
    suggest_insertion_sites,
)
from bt5.vector.findings import (
    hairpin_liabilities,
    identical_exempt_repeats,
    repeat_liabilities,
)
from bt5.vector.gibson import (
    MAX_ARM_BP,
    MIN_ARM_BP,
    MIN_ARM_TM_C,
    HomologyArm,
    JunctionPlan,
    TmConditions,
    insert_shared_repeats,
    junction_points,
    melting_temperature,
    plan_junctions,
    shortest_usable_arm,
)
from bt5.vector.io import (
    backbone_from_record,
    backbone_to_record,
    construct_to_record,
    read_fasta,
    read_genbank,
    read_snapgene,
    read_vector,
    write_genbank,
)
from bt5.vector.kmers import (
    ConstructKmerIndex,
    InvertedRepeat,
    RepeatPair,
    RiskBand,
    repeat_risk,
)
from bt5.vector.locations import (
    LocationError,
    ParsedLocation,
    interval_to_location,
    location_to_interval,
    parts_to_location,
)
from bt5.vector.markers import is_marker, is_recombination_site
from bt5.vector.notes import DesignNote, NoteKind, format_span
from bt5.vector.remap import IntervalRemapper
from bt5.vector.survey import Survey, construct_from_backbone, survey

__all__ = [
    "MIN_ARM_TM_C",
    "MIN_ARM_BP",
    "MAX_ARM_BP",
    "COMMENT_WIDTH",
    "DEFAULT_EXEMPT_KINDS",
    "DEFAULT_EXEMPT_LABELS",
    "ORIGIN_QUALIFIER",
    "Assembly",
    "ConstructKmerIndex",
    "DesignNote",
    "DesignReport",
    "HomologyArm",
    "JunctionPlan",
    "NoteKind",
    "SiteCandidate",
    "InsertionSite",
    "InvertedRepeat",
    "Survey",
    "IntervalRemapper",
    "LocationError",
    "ParsedLocation",
    "RepeatPair",
    "RiskBand",
    "TmConditions",
    "UtrContext",
    "VectorBackbone",
    "VectorError",
    "annotate",
    "assemble",
    "cloning_sites",
    "backbone_from_record",
    "backbone_to_record",
    "construct_from_backbone",
    "construct_to_record",
    "insert_shared_repeats",
    "insertion_site_from_interval",
    "junction_points",
    "find_orfs",
    "format_span",
    "hairpin_liabilities",
    "identical_exempt_repeats",
    "interval_to_location",
    "is_marker",
    "is_recombination_site",
    "location_to_interval",
    "melting_temperature",
    "parts_to_location",
    "plan_junctions",
    "read_fasta",
    "read_genbank",
    "read_snapgene",
    "read_vector",
    "repeat_liabilities",
    "repeat_risk",
    "rotate_interval",
    "shortest_usable_arm",
    "suggest_insertion_sites",
    "survey",
    "write_genbank",
]
