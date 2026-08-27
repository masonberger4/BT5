"""M2 vector lane: vector I/O, CDS/backbone partitioning, and construct assembly.

The lane's job is to turn "a plasmid map plus a designed CDS" into the single
`Construct` every rule is evaluated against, without ever letting a backbone base
change and without silently guessing anything the input file did not say --
topology, genetic code, and the 5'UTR are all required-or-reported, never
defaulted.
"""

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
from bt5.vector.locations import (
    LocationError,
    ParsedLocation,
    interval_to_location,
    location_to_interval,
    parts_to_location,
)
from bt5.vector.remap import IntervalRemapper

__all__ = [
    "DEFAULT_EXEMPT_KINDS",
    "DEFAULT_EXEMPT_LABELS",
    "Assembly",
    "InsertionSite",
    "IntervalRemapper",
    "LocationError",
    "ParsedLocation",
    "UtrContext",
    "VectorBackbone",
    "VectorError",
    "assemble",
    "backbone_from_record",
    "backbone_to_record",
    "construct_to_record",
    "insertion_site_from_interval",
    "interval_to_location",
    "location_to_interval",
    "parts_to_location",
    "read_fasta",
    "read_genbank",
    "read_snapgene",
    "read_vector",
    "rotate_interval",
    "write_genbank",
]
