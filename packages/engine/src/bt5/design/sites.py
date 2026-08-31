"""Choose where the new coding sequence goes -- multiple cloning site first.

The common case is NOT re-optimising a CDS in place. It is a user with a
backbone that carries a promoter and a polylinker (a multiple cloning site)
downstream of it, into which a brand-new protein's DNA is inserted. So site
selection is MCS-first:

1. An explicit `InsertionSite`, or an explicit interval, or an explicit label --
   the caller knows exactly where.
2. `cloning_sites(backbone)` -- the promoter-plus-polylinker case, which is the
   common one and which had no production caller before this module.
3. `find_insertion_site()` -- an annotated CDS to replace, the re-optimise case.
4. Otherwise raise, naming the ranked `suggest_insertion_sites` candidates so the
   error tells the user where a site COULD go rather than only that none was
   found.

A note on what E1 can and cannot express: a zero-length pure insertion between
two bases is not representable, because `Interval` rejects `end <= start`. The
MCS case therefore REPLACES the polylinker span with the insert, which is the
real cloning operation anyway -- the polylinker exists to be cut out and
replaced. A backbone whose only site is a zero-length point is out of scope.
"""

from __future__ import annotations

from bt5.core.types import Interval
from bt5.design.errors import DesignError
from bt5.vector.backbone import (
    InsertionSite,
    VectorBackbone,
    VectorError,
    insertion_site_from_interval,
)
from bt5.vector.candidates import cloning_sites, suggest_insertion_sites


def choose_site(
    backbone: VectorBackbone,
    *,
    table_id: int,
    site: InsertionSite | None = None,
    site_interval: Interval | None = None,
    site_label: str | None = None,
) -> InsertionSite:
    """Resolve the insertion site, MCS-first. See the module docstring for order.

    `table_id` is stamped onto an interval- or MCS-derived site so the assembler
    can cross-check it against any `/transl_table` the vector annotates.
    """
    if site is not None:
        return site
    if site_interval is not None:
        return insertion_site_from_interval(
            site_interval, label=site_label or "insert", table_id=table_id
        )
    if site_label is not None:
        # A named annotated CDS to replace -- the caller pointed at it by label.
        return backbone.find_insertion_site(label=site_label)

    mcs = cloning_sites(backbone)
    if mcs:
        # cloning_sites() returns sites sorted by descending span, so the first is
        # the widest polylinker -- the most likely intended cloning site.
        return mcs[0].as_site(table_id=table_id)

    try:
        return backbone.find_insertion_site()
    except VectorError as exc:
        suggestions = suggest_insertion_sites(backbone, table_id=table_id)
        if suggestions:
            named = "; ".join(f"{c.label} ({c.kind}, score {c.score})" for c in suggestions)
            raise DesignError(
                f"no cloning site or annotated CDS to insert into: {exc}. Candidate "
                f"sites, best first: {named}. Pass one as site_interval= or site_label=."
            ) from exc
        raise DesignError(
            f"no insertion site could be found and none could be suggested: {exc}. "
            f"Give an explicit site_interval= (the span to replace) or annotate a CDS "
            f"or a promoter with a downstream polylinker."
        ) from exc
