"""Shared type aliases and constants."""

from typing import Literal

Dialect = Literal["gtf", "gff3"]

GFF_COLUMNS: tuple[str, ...] = (
    "seqid",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
)

COMMON_ATTRIBUTE_COLUMNS: tuple[str, ...] = (
    "gene_id",
    "transcript_id",
    "gene_name",
    "transcript_name",
    "exon_number",
    "gene_biotype",
    "transcript_biotype",
    "protein_id",
    "ID",
    "Parent",
    "Name",
    "biotype",
)

GTF_TRANSCRIPT_FEATURES: tuple[str, ...] = (
    "transcript",
    "mRNA",
)

GTF_TRANSCRIPT_CHILD_FEATURES: tuple[str, ...] = (
    "exon",
    "CDS",
    "UTR",
    "five_prime_UTR",
    "three_prime_UTR",
    "start_codon",
    "stop_codon",
)
