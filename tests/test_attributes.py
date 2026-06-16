from __future__ import annotations

import polars as pl

from griffith.attributes import (
    attribute_extract_expr,
    build_attributes_expr,
    format_attributes_py,
)


def test_format_gtf_attributes() -> None:
    attrs = format_attributes_py(
        {"gene_id": "gene1", "transcript_id": "tx1", "empty": None},
        "gtf",
    )

    assert attrs == 'gene_id "gene1"; transcript_id "tx1";'


def test_format_gff3_attributes() -> None:
    attrs = format_attributes_py(
        {"ID": "tx1", "Parent": "gene1"},
        "gff3",
    )

    assert attrs == "ID=tx1;Parent=gene1"


def test_extract_gtf_attribute() -> None:
    df = pl.DataFrame(
        {
            "attributes": [
                'gene_id "gene1"; transcript_id "tx1"; exon_number "1";'
            ]
        }
    )

    out = df.lazy().with_columns(attribute_extract_expr("transcript_id")).collect()

    assert out["transcript_id"][0] == "tx1"


def test_extract_gff3_attribute() -> None:
    df = pl.DataFrame({"attributes": ["ID=exon1;Parent=tx1;Name=exon-A"]})

    out = df.lazy().with_columns(attribute_extract_expr("Parent")).collect()

    assert out["Parent"][0] == "tx1"


def test_build_gtf_attributes_expr() -> None:
    df = pl.DataFrame(
        {
            "gene_id": ["gene1"],
            "transcript_id": ["tx1"],
            "exon_number": ["1"],
        }
    )

    out = df.lazy().with_columns(
        build_attributes_expr(
            ["gene_id", "transcript_id", "exon_number"],
            "gtf",
        )
    ).collect()

    assert out["attributes"][0] == (
        'gene_id "gene1"; transcript_id "tx1"; exon_number "1";'
    )
