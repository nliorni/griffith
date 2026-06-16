from __future__ import annotations

from pathlib import Path

import polars as pl

from griffith import Griffith

DATA = Path(__file__).parent / "data"


def test_from_file_parses_gtf_attributes() -> None:
    gf = Griffith.from_file(DATA / "toy.gtf")
    df = gf.collect()

    assert df.height == 6
    assert "gene_id" in df.columns
    assert "transcript_id" in df.columns
    assert df.filter(pl.col("feature") == "transcript")["transcript_id"][0] == "tx1"


def test_subset_by_seqid_and_feature() -> None:
    gf = Griffith.from_file(DATA / "toy.gtf")
    exons_chr1 = gf.subset(seqid="chr1", feature="exon").collect()

    assert exons_chr1.height == 2
    assert exons_chr1["feature"].to_list() == ["exon", "exon"]


def test_add_feature_and_export_gtf() -> None:
    gf = Griffith.from_file(DATA / "toy.gtf")
    gf2 = gf.add_feature(
        seqid="chr1",
        source="griffith",
        feature="transcript",
        start=500,
        end=800,
        strand="+",
        attributes={
            "gene_id": "gene1",
            "transcript_id": "tx_alt",
            "transcript_name": "GENE1-ALT",
        },
    )

    df = gf2.by_attribute("transcript_id", "tx_alt").collect()

    assert df.height == 1
    assert df["attributes"][0] == (
        'gene_id "gene1"; transcript_id "tx_alt"; '
        'transcript_name "GENE1-ALT";'
    )

    exported = gf2.to_gff().collect()
    alt_attrs = exported.filter(
        pl.col("attributes").str.contains("tx_alt")
    )["attributes"][0]

    assert 'transcript_id "tx_alt";' in alt_attrs


def test_validate_gtf_parent_child_ok() -> None:
    gf = Griffith.from_file(DATA / "toy.gtf")
    errors = gf.validate_parent_child()

    assert errors.height == 0


def test_validate_gtf_parent_child_broken() -> None:
    gf = Griffith.from_file(DATA / "broken.gtf")
    errors = gf.validate_parent_child()

    assert errors.height == 2
    assert set(errors["validation_error"].to_list()) == {
        "missing_transcript",
        "missing_gene",
    }


def test_parse_and_validate_gff3_ok() -> None:
    gf = Griffith.from_file(
        DATA / "toy.gff3",
        dialect="gff3",
        attribute_columns=["ID", "Parent", "Name"],
    )
    df = gf.collect()
    errors = gf.validate_parent_child()

    assert df.height == 4
    assert errors.height == 0


def test_validate_gff3_broken() -> None:
    gf = Griffith.from_file(
        DATA / "broken.gff3",
        dialect="gff3",
        attribute_columns=["ID", "Parent", "Name"],
    )
    errors = gf.validate_parent_child()

    assert errors.height == 1
    assert errors["validation_error"][0] == "missing_parent_id"


def test_write_output(tmp_path: Path) -> None:
    gf = Griffith.from_file(DATA / "toy.gtf")
    output = tmp_path / "out.gtf"

    gf.by_seqid("chr1").write(output)

    text = output.read_text()
    assert "chr1" in text
    assert "chr2" not in text
