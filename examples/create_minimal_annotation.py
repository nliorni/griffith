"""Create a small annotation from scratch."""

from __future__ import annotations

import polars as pl

from griffith import Griffith


def main() -> None:
    empty = pl.DataFrame(
        schema={
            "seqid": pl.Utf8,
            "source": pl.Utf8,
            "feature": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "score": pl.Utf8,
            "strand": pl.Utf8,
            "phase": pl.Utf8,
            "attributes": pl.Utf8,
            "gene_id": pl.Utf8,
            "transcript_id": pl.Utf8,
        }
    )

    gf = Griffith.from_frame(
        empty,
        dialect="gtf",
        attribute_columns=["gene_id", "transcript_id"],
    )

    gf = (
        gf.add_feature(
            seqid="chr1",
            source="griffith",
            feature="gene",
            start=100,
            end=500,
            strand="+",
            attributes={"gene_id": "gene1"},
        )
        .add_feature(
            seqid="chr1",
            source="griffith",
            feature="transcript",
            start=100,
            end=500,
            strand="+",
            attributes={"gene_id": "gene1", "transcript_id": "tx1"},
        )
        .add_feature(
            seqid="chr1",
            source="griffith",
            feature="exon",
            start=100,
            end=200,
            strand="+",
            attributes={"gene_id": "gene1", "transcript_id": "tx1"},
        )
    )

    gf.write("minimal.gtf")


if __name__ == "__main__":
    main()
