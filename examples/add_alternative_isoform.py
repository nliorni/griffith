"""Add an alternative isoform to an existing GTF file."""

from __future__ import annotations

from griffith import Griffith


def main() -> None:
    gf = Griffith.from_file(
        "annotation.gtf",
        dialect="gtf",
        attribute_columns=[
            "gene_id",
            "transcript_id",
            "gene_name",
            "transcript_name",
            "exon_number",
            "gene_biotype",
            "transcript_biotype",
        ],
    )

    gf_alt = (
        gf.add_feature(
            seqid="chr1",
            source="griffith",
            feature="transcript",
            start=12000,
            end=18000,
            strand="+",
            attributes={
                "gene_id": "ENSG00000123456",
                "transcript_id": "ENST_NEW_ISOFORM_001",
                "gene_name": "MYGENE",
                "transcript_name": "MYGENE-ALT1",
                "gene_biotype": "protein_coding",
                "transcript_biotype": "protein_coding",
            },
        )
        .add_feature(
            seqid="chr1",
            source="griffith",
            feature="exon",
            start=12000,
            end=12500,
            strand="+",
            attributes={
                "gene_id": "ENSG00000123456",
                "transcript_id": "ENST_NEW_ISOFORM_001",
                "exon_number": "1",
            },
        )
        .add_feature(
            seqid="chr1",
            source="griffith",
            feature="exon",
            start=15000,
            end=18000,
            strand="+",
            attributes={
                "gene_id": "ENSG00000123456",
                "transcript_id": "ENST_NEW_ISOFORM_001",
                "exon_number": "2",
            },
        )
        .sort()
    )

    errors = gf_alt.validate_parent_child()
    print(errors)

    gf_alt.write("annotation.with_alt_isoform.gtf")


if __name__ == "__main__":
    main()
