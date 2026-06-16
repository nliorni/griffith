"""Filter a GTF file, validate relationships, and write output."""

from __future__ import annotations

from griffith import Griffith


def main() -> None:
    gf = Griffith.from_file("annotation.gtf", dialect="gtf")

    result = gf.by_seqid("chr1").by_feature(["gene", "transcript", "exon"])

    errors = result.validate_parent_child()
    if errors.height:
        print("Validation errors:")
        print(errors)

    result.write("chr1.gene_transcript_exon.gtf")


if __name__ == "__main__":
    main()
