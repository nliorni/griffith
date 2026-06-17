# Griffith

**Griffith** is a fast, dataframe-oriented Python library for manipulating
GFF/GTF genome annotation files using [Polars](https://pola.rs/).

It is designed for large annotations where object-per-feature approaches can
become slow or memory-heavy. Griffith keeps the annotation as a flat Polars
`LazyFrame` and represents the biological hierarchy through attributes such as
`gene_id`, `transcript_id`, `ID`, and `Parent`.

> Status: early alpha. The core API is intentionally small and suitable for
> extension.

---

## Why Griffith?

GFF/GTF files are tabular, but many tools expose them as nested Python objects.
That is convenient for small files, but can be expensive for large annotations.
Griffith instead treats annotations as a columnar table:

```text
seqid | source | feature | start | end | score | strand | phase | attributes | gene_id | transcript_id | ID | Parent | ...
```

This makes common operations such as filtering, joining, sorting, validating,
and exporting expressible as Polars operations.

---

## Installation

### With Conda

```bash
git clone https://github.com/nliorni/griffith.git
conda create -n griffith python=3.11 -c conda-forge
conda activate griffith
conda install -c conda-forge polars
cd griffith
pip install -e .
python -c "from griffith import Griffith; print('Griffith installed correctly')"
```

### Development install

```bash
pip install -e '.[dev]'
pytest
```

---

## Quick start

```python
import polars as pl

from griffith import Griffith

# Lazy read; attributes are parsed into selected columns.
gf = Griffith.from_file(
    "annotation.gtf",
    dialect="gtf",
    attribute_columns=[
        "gene_id",
        "transcript_id",
        "gene_name",
        "transcript_name",
        "exon_number",
    ],
)

# Fluent filtering.
chr1_exons = (
    gf
    .by_seqid("chr1")
    .by_feature("exon")
)

# Write back to GTF.
chr1_exons.write("chr1.exons.gtf")
```

---

## Add an alternative isoform

```python
from griffith import Griffith


gf = Griffith.from_file("annotation.gtf", dialect="gtf")

gf_alt = (
    gf
    .add_feature(
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
```

---

## GFF3 example

```python
from griffith import Griffith


gf = Griffith.from_file(
    "annotation.gff3",
    dialect="gff3",
    attribute_columns=["ID", "Parent", "Name", "gene_id"],
)

gf2 = gf.add_feature(
    seqid="chr1",
    source="griffith",
    feature="mRNA",
    start=12000,
    end=18000,
    strand="+",
    attributes={
        "ID": "transcript:NEW001",
        "Parent": "gene:GENE001",
        "Name": "GENE001-ALT1",
    },
)

gf2.write("annotation.with_new_mrna.gff3")
```

---

## Command-line usage

After installation, the `griffith` command is available.

### Validate parent-child relationships

```bash
griffith validate annotation.gtf --dialect gtf
```

### Subset an annotation

```bash
griffith subset annotation.gtf chr1.exons.gtf \
  --dialect gtf \
  --seqid chr1 \
  --feature exon
```

### Flatten selected attributes to TSV

```bash
griffith flatten annotation.gtf annotation.flattened.tsv \
  --dialect gtf \
  --attr gene_id \
  --attr transcript_id \
  --attr gene_name
```

### Feature counts

```bash
griffith stats annotation.gtf --dialect gtf
```

---

## Design principles

1. **No Pandas.** All dataframe operations use Polars.
2. **Lazy by default.** File parsing uses `pl.scan_csv()`.
3. **Flat hierarchy.** Parent-child relationships are represented through
   attributes, not nested Python objects.
4. **Controlled attribute flattening.** Only selected attributes are parsed into
   columns to avoid memory blow-up.
5. **Round-trip export.** Flattened attributes can be reconstructed into a valid
   GTF/GFF3 9th column.

---

## Current limitations

- Attribute parsing is intentionally pragmatic and optimized for common GTF/GFF3
  patterns. Highly irregular attribute strings may require custom preprocessing.
- Global sorting may require more memory than streaming filters because sorting
  is a global operation.
- GFF3 URL-decoding is not yet applied during parsing.
- Full ontology-level validation is not yet implemented.

---

## Roadmap

- Region queries by interval overlap.
- Optional attribute URL decoding/encoding for GFF3.
- Transcript model constructors.
- Streaming-safe batch append helpers.
- Parquet-backed intermediate cache.
- Richer validation rules for gene/transcript/exon/CDS structure.
