# Griffith tutorial using the toy annotation files

This tutorial walks through the basic Griffith workflow using the small GTF and GFF3 files distributed with the test suite.

The goal is to show the core operations on a safe toy dataset before moving to large genome annotations:

- lazy reading of GTF/GFF3 files;
- flattening the 9th attribute column into normal dataframe columns;
- filtering by chromosome, feature type, and attributes;
- validating parent-child relationships;
- adding a new alternative isoform;
- writing valid GTF/GFF3 output;
- using the command-line interface.

All examples assume you are running commands from the root of the Griffith repository.

---

## 1. Install Griffith for development

Create and activate a dedicated conda environment:

```bash
conda create -n griffith python=3.11 -c conda-forge
conda activate griffith
```

Install Polars:

```bash
conda install -c conda-forge polars
```

Install Griffith in editable development mode:

```bash
pip install -e ".[dev]"
```

Check that the package is importable:

```bash
python -c "from griffith import Griffith; print('Griffith installed correctly')"
```

Run the tests:

```bash
pytest -q
```

---

## 2. Toy files used in this tutorial

The repository contains four small annotation files:

```text
tests/data/toy.gtf
tests/data/broken.gtf
tests/data/toy.gff3
tests/data/broken.gff3
```

The valid toy GTF contains two genes:

```text
chr1  gene        gene1
chr1  transcript  tx1
chr1  exon        tx1 exon 1
chr1  exon        tx1 exon 2
chr2  gene        gene2
chr2  transcript  tx2
```

The valid toy GFF3 contains one gene, one mRNA, and two exons:

```text
chr1  gene  ID=gene1
chr1  mRNA  ID=tx1;Parent=gene1
chr1  exon  ID=exon1;Parent=tx1
chr1  exon  ID=exon2;Parent=tx1
```

The `broken.*` files are intentionally invalid and are used to demonstrate hierarchy validation.

---

## 3. Load a GTF file lazily

```python
from griffith import Griffith


gf = Griffith.from_file(
    "tests/data/toy.gtf",
    dialect="gtf",
)
```

Griffith uses a Polars `LazyFrame` internally, so loading the file does not immediately materialize the full annotation into memory.

Inspect the first rows:

```python
print(gf.head(10))
```

Inspect the schema:

```python
print(gf.schema())
```

You should see the nine canonical GTF/GFF columns plus flattened attribute columns such as:

```text
gene_id
transcript_id
gene_name
transcript_name
exon_number
ID
Parent
Name
```

The extra columns are parsed from the 9th GTF attribute column.

---

## 4. Count features

Because `gf.frame` is a Polars `LazyFrame`, you can use normal Polars expressions directly.

```python
import polars as pl


feature_counts = (
    gf.frame
    .group_by("feature")
    .agg(pl.len().alias("n"))
    .sort("feature")
    .collect()
)

print(feature_counts)
```

Expected result:

```text
shape: (3, 2)
┌────────────┬─────┐
│ feature    ┆ n   │
╞════════════╪═════╡
│ exon       ┆ 2   │
│ gene       ┆ 2   │
│ transcript ┆ 2   │
└────────────┴─────┘
```

---

## 5. Filter by chromosome and feature type

Extract only exons on `chr1`:

```python
chr1_exons = (
    gf
    .by_seqid("chr1")
    .by_feature("exon")
)

result = chr1_exons.select(
    "seqid",
    "feature",
    "start",
    "end",
    "gene_id",
    "transcript_id",
    "exon_number",
).collect()

print(result)
```

Expected result:

```text
shape: (2, 7)
┌───────┬─────────┬───────┬─────┬─────────┬───────────────┬─────────────┐
│ seqid ┆ feature ┆ start ┆ end ┆ gene_id ┆ transcript_id ┆ exon_number │
╞═══════╪═════════╪═══════╪═════╪═════════╪═══════════════╪═════════════╡
│ chr1  ┆ exon    ┆ 100   ┆ 200 ┆ gene1   ┆ tx1           ┆ 1           │
│ chr1  ┆ exon    ┆ 300   ┆ 400 ┆ gene1   ┆ tx1           ┆ 2           │
└───────┴─────────┴───────┴─────┴─────────┴───────────────┴─────────────┘
```

Write this subset back to a valid GTF file:

```python
chr1_exons.write("chr1.exons.gtf")
```

The output contains the standard nine GTF columns, with the 9th attribute column rebuilt from the flattened attributes.

---

## 6. Filter by attribute

Select all rows belonging to `gene1`:

```python
gene1 = gf.by_attribute("gene_id", "gene1")
print(gene1.collect())
```

Select all rows belonging to transcript `tx1`:

```python
tx1 = gf.by_attribute("transcript_id", "tx1")
print(tx1.collect())
```

For advanced filtering, use native Polars expressions:

```python
import polars as pl


long_features = gf.filter(
    (pl.col("end") - pl.col("start") + 1) > 500
)

print(long_features.collect())
```

---

## 7. Validate parent-child relationships in GTF

In GTF mode, Griffith checks that:

- transcript-level child features such as `exon`, `CDS`, and `UTR` have a corresponding `transcript_id` in a transcript row;
- transcript rows have a `gene_id` matching an existing gene row.

Validate the correct toy GTF:

```python
errors = gf.validate_parent_child()
print(errors)
```

Expected result: zero validation errors.

Now validate the intentionally broken GTF:

```python
broken = Griffith.from_file(
    "tests/data/broken.gtf",
    dialect="gtf",
)

errors = broken.validate_parent_child()

print(
    errors.select(
        "seqid",
        "feature",
        "gene_id",
        "transcript_id",
        "validation_error",
    )
)
```

Expected result:

```text
shape: (2, 5)
┌───────┬────────────┬──────────────┬───────────────┬────────────────────┐
│ seqid ┆ feature    ┆ gene_id      ┆ transcript_id ┆ validation_error   │
╞═══════╪════════════╪══════════════╪═══════════════╪════════════════════╡
│ chr1  ┆ exon       ┆ gene1        ┆ missing_tx    ┆ missing_transcript │
│ chr1  ┆ transcript ┆ missing_gene ┆ tx1           ┆ missing_gene       │
└───────┴────────────┴──────────────┴───────────────┴────────────────────┘
```

This is useful after editing an annotation or appending new isoforms.

---

## 8. Add an alternative isoform to an existing GTF gene

The toy GTF already contains `gene1` and transcript `tx1`. Here we add a second transcript, `tx1_alt`, with two exons.

```python
from griffith import Griffith


gf = Griffith.from_file(
    "tests/data/toy.gtf",
    dialect="gtf",
)

with_alt_isoform = (
    gf
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="transcript",
        start=500,
        end=850,
        strand="+",
        attributes={
            "gene_id": "gene1",
            "transcript_id": "tx1_alt",
            "gene_name": "GENE1",
            "transcript_name": "GENE1-ALT",
        },
    )
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="exon",
        start=500,
        end=550,
        strand="+",
        attributes={
            "gene_id": "gene1",
            "transcript_id": "tx1_alt",
            "exon_number": "1",
        },
    )
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="exon",
        start=700,
        end=850,
        strand="+",
        attributes={
            "gene_id": "gene1",
            "transcript_id": "tx1_alt",
            "exon_number": "2",
        },
    )
    .sort()
)
```

Check that the edited annotation is still valid:

```python
errors = with_alt_isoform.validate_parent_child()
print(errors.height)
```

Expected result:

```text
0
```

Write the updated annotation:

```python
with_alt_isoform.write("toy.with_alt_isoform.gtf")
```

Preview the added transcript:

```python
preview = (
    with_alt_isoform
    .by_attribute("transcript_id", "tx1_alt")
    .select(
        "seqid",
        "source",
        "feature",
        "start",
        "end",
        "strand",
        "attributes",
    )
    .collect()
)

print(preview)
```

---

## 9. Append many features at once

For a few rows, `add_feature()` is convenient. For thousands or millions of rows, build a Polars `DataFrame` or `LazyFrame` and append it with `add_features()`.

```python
import polars as pl


new_features = pl.DataFrame(
    [
        {
            "seqid": "chr1",
            "source": "griffith",
            "feature": "transcript",
            "start": 500,
            "end": 850,
            "score": None,
            "strand": "+",
            "phase": None,
            "gene_id": "gene1",
            "transcript_id": "tx1_batch",
            "gene_name": "GENE1",
            "transcript_name": "GENE1-BATCH",
        },
        {
            "seqid": "chr1",
            "source": "griffith",
            "feature": "exon",
            "start": 500,
            "end": 550,
            "score": None,
            "strand": "+",
            "phase": None,
            "gene_id": "gene1",
            "transcript_id": "tx1_batch",
            "exon_number": "1",
        },
    ]
)

updated = gf.add_features(
    new_features,
    attribute_columns=[
        "gene_id",
        "transcript_id",
        "gene_name",
        "transcript_name",
        "exon_number",
    ],
)

updated.write("toy.with_batch_features.gtf")
```

`add_features()` is preferred for batch operations because it avoids repeatedly concatenating one row at a time.

---

## 10. Work with GFF3 files

GFF3 uses `ID` and `Parent` instead of the GTF-style `gene_id` and `transcript_id` hierarchy.

Load the toy GFF3:

```python
gff3 = Griffith.from_file(
    "tests/data/toy.gff3",
    dialect="gff3",
    attribute_columns=["ID", "Parent", "Name"],
)

print(gff3.head())
```

Validate parent-child relationships:

```python
errors = gff3.validate_parent_child()
print(errors.height)
```

Expected result:

```text
0
```

Validate the intentionally broken GFF3:

```python
broken_gff3 = Griffith.from_file(
    "tests/data/broken.gff3",
    dialect="gff3",
    attribute_columns=["ID", "Parent", "Name"],
)

errors = broken_gff3.validate_parent_child()

print(
    errors.select(
        "seqid",
        "feature",
        "ID",
        "Parent",
        "validation_error",
    )
)
```

Expected result:

```text
shape: (1, 5)
┌───────┬─────────┬───────┬────────────┬───────────────────┐
│ seqid ┆ feature ┆ ID    ┆ Parent     ┆ validation_error  │
╞═══════╪═════════╪═══════╪════════════╪═══════════════════╡
│ chr1  ┆ exon    ┆ exon1 ┆ missing_tx ┆ missing_parent_id │
└───────┴─────────┴───────┴────────────┴───────────────────┘
```

---

## 11. Add a new GFF3 mRNA and exon

```python
gff3_updated = (
    gff3
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="mRNA",
        start=500,
        end=850,
        strand="+",
        attributes={
            "ID": "tx1_alt",
            "Parent": "gene1",
            "Name": "GENE1-ALT",
        },
    )
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="exon",
        start=500,
        end=850,
        strand="+",
        attributes={
            "ID": "exon_alt_1",
            "Parent": "tx1_alt",
        },
    )
    .sort()
)

print(gff3_updated.validate_parent_child().height)

gff3_updated.write("toy.with_alt_isoform.gff3")
```

Expected result:

```text
0
```

---

## 12. Export a flattened table

For inspection or downstream analysis, write the full flattened table as TSV:

```python
gf.write_table("toy.flattened.tsv")
```

This keeps the canonical GTF/GFF columns and the parsed attributes as additional columns.

For example, the flattened table contains columns such as:

```text
seqid
source
feature
start
end
score
strand
phase
attributes
gene_id
transcript_id
gene_name
transcript_name
exon_number
```

---

## 13. Command-line usage

After installing Griffith, the `griffith` command is available.

### Count feature types

```bash
griffith stats tests/data/toy.gtf --dialect gtf
```

Expected result:

```text
shape: (3, 2)
┌────────────┬─────┐
│ feature    ┆ n   │
╞════════════╪═════╡
│ exon       ┆ 2   │
│ gene       ┆ 2   │
│ transcript ┆ 2   │
└────────────┴─────┘
```

The exact row order can vary depending on sorting options and Polars version.

### Validate a correct GTF

```bash
griffith validate tests/data/toy.gtf --dialect gtf
```

Expected result:

```text
No parent-child validation errors found.
```

### Validate a broken GTF

```bash
griffith validate tests/data/broken.gtf --dialect gtf
```

Expected result: the command prints the invalid rows and exits with status code `1`.

### Subset chr1 exons

```bash
griffith subset tests/data/toy.gtf toy.chr1.exons.gtf \
  --dialect gtf \
  --seqid chr1 \
  --feature exon
```

### Filter by attribute

```bash
griffith subset tests/data/toy.gtf toy.gene1.gtf \
  --dialect gtf \
  --where-attr gene_id=gene1
```

### Flatten attributes to TSV

```bash
griffith flatten tests/data/toy.gtf toy.flattened.tsv \
  --dialect gtf \
  --attr gene_id \
  --attr transcript_id \
  --attr gene_name \
  --attr transcript_name \
  --attr exon_number
```

---

## 14. Practical notes before using large files

For large annotations, avoid this pattern unless you really need all rows in memory:

```python
full_df = gf.collect()
```

Prefer lazy filters followed by direct writing:

```python
(
    gf
    .by_seqid("chr1")
    .by_feature("exon")
    .write("chr1.exons.gtf")
)
```

Also avoid flattening every possible attribute key. Pass only the attributes you need:

```python
gf = Griffith.from_file(
    "large.annotation.gtf.gz",
    dialect="gtf",
    attribute_columns=[
        "gene_id",
        "transcript_id",
        "gene_name",
        "exon_number",
    ],
)
```

This is one of the main design choices that keeps Griffith memory-efficient.

---

## 15. Minimal complete example

```python
from griffith import Griffith


gf = Griffith.from_file(
    "tests/data/toy.gtf",
    dialect="gtf",
)

edited = (
    gf
    .by_attribute("gene_id", "gene1")
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="transcript",
        start=500,
        end=850,
        strand="+",
        attributes={
            "gene_id": "gene1",
            "transcript_id": "tx1_alt",
            "gene_name": "GENE1",
        },
    )
    .add_feature(
        seqid="chr1",
        source="griffith",
        feature="exon",
        start=500,
        end=850,
        strand="+",
        attributes={
            "gene_id": "gene1",
            "transcript_id": "tx1_alt",
            "exon_number": "1",
        },
    )
    .sort()
)

errors = edited.validate_parent_child()
assert errors.height == 0

edited.write("gene1.with_alt_isoform.gtf")
```