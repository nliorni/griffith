# API overview

## `Griffith.from_file()`

```python
Griffith.from_file(
    path,
    dialect="gtf",
    attribute_columns=(...),
    parse_attributes=True,
)
```

Lazy-read a GFF/GTF file and optionally flatten selected attributes.

## `Griffith.subset()`

```python
gf.subset(
    seqid="chr1",
    feature="exon",
    attributes={"gene_id": "ENSG00000123456"},
)
```

Subset by chromosome/contig, feature type, and attributes.

## `Griffith.add_feature()`

```python
gf.add_feature(
    seqid="chr1",
    source="griffith",
    feature="exon",
    start=100,
    end=200,
    strand="+",
    attributes={"gene_id": "gene1", "transcript_id": "tx1"},
)
```

Append one feature row.

## `Griffith.add_features()`

Append many rows from an existing Polars DataFrame or LazyFrame.

## `Griffith.validate_parent_child()`

Return a Polars DataFrame containing broken parent-child relationships.

## `Griffith.to_gff()`

Return a LazyFrame with exactly the 9 canonical columns.

## `Griffith.write()`

Write GFF/GTF output to disk.

## `Griffith.write_table()`

Write the full flattened table to TSV.
