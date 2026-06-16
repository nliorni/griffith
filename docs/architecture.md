# Griffith architecture

## Data model

Griffith represents every annotation as a flat Polars `LazyFrame`.

The mandatory GFF/GTF columns are:

```text
seqid, source, feature, start, end, score, strand, phase, attributes
```

Selected attributes from column 9 are flattened into additional columns, for
example:

```text
gene_id, transcript_id, exon_number, ID, Parent, Name
```

## Hierarchy handling

Griffith does not store genes, transcripts, exons, and CDS entries as nested
Python objects. The hierarchy is represented through dataframe columns.

For GTF:

```text
gene row:       gene_id
transcript row: gene_id + transcript_id
exon row:       gene_id + transcript_id
```

For GFF3:

```text
parent row: ID
child row:  Parent
```

Validation is implemented as dataframe joins:

- GTF exon/CDS/UTR rows are anti-joined against transcript rows by
  `transcript_id`.
- GTF transcript rows are anti-joined against gene rows by `gene_id`.
- GFF3 child `Parent` values are exploded and anti-joined against parent `ID`
  values.

## Attribute flattening

Column 9 is preserved as `attributes`, but selected keys are parsed into normal
columns.

This avoids the most common GFF/GTF pain point while preventing excessive memory
use. Flattening every attribute in a 1GB+ annotation can create thousands of
sparse columns and should be avoided.

## Export

`to_gff()` reconstructs the canonical 9-column annotation table. If
`rebuild_attributes=True`, the 9th column is regenerated from flattened
attribute columns.

## Performance model

The default path is:

```text
pl.scan_csv(...) -> lazy filters/transforms -> sink_csv(...)
```

This keeps large files out of memory whenever the query plan supports streaming.
Operations such as global sort are naturally more memory-intensive because they
require comparing rows across the full dataset.
