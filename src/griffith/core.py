"""Core Griffith dataframe object."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import polars as pl

from griffith.attributes import (
    attribute_extract_expr,
    build_attributes_expr,
    format_attributes_py,
    normalise_attribute_value,
)
from griffith.exceptions import GriffithValidationError
from griffith.types import (
    COMMON_ATTRIBUTE_COLUMNS,
    GFF_COLUMNS,
    GTF_TRANSCRIPT_CHILD_FEATURES,
    GTF_TRANSCRIPT_FEATURES,
    Dialect,
)


def _as_list(value: str | Sequence[str]) -> list[str]:
    if isinstance(value, str):
        return [value]
    return list(value)


@dataclass(frozen=True)
class Griffith:
    """
    Polars-backed GFF/GTF annotation table.

    Parameters
    ----------
    frame:
        A Polars LazyFrame containing the 9 canonical GFF/GTF columns and,
        optionally, flattened attribute columns.
    dialect:
        Either ``"gtf"`` or ``"gff3"``.
    attribute_columns:
        Attribute columns used when reconstructing the 9th GFF/GTF field.

    Notes
    -----
    Griffith intentionally keeps the biological hierarchy flat. Genes,
    transcripts, exons, CDS rows, and UTR rows are represented as rows in one
    table. Parent-child relationships are validated through dataframe joins on
    ``gene_id``/``transcript_id`` for GTF or ``ID``/``Parent`` for GFF3.
    """

    frame: pl.LazyFrame
    dialect: Dialect = "gtf"
    attribute_columns: tuple[str, ...] = COMMON_ATTRIBUTE_COLUMNS

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        dialect: Dialect = "gtf",
        attribute_columns: Sequence[str] = COMMON_ATTRIBUTE_COLUMNS,
        parse_attributes: bool = True,
        infer_schema_length: int | None = 1_000,
    ) -> Griffith:
        """
        Lazily read a GFF/GTF file.

        Parameters
        ----------
        path:
            Input GFF/GTF file. Plain text and compressed files supported by
            Polars are accepted.
        dialect:
            ``"gtf"`` or ``"gff3"``.
        attribute_columns:
            Attribute keys to flatten into dataframe columns. Avoid flattening
            every possible attribute in large files.
        parse_attributes:
            If True, selected keys from column 9 are extracted as columns.
        infer_schema_length:
            Number of rows used by Polars for schema inference. Core GFF/GTF
            columns are explicitly typed through schema_overrides.
        """
        attrs = tuple(dict.fromkeys(attribute_columns))

        schema_overrides = {
            "seqid": pl.Utf8,
            "source": pl.Utf8,
            "feature": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "score": pl.Utf8,
            "strand": pl.Utf8,
            "phase": pl.Utf8,
            "attributes": pl.Utf8,
        }

        frame = pl.scan_csv(
            path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=list(GFF_COLUMNS),
            schema_overrides=schema_overrides,
            null_values=["."],
            infer_schema_length=infer_schema_length,
        )

        if parse_attributes and attrs:
            frame = frame.with_columns([attribute_extract_expr(key) for key in attrs])

        return cls(frame=frame, dialect=dialect, attribute_columns=attrs)

    @classmethod
    def from_frame(
        cls,
        frame: pl.DataFrame | pl.LazyFrame,
        *,
        dialect: Dialect = "gtf",
        attribute_columns: Sequence[str] = COMMON_ATTRIBUTE_COLUMNS,
    ) -> Griffith:
        """Create a Griffith object from an existing Polars frame."""
        lazy_frame = frame.lazy() if isinstance(frame, pl.DataFrame) else frame
        return cls(
            frame=lazy_frame,
            dialect=dialect,
            attribute_columns=tuple(dict.fromkeys(attribute_columns)),
        )

    def collect(self, *, engine: str = "auto") -> pl.DataFrame:
        """
        Materialise the current LazyFrame.

        For large files, prefer ``write()`` over ``collect()``.
        """
        return self.frame.collect(engine=engine)

    def head(self, n: int = 5) -> pl.DataFrame:
        """Return the first n rows as a Polars DataFrame."""
        return self.frame.head(n).collect()

    def schema(self) -> dict[str, pl.DataType]:
        """Return the current lazy schema without materialising the data."""
        schema = self.frame.collect_schema()
        return dict(zip(schema.names(), schema.dtypes(), strict=True))

    def _schema_names(self) -> list[str]:
        return list(self.frame.collect_schema().names())

    def _copy(
        self,
        frame: pl.LazyFrame,
        *,
        attrs: Sequence[str] | None = None,
    ) -> Griffith:
        return Griffith(
            frame=frame,
            dialect=self.dialect,
            attribute_columns=tuple(dict.fromkeys(attrs or self.attribute_columns)),
        )

    def _ensure_columns(self, columns: Iterable[str]) -> pl.LazyFrame:
        existing = set(self._schema_names())
        additions = [
            pl.lit(None, dtype=pl.Utf8).alias(column)
            for column in columns
            if column not in existing
        ]

        if not additions:
            return self.frame

        return self.frame.with_columns(additions)

    def with_columns(self, *exprs: pl.Expr, **named_exprs: pl.Expr) -> Griffith:
        """Polars-style column addition/transformation."""
        return self._copy(self.frame.with_columns(*exprs, **named_exprs))

    def select(self, *exprs: str | pl.Expr) -> pl.LazyFrame:
        """Return a Polars LazyFrame selection for advanced users."""
        return self.frame.select(*exprs)

    def filter(self, predicate: pl.Expr) -> Griffith:
        """
        Generic Polars-expression filter.

        Examples
        --------
        >>> gf.filter(pl.col("feature") == "exon")
        """
        return self._copy(self.frame.filter(predicate))

    def subset(
        self,
        *,
        seqid: str | Sequence[str] | None = None,
        feature: str | Sequence[str] | None = None,
        attributes: Mapping[str, str | Sequence[str]] | None = None,
    ) -> Griffith:
        """Subset by chromosome/contig, feature type, or flattened attributes."""
        predicate = pl.lit(True)

        if seqid is not None:
            predicate = predicate & pl.col("seqid").is_in(_as_list(seqid))

        if feature is not None:
            predicate = predicate & pl.col("feature").is_in(_as_list(feature))

        if attributes is not None:
            for key, value in attributes.items():
                predicate = predicate & pl.col(key).is_in(_as_list(value))

        return self.filter(predicate)

    def by_seqid(self, seqid: str | Sequence[str]) -> Griffith:
        """Subset by chromosome/contig name."""
        return self.subset(seqid=seqid)

    def by_feature(self, feature: str | Sequence[str]) -> Griffith:
        """Subset by feature type."""
        return self.subset(feature=feature)

    def by_attribute(self, key: str, value: str | Sequence[str]) -> Griffith:
        """Subset by one flattened attribute."""
        return self.subset(attributes={key: value})

    def sort(
        self,
        by: Sequence[str] = ("seqid", "start", "end", "feature"),
    ) -> Griffith:
        """Sort annotation rows."""
        return self._copy(self.frame.sort(list(by)))

    def add_feature(
        self,
        *,
        seqid: str,
        source: str,
        feature: str,
        start: int,
        end: int,
        score: str | float | int | None = None,
        strand: str | None = None,
        phase: str | int | None = None,
        attributes: Mapping[str, Any] | None = None,
    ) -> Griffith:
        """
        Add one feature row.

        The attributes mapping is written both to the flattened attribute
        columns and to the raw ``attributes`` column. For adding thousands of
        features, prefer ``add_features()`` with a Polars DataFrame/LazyFrame.
        """
        if start > end:
            raise ValueError("Feature start must be <= end.")

        attrs = dict(attributes or {})
        new_attribute_columns = tuple(
            dict.fromkeys([*self.attribute_columns, *attrs.keys()])
        )

        expanded_frame = self._ensure_columns(new_attribute_columns)
        existing_columns = list(expanded_frame.collect_schema().names())
        row: dict[str, Any] = dict.fromkeys(existing_columns)

        row.update(
            {
                "seqid": seqid,
                "source": source,
                "feature": feature,
                "start": start,
                "end": end,
                "score": score,
                "strand": strand,
                "phase": phase,
                "attributes": format_attributes_py(attrs, self.dialect),
            }
        )

        for key, value in attrs.items():
            row[key] = normalise_attribute_value(value)

        new_row = pl.DataFrame([row]).lazy()
        combined = pl.concat([expanded_frame, new_row], how="diagonal_relaxed")

        return self._copy(combined, attrs=new_attribute_columns)

    def add_features(
        self,
        frame: pl.DataFrame | pl.LazyFrame,
        *,
        attribute_columns: Sequence[str] | None = None,
        rebuild_attributes: bool = True,
    ) -> Griffith:
        """
        Append multiple rows from a Polars DataFrame/LazyFrame.

        Parameters
        ----------
        frame:
            Incoming features. It should contain the 9 core GFF/GTF columns or
            enough flattened attribute columns to rebuild ``attributes``.
        attribute_columns:
            Additional attribute columns present in ``frame``.
        rebuild_attributes:
            If True, rebuild the incoming ``attributes`` column before appending.
        """
        incoming = frame.lazy() if isinstance(frame, pl.DataFrame) else frame
        incoming_schema = incoming.collect_schema()
        incoming_names = set(incoming_schema.names())

        attrs = tuple(
            dict.fromkeys([*self.attribute_columns, *(attribute_columns or ())])
        )

        expanded_frame = self._ensure_columns(attrs)

        if rebuild_attributes:
            usable_attrs = [column for column in attrs if column in incoming_names]
            if usable_attrs:
                incoming = incoming.with_columns(
                    build_attributes_expr(usable_attrs, self.dialect)
                )

        combined = pl.concat([expanded_frame, incoming], how="diagonal_relaxed")
        return self._copy(combined, attrs=attrs)

    def remove_by_attribute(self, key: str, values: str | Sequence[str]) -> Griffith:
        """Remove rows where a flattened attribute column matches given values."""
        return self.filter(~pl.col(key).is_in(_as_list(values)))

    def remove_feature(self, predicate: pl.Expr) -> Griffith:
        """Remove rows matching a Polars predicate."""
        return self.filter(~predicate)

    def validate_parent_child(self) -> pl.DataFrame:
        """
        Return rows with broken parent-child relationships.

        For GFF3, each non-null ``Parent`` value must match an existing ``ID``.
        For GTF, transcript children must have a corresponding transcript row,
        and transcript rows must have a corresponding gene row.
        """
        if self.dialect == "gff3":
            return self._validate_gff3_parent_child()

        return self._validate_gtf_parent_child()

    def _validate_gff3_parent_child(self) -> pl.DataFrame:
        required = {"ID", "Parent"}
        missing = required.difference(self._schema_names())

        if missing:
            raise GriffithValidationError(
                f"Cannot validate GFF3 hierarchy. Missing columns: {sorted(missing)}"
            )

        child_columns = [
            column
            for column in ("seqid", "feature", "start", "end", "ID", "Parent")
            if column in self._schema_names()
        ]

        children = (
            self.frame.select(child_columns)
            .filter(pl.col("Parent").is_not_null())
            .with_columns(
                pl.col("Parent")
                .cast(pl.Utf8)
                .str.split(",")
                .alias("_parent_id")
            )
            .explode("_parent_id")
            .with_columns(pl.col("_parent_id").str.strip_chars().alias("_parent_id"))
        )

        parents = (
            self.frame.filter(pl.col("ID").is_not_null())
            .select(pl.col("ID").cast(pl.Utf8).alias("_parent_id"))
            .unique()
        )

        return (
            children.join(parents, on="_parent_id", how="anti")
            .with_columns(pl.lit("missing_parent_id").alias("validation_error"))
            .collect()
        )

    def _validate_gtf_parent_child(self) -> pl.DataFrame:
        required = {"gene_id", "transcript_id"}
        missing = required.difference(self._schema_names())

        if missing:
            raise GriffithValidationError(
                f"Cannot validate GTF hierarchy. Missing columns: {sorted(missing)}"
            )

        transcript_children = self.frame.filter(
            pl.col("feature").is_in(GTF_TRANSCRIPT_CHILD_FEATURES)
            & pl.col("transcript_id").is_not_null()
        )

        transcripts = (
            self.frame.filter(pl.col("feature").is_in(GTF_TRANSCRIPT_FEATURES))
            .select("transcript_id")
            .filter(pl.col("transcript_id").is_not_null())
            .unique()
        )

        missing_transcripts = transcript_children.join(
            transcripts,
            on="transcript_id",
            how="anti",
        ).with_columns(pl.lit("missing_transcript").alias("validation_error"))

        transcript_rows = self.frame.filter(
            pl.col("feature").is_in(GTF_TRANSCRIPT_FEATURES)
            & pl.col("gene_id").is_not_null()
        )

        genes = (
            self.frame.filter(pl.col("feature") == "gene")
            .select("gene_id")
            .filter(pl.col("gene_id").is_not_null())
            .unique()
        )

        missing_genes = transcript_rows.join(
            genes,
            on="gene_id",
            how="anti",
        ).with_columns(pl.lit("missing_gene").alias("validation_error"))

        return pl.concat(
            [missing_transcripts, missing_genes],
            how="diagonal_relaxed",
        ).collect()

    def to_gff(
        self,
        *,
        dialect: Dialect | None = None,
        attribute_columns: Sequence[str] | None = None,
        rebuild_attributes: bool = True,
    ) -> pl.LazyFrame:
        """
        Return a LazyFrame with exactly the 9 canonical GFF/GTF columns.
        """
        selected_dialect = dialect or self.dialect
        schema_names = set(self._schema_names())

        attrs = tuple(
            column
            for column in (attribute_columns or self.attribute_columns)
            if column in schema_names and column not in GFF_COLUMNS
        )

        frame = self.frame

        if rebuild_attributes:
            frame = frame.with_columns(build_attributes_expr(attrs, selected_dialect))

        return frame.select(
            [
                pl.col("seqid").cast(pl.Utf8).fill_null("."),
                pl.col("source").cast(pl.Utf8).fill_null("."),
                pl.col("feature").cast(pl.Utf8).fill_null("."),
                pl.col("start"),
                pl.col("end"),
                pl.col("score").cast(pl.Utf8).fill_null("."),
                pl.col("strand").cast(pl.Utf8).fill_null("."),
                pl.col("phase").cast(pl.Utf8).fill_null("."),
                pl.col("attributes").cast(pl.Utf8).fill_null("."),
            ]
        )

    def write(
        self,
        path: str | Path,
        *,
        dialect: Dialect | None = None,
        attribute_columns: Sequence[str] | None = None,
        rebuild_attributes: bool = True,
    ) -> None:
        """
        Write the annotation as GFF/GTF with no header and tab separator.
        """
        self.to_gff(
            dialect=dialect,
            attribute_columns=attribute_columns,
            rebuild_attributes=rebuild_attributes,
        ).sink_csv(
            path,
            separator="\t",
            include_header=False,
        )

    def write_table(self, path: str | Path) -> None:
        """
        Write the full flattened table, including parsed attributes, as TSV.
        """
        self.frame.sink_csv(path, separator="\t", include_header=True)
