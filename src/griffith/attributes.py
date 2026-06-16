"""Attribute parsing and serialization helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any
from urllib.parse import quote

import polars as pl

from griffith.types import GFF_COLUMNS, Dialect


def normalise_attribute_value(value: Any) -> str | None:
    """Convert an attribute value to a clean string or None."""
    if value is None:
        return None

    value_str = str(value).strip()

    if value_str in {"", "."}:
        return None

    return value_str


def format_attributes_py(
    attributes: Mapping[str, Any],
    dialect: Dialect,
) -> str:
    """
    Format a Python mapping as a GTF or GFF3 attribute string.

    This function is used for newly added rows. Bulk export is performed with
    Polars expressions via build_attributes_expr().
    """
    clean_items: list[tuple[str, str]] = []

    for key, value in attributes.items():
        clean_value = normalise_attribute_value(value)
        if clean_value is not None:
            clean_items.append((key, clean_value))

    if not clean_items:
        return "."

    if dialect == "gtf":
        return " ".join(
            f'{key} "{value.replace(chr(34), chr(92) + chr(34))}";'
            for key, value in clean_items
        )

    return ";".join(
        f"{key}={quote(value, safe=',.:_-/')}" for key, value in clean_items
    )


def attribute_extract_expr(key: str) -> pl.Expr:
    """
    Extract one attribute key from either GTF-style or GFF3-style column 9.

    Supported examples
    ------------------
    gene_id "ENSG000001";
    transcript_id "ENST000001";
    ID=gene:ABC;
    Parent=transcript:XYZ;
    """
    escaped_key = key.replace(".", r"\.")

    pattern = (
        rf"(?:^|;\s*){escaped_key}"
        rf"\s*(?:=|\s+)"
        rf"\s*\"?([^\";]+)\"?"
    )

    return (
        pl.col("attributes")
        .cast(pl.Utf8)
        .str.extract(pattern, group_index=1)
        .str.strip_chars('" ')
        .alias(key)
    )


def _attribute_piece_expr(column: str, dialect: Dialect) -> pl.Expr:
    value = pl.col(column).cast(pl.Utf8)
    valid = value.is_not_null() & (value != "") & (value != ".")

    if dialect == "gtf":
        return (
            pl.when(valid)
            .then(pl.format(f'{column} "{{}}";', value))
            .otherwise(None)
        )

    return (
        pl.when(valid)
        .then(pl.format(f"{column}={{}}", value))
        .otherwise(None)
    )


def build_attributes_expr(
    attribute_columns: Sequence[str],
    dialect: Dialect,
) -> pl.Expr:
    """
    Reconstruct column 9 from flattened attribute columns as a Polars expression.
    """
    export_columns = [
        column for column in attribute_columns if column not in set(GFF_COLUMNS)
    ]

    if not export_columns:
        return pl.col("attributes").fill_null(".").alias("attributes")

    pieces = [_attribute_piece_expr(column, dialect) for column in export_columns]
    separator = " " if dialect == "gtf" else ";"

    combined = pl.concat_str(
        pieces,
        separator=separator,
        ignore_nulls=True,
    )

    return (
        pl.when(combined.is_null() | (combined == ""))
        .then(pl.lit("."))
        .otherwise(combined)
        .alias("attributes")
    )
