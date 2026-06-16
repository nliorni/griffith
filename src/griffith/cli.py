"""Command-line interface for Griffith."""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from griffith.core import Griffith
from griffith.types import COMMON_ATTRIBUTE_COLUMNS


def _parse_attributes(raw: list[str] | None) -> dict[str, list[str]]:
    parsed: dict[str, list[str]] = {}

    for item in raw or []:
        if "=" not in item:
            raise argparse.ArgumentTypeError(
                f"Invalid attribute filter '{item}'. Expected KEY=VALUE."
            )

        key, value = item.split("=", 1)
        parsed.setdefault(key, []).append(value)

    return parsed


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", type=Path, help="Input GFF/GTF file.")
    parser.add_argument(
        "--dialect",
        choices=["gtf", "gff3"],
        default="gtf",
        help="Annotation dialect. Default: gtf.",
    )
    parser.add_argument(
        "--attr",
        action="append",
        default=[],
        help=(
            "Attribute key to flatten. Can be repeated. "
            "Defaults to common GTF/GFF3 keys."
        ),
    )


def _load(args: argparse.Namespace) -> Griffith:
    attrs = tuple(dict.fromkeys([*COMMON_ATTRIBUTE_COLUMNS, *args.attr]))
    return Griffith.from_file(
        args.input,
        dialect=args.dialect,
        attribute_columns=attrs,
    )


def cmd_validate(args: argparse.Namespace) -> int:
    gf = _load(args)
    errors = gf.validate_parent_child()

    if errors.height == 0:
        print("No parent-child validation errors found.")
        return 0

    print(errors)
    return 1


def cmd_subset(args: argparse.Namespace) -> int:
    gf = _load(args)
    attribute_filters = _parse_attributes(args.where_attr)

    subset = gf.subset(
        seqid=args.seqid,
        feature=args.feature,
        attributes=attribute_filters,
    )
    subset.write(args.output, dialect=args.dialect)
    return 0


def cmd_flatten(args: argparse.Namespace) -> int:
    gf = _load(args)
    gf.write_table(args.output)
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    gf = _load(args)
    stats = (
        gf.frame.group_by("feature")
        .agg(pl.len().alias("n"))
        .sort("n", descending=True)
        .collect()
    )
    print(stats)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="griffith",
        description="Fast Polars-backed utilities for GFF/GTF files.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser(
        "validate",
        help="Validate parent-child relationships.",
    )
    _add_common_args(validate)
    validate.set_defaults(func=cmd_validate)

    subset = subparsers.add_parser(
        "subset",
        help="Subset an annotation and write GFF/GTF output.",
    )
    _add_common_args(subset)
    subset.add_argument("output", type=Path, help="Output GFF/GTF file.")
    subset.add_argument(
        "--seqid",
        action="append",
        help="Chromosome/contig to keep. Can be repeated.",
    )
    subset.add_argument(
        "--feature",
        action="append",
        help="Feature type to keep. Can be repeated.",
    )
    subset.add_argument(
        "--where-attr",
        action="append",
        help="Attribute filter formatted as KEY=VALUE. Can be repeated.",
    )
    subset.set_defaults(func=cmd_subset)

    flatten = subparsers.add_parser(
        "flatten",
        help="Write a flattened TSV table with selected attributes as columns.",
    )
    _add_common_args(flatten)
    flatten.add_argument("output", type=Path, help="Output TSV file.")
    flatten.set_defaults(func=cmd_flatten)

    stats = subparsers.add_parser(
        "stats",
        help="Print simple feature counts.",
    )
    _add_common_args(stats)
    stats.set_defaults(func=cmd_stats)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
