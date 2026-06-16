"""Griffith: fast, Polars-backed GFF/GTF manipulation."""

from griffith.core import Griffith
from griffith.exceptions import GriffithError, GriffithValidationError

__all__ = [
    "Griffith",
    "GriffithError",
    "GriffithValidationError",
]

__version__ = "0.1.0"
