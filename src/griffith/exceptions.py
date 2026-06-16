"""Custom exceptions for Griffith."""


class GriffithError(Exception):
    """Base exception for Griffith."""


class GriffithValidationError(GriffithError):
    """Raised when a validation operation cannot be performed."""
