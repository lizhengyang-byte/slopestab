"""Package exceptions."""

from __future__ import annotations


class SlopestabError(Exception):
    """Base class for every error raised by slopestab."""


class GeometryError(SlopestabError, ValueError):
    """Invalid or degenerate geometry (profile, surface, slices)."""


class ConvergenceError(SlopestabError, RuntimeError):
    """A limit-equilibrium solver failed to converge.

    The message always contains the last iterate and the convergence criterion,
    because a silently wrong factor of safety is worse than no answer.
    """

    def __init__(self, method: str, message: str, last_fs: float | None = None) -> None:
        self.method = method
        self.last_fs = last_fs
        detail = f"{method}: {message}"
        if last_fs is not None:
            detail += f" (last iterate Fs={last_fs:.6f})"
        super().__init__(detail)
