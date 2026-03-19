"""
CheckU: UNI56 marker completeness profiling for microbial genomes.
"""

__version__ = "0.2.0"
__author__ = "Juan C. Villada"
__email__ = "jvillada@lbl.gov"
__title__ = "checku"
__description__ = "UNI56 marker completeness profiling for microbial genomes"
__license__ = "BSD-3-Clause"

from .checku import app

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__title__",
    "__description__",
    "__license__",
    "app",
]


def get_version() -> str:
    """Return the package version."""
    return __version__


def main() -> None:
    """Entry point for the CheckU command line interface."""
    app()
