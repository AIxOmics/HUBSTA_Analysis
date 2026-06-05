"""Small runtime helpers for standalone scripts."""

from __future__ import annotations

import logging
import sys
from pathlib import Path


def add_analysis_dir_to_path(file: str) -> Path:
    """Make 10.gsmap_analysis importable when a script is run directly."""
    analysis_dir = Path(file).resolve().parents[1]
    if str(analysis_dir) not in sys.path:
        sys.path.insert(0, str(analysis_dir))
    return analysis_dir


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

