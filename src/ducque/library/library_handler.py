from pathlib import Path
from typing import Any

from yaml import safe_load


class LibraryHandler:
    """Handles IO operations of the library."""

    def __init__(self):
        self.initialised = ""

    def get_all_entries(self, path: Path) -> dict[str, Any]:
        """Read all entries from the library"""

        d = safe_load(path)
        return d

    def add_entry(self, path: Path) -> None:
        """Handles writing (addition of entry) to the library"""

    def remove_entry(self, path: Path) -> None:
        """Handles writing (removal of entry)to the library"""
