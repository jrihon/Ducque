from importlib.abc import Traversable
from importlib.resources import files
from pathlib import Path
from typing import Any

from platformdirs import user_data_path
from yaml import safe_load


class LibraryHandler:
    """Handles IO operations of the library.

    This class should be kept simple.
    LibraryHandler is plumbing at the edge.
    It knows how entries are stored (YAML files, which directories, how to read, write and delete them)
    but nothing about what they mean.
    -> It's an implementation detail of Library

    ducque_lib : the library that comes with installing ducque.
            let's ducque know where it installed on the machine.
    user_lib : a location, defined by your OS, for a user-defined library.
            https://platformdirs.readthedocs.io/en/latest/platforms.html#user-directories
    """

    def __init__(self):
        self.initialised = ""
        self.ducque_lib: Traversable = files("ducque") / "assets"
        self.user_lib = user_data_path("ducque") / "user_assets"

    def get_all_entries(self, path: Path) -> dict[str, Any]:
        """Read all entries from the library"""

        d = safe_load(path)
        return d

    def add_entry(self, path: Path) -> None:
        """Handles writing (addition of entry) to the library"""

    def remove_entry(self, path: Path) -> None:
        """Handles writing (removal of entry)to the library"""
