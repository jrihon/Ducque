from importlib.abc import Traversable
from importlib.resources import files
from typing import Any

from platformdirs import user_data_path


class LibraryRepository:
    """Handles IO operations of the library.

    This class should be kept simple.
    LibraryRepository is plumbing at the edge.
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
        self.ducque_lib: Traversable = files("ducque") / "assets" / "ducque_library.yml"
        self.user_lib = user_data_path("ducque") / "user_assets" / "ducque_library.yml"

    def fetch_entries(self) -> dict[str, Any]:
        """Read all entries from the library"""

        # Mock entries for now
        return {
            "DA": "DNA adenosine",
            "DT": "DNA thymidine",
            "DG": "DNA guanosine",
            "DC": "DNA cytidine",
        }
        # d = safe_load(path)
        # return d

    def add_entry(self, py_obj: dict[str, Any]) -> None:
        """Handles writing (addition of entry) to the library"""
        print(f"Adding {py_obj}")

    def remove_entry(self, entry: str) -> None:
        """Handles writing (removal of entry)to the library"""
        print(f"Deleting {entry}")
