from importlib.abc import Traversable
from importlib.resources import files
from pathlib import Path
from typing import Any

from platformdirs import user_data_path
from uaclient.yaml import safe_dump


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
        self.ducque_db: Traversable = files("ducque") / "assets" / "ducque_library.yml"
        self.user_db = user_data_path("ducque") / "user_assets" / "ducque_library.yml"

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

    @staticmethod
    def _write_to_db(database_path: Path, content: str):
        """Apply the change to the targetted database."""
        with database_path.open("w") as db:
            db.write(content)

    def _write_to_ducque_lib(self, content: dict[str, Any]):
        """Write to the base database"""

        content_yaml = safe_dump(content)
        ducque_db = Path(str(self.ducque_db))
        self._write_to_db(ducque_db, content_yaml)

    def _write_to_user_lib(self, content: dict[str, Any]):
        """Write to the user defined database."""
        content_yaml = safe_dump(content)
        self._write_to_db(self.user_db, content_yaml)

    def add_entry(self, entry: dict[str, Any], all_entries: dict[str, Any]) -> None:
        """Handles writing (addition of entry) to the library.

        Merge the two objects together, then write to disk.
        """

        print(f"Adding {entry}")

    def remove_entry(self, entry: str, all_entries: dict[str, Any]) -> None:
        """Handles writing (removal of entry)to the library.

        Remove the entry from the object, then write to disk.
        """
        print(f"Deleting {entry}")

        print("Writing deletion into library.")
