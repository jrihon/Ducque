# let the package know where it is installed in the system
from importlib.resources import files
from pathlib import Path
from typing import Any

# For user-defined inputs, assigns a `local` directory automatically.
# Different per major OS
from platformdirs import user_data_path

from ducque.library.library_handler import LibraryHandler

BUILTIN_LIB = files("ducque") / "assets"
USER_LIB = user_data_path("ducque") / "user_assets"
RESERVED_DIR = ["sugar", "linker", "base"]


def _mock_entries() -> dict[str, Any]:
    return {
        "DA": "DNA adenosine",
        "DT": "DNA thymidine",
        "DG": "DNA guanosine",
        "DC": "DNA cytidine",
    }


class Library:
    def __init__(self):
        self.entries: dict[str, Any] = {}

    def fetch_entries(self) -> None:
        """Populate ducque's memory with all possible entries."""
        if not self.entries:
            self.entries = _mock_entries()

    def add(self, use_force: bool = False):
        """Add residue to the library.

        Always add residue
        """
        lh = LibraryHandler()
        lh.add_entry(Path())

    def remove(self, use_force: bool = False):
        """Remove residue to the library."""
        lh = LibraryHandler()
        lh.remove_entry(Path())

    def list_entries(self, group_by: str) -> dict[str, Any]:
        """List all entries known by ducque.

        Allow to filter by group and collect a structured output.
        """
        self.fetch_entries()
        return self.entries

    def where(self):
        """Print the locations of the library, if they exist."""
        void = "⭕"
        correct = "✅"
        assets = Path(str(BUILTIN_LIB)).absolute().as_posix()
        user_assets = USER_LIB.absolute().as_posix()

        print(f"{correct} directory exists at '{assets}'")

        if USER_LIB.is_dir():
            print(f"{correct} directory created at '{user_assets}'")
        else:
            print(f"{void} directory not yet created.")
