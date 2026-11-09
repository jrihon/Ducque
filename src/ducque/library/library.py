from pathlib import Path
from typing import Any

from ducque.library.library_handler import LibraryRepository

RESERVED_DIR = ["sugar", "linker", "base"]


def _mock_entries() -> dict[str, Any]:
    return {
        "DA": "DNA adenosine",
        "DT": "DNA thymidine",
        "DG": "DNA guanosine",
        "DC": "DNA cytidine",
    }


class Library:
    """Library is the core logic and the public interface.

    It holds the rules of the library: which entries exist, merging the built-in and user libraries,
    refusing to overwrite reserved names like sugar, linker and base unless force is set, filtering by group.
    The CLI, builder and transmute should only ever call Library.

    handler : the class that is aware of the underlying operating system, file operations etc.
        Can be easily mocked for testing purposes.
    """

    def __init__(self, handler: LibraryRepository | None = None):
        self._repository = handler or LibraryRepository()
        self.entries: dict[str, Any] = {}

    def _fetch_entries(self) -> None:
        """Populate ducque's memory with all possible entries."""
        if not self.entries:
            self.entries = _mock_entries()

    def add(self, use_force: bool = False):
        """Add residue to the library.

        Always add residue
        """
        a = Path()
        self._repository.add_entry(a)

    def remove(self, use_force: bool = False):
        """Remove residue to the library."""
        a = Path()
        self._repository.remove_entry(a)

    def list_entries(self, group_by: str) -> dict[str, Any]:
        """List all entries known by ducque.

        Allow to filter by group and collect a structured output.
        """
        self._fetch_entries()
        return self.entries

    def where(self) -> tuple[str, str | None]:
        """Return the locations of the library, if they exist."""

        return (
            Path(str(self._repository.ducque_lib)).absolute().as_posix(),
            (
                self._repository.user_lib.absolute().as_posix()
                if self._repository.user_lib.is_dir()
                else None
            ),
        )
