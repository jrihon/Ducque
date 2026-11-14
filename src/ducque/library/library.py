from pathlib import Path
from typing import Any

from ducque.library.errors import DucqueLibraryNotFound
from ducque.library.library_handler import LibraryRepository
from ducque.library.validator import LibaryValidator

RESERVED_DIR = ["sugar", "linker", "base"]


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
        """Populate the runtime memory with all possible entries."""
        if not self.entries:
            self.entries = self._repository.fetch_entries()

    def add(self, yaml_file: Path, force_push: bool = False):
        """Add residue to the library."""
        validator = LibaryValidator(yaml_file)
        self._repository.add_entry(validator.py_obj)

    def remove(self, entry: str, force_push: bool = False):
        """Remove residue to the library."""
        self._repository.remove_entry(entry)

    def list_entries(self, group_by: str | None = None) -> dict[str, Any]:
        """List all entries within the library.

        Allow to filter by group and collect a structured output.
        """
        self._fetch_entries()
        return self.entries

    def where(self) -> tuple[str, str | None]:
        """Return the locations of the library, if they exist."""

        ducque_lib = Path(str(self._repository.ducque_lib)).absolute()

        if not ducque_lib.exists():
            raise DucqueLibraryNotFound(ducque_lib)

        return (
            ducque_lib.as_posix(),
            (
                self._repository.user_lib.absolute().as_posix()
                if self._repository.user_lib.is_dir()
                else None
            ),
        )
