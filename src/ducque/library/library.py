from pathlib import Path
from typing import Any

from ducque.library.errors import DucqueLibraryNotFound
from ducque.library.library_handler import AbstractLibraryRepository, LibraryRepository
from ducque.library.validator import LibaryValidator

RESERVED_DIR = ["sugar", "linker", "base"]


class Library:
    """Library is the core logic and the public interface.

    It holds the rules of the library: which entries exist, merging the built-in and user libraries,
    refusing to overwrite reserved names like sugar, linker and base unless force is set, filtering by group.
    The CLI, builder and transmute should only ever call Library.

    handler : the class that is aware of the underlying operating system, file operations etc.
        Can be easily mocked for testing purposes.
        Defaults to the LibraryRepository class.
    """

    def __init__(self, handler: AbstractLibraryRepository | None = None):
        self._repository = handler or LibraryRepository()
        self.entries: dict[str, Any] = self._repository.fetch_entries()

    def list_entries(self, _group_by: str | None = None) -> str:
        """List all entries within the library.

        Allow to filter by group and collect a structured output.
        """
        collected_keys = list(self.entries.keys())

        return "\n".join(collected_keys)

    def add(self, yaml_file: str, _force_push: bool = False) -> None:
        """Add residue to the library."""
        validator = LibaryValidator(yaml_file)

        if validator.model.entry_name not in self.entries:
            self._repository.add_entry(validator.py_obj, self.entries)

    def remove(self, entry: str, _force_push: bool = False) -> None:
        """Remove residue from the library."""

        # Filter for entry, to check it exists.
        try:
            self.entries[entry]
            self._repository.remove_entry(entry, self.entries)
        except KeyError:
            print(f"Entry {entry} does not exist in the ducque library.")

    def recreate(self, entry: str) -> None:
        """Recreate residue from the library."""

        # Filter for entry
        try:
            self.entries[entry]
        except KeyError:
            print(f"Entry {entry} does not exist.")

        # Recreate pdb from library

    def where(self) -> tuple[str, str | None]:
        """Return the locations of the library, if they exist."""

        ducque_lib = Path(str(self._repository.ducque_db)).absolute()

        if not ducque_lib.exists():
            raise DucqueLibraryNotFound(ducque_lib)

        return (
            ducque_lib.as_posix(),
            (
                self._repository.user_db.absolute().as_posix()
                if self._repository.user_db.is_dir()
                else None
            ),
        )
