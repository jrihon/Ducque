from ducque.library.errors import DucqueLibraryNotFound
from ducque.library.library_handler import AbstractLibraryRepository, LibraryRepository
from ducque.library.validator import MoietyTypeValidator

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

    def list_entries(self, _group_by: str | None = None) -> str:
        """List all entries within the library.

        Allow to filter by group and collect a structured output.
        """
        entries = self._repository.fetch_entries()
        collected_keys = list(entries.keys())

        return "\n".join(collected_keys)

    def add(self, yaml_file: str, override: bool) -> None:
        """Add residue to the library."""
        validated_molecule = MoietyTypeValidator(yaml_file)

        self._repository.add_entry(validated_molecule.py_obj, override)

    def remove(self, entry: str) -> None:
        """Remove residue from the library."""
        self._repository.remove_entry(entry)

    def recreate(self, entry: str) -> None:
        """Recreate residue from the library."""

        entries = self._repository.fetch_entries()
        if not entry in entries:
            print(f"Entry {entry} does not exist.")
            return

        # Recreate pdb from library
        print(f"Creating entry {entry}")

    def where(self) -> tuple[str, str | None]:
        """Return the locations of the library, if they exist."""

        ducque_lib = self._repository.ducque_db.absolute()

        if not ducque_lib.exists():
            raise DucqueLibraryNotFound(ducque_lib)

        return (
            ducque_lib.as_posix(),
            (
                self._repository.user_db.absolute().as_posix()
                if self._repository.user_db.exists()
                else None
            ),
        )
