from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, override

from platformdirs import user_data_path
from yaml import safe_dump, safe_load

_DEFAULT_DUCQUE_DB = (
    Path(__file__).resolve().parent.parent.parent / "assets" / "ducque_library.yml"
)


_SECRET_DUCQUE_TOGGLE: bool = True


class AbstractLibraryRepository(ABC):
    ducque_db: Path
    user_db: Path

    @abstractmethod
    def fetch_entries(self) -> dict[str, Any]: ...

    @abstractmethod
    def add_entry(self, entry: dict[str, Any], override: bool) -> None: ...

    @abstractmethod
    def remove_entry(self, entry_key: str) -> None: ...


class LibraryRepository(AbstractLibraryRepository):
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

    def __init__(self, ducque_db: Path | None = None, user_db: Path | None = None):
        self.ducque_db: Path = ducque_db or _DEFAULT_DUCQUE_DB
        self.user_db = (
            user_db or user_data_path("ducque") / "user_assets" / "ducque_library.yml"
        )

    def _fetch_ducque_db(self) -> dict[str, Any]:
        """Fetch only entries from the ducque db"""
        with open(self.ducque_db, "r") as db:
            return safe_load(db.read())

    def _fetch_user_db(self) -> dict[str, Any]:
        """Fetch only entries from the user-defined db"""
        if self.user_db.exists():
            with open(self.user_db, "r") as db:
                return safe_load(db.read())
        else:
            return {}

    @override
    def fetch_entries(self) -> dict[str, Any]:
        """Read all entries from the library"""
        return {**self._fetch_ducque_db(), **self._fetch_user_db()}

    @staticmethod
    def _write_to_db(database_path: Path, content: str):
        """Apply the change to the targetted database."""
        with database_path.open("w") as db:
            db.write(content)

    def _write_to_ducque_lib(self, content: dict[str, Any]):
        """Write to the ducque database"""
        content_yaml = safe_dump(content)
        self._write_to_db(self.ducque_db, content_yaml)

    def _write_to_user_lib(self, content: dict[str, Any]):
        """Write to the user-defined database."""
        content_yaml = safe_dump(content)

        if not self.user_db.exists():
            print(f"Creating user-defined library at location '{self.user_db}'.")
            self.user_db.parent.mkdir(parents=True, exist_ok=True)
            self.user_db.touch()

        self._write_to_db(self.user_db, content_yaml)

    @override
    def add_entry(self, entry: dict[str, Any], override: bool = False) -> None:
        """Handles writing (addition of entry) to the library.

        Merge the two objects together, then write to disk.
        """
        print(f"Adding {entry}")

        force_msg = "Use the '-f/--force' to override the entry."
        user_db = self._fetch_user_db()
        ducque_db = self._fetch_ducque_db()

        if _SECRET_DUCQUE_TOGGLE:
            if not entry in ducque_db or override:
                self._write_to_ducque_lib(entry)
            else:
                print(f"Entry already exists in the Ducque's  library. {force_msg}")
        else:
            if not entry in user_db or override:
                self._write_to_user_lib(entry)
            else:
                print(
                    f"Entry already exists in the User's defined library. {force_msg}"
                )

    @override
    def remove_entry(self, entry_key: str) -> None:
        """Handles writing (removal of entry)to the library.

        Remove the entry from the object, then write to disk.
        """
        print(f"Deleting {entry_key}")

        user_db = self._fetch_user_db()
        ducque_db = self._fetch_ducque_db()

        # Try the user library
        if entry_key in user_db:
            _ = user_db.pop(entry_key)
            self._write_to_user_lib(user_db)
            return

        # If developers are updating the library, a toggle needs to be active.
        # Users cannot delete entries from the original database.
        elif entry_key in ducque_db and _SECRET_DUCQUE_TOGGLE:
            _ = ducque_db.pop(entry_key)
            self._write_to_user_lib(ducque_db)
            return

        else:
            print(f"Entry with key {entry_key} does not exist in the ducque library.")
