from pathlib import Path
from typing import Any, override

import pytest

from ducque.library import Library
from ducque.library.library_handler import AbstractLibraryRepository


class TestLibraryRepository(AbstractLibraryRepository):
    """Mock library that will work with an internal dictionary."""

    def __init__(self, tmp_path):
        # Create lib paths so the 'where' command comes through
        self.ducque_db = Path(tmp_path) / "ducque_library.yml"
        self.user_db = Path(tmp_path) / "user_library.yml"
        self._mock_entries = {
            "DA": "DNA adenosine",
            "DT": "DNA thymidine",
            "DG": "DNA guanosine",
            "DC": "DNA cytidine",
        }

        self.ducque_db.touch()
        self.user_db.touch()

    @override
    def fetch_entries(self) -> dict[str, Any]:
        """Fetch entries from the mock database."""
        return self._mock_entries

    @override
    def add_entry(self, entry: dict[str, Any], override: bool) -> None:
        """Add entry to mock database."""
        self._mock_entries.update(entry)

    @override
    def remove_entry(self, entry_key: str) -> None:
        """Remove entry from mock database."""
        self._mock_entries.pop(entry_key)


@pytest.fixture
def snapshot_main_directory() -> Path:
    return Path("tests/test-data/snapshots")


@pytest.fixture
def library_fixture(tmp_path) -> Library:
    return Library(TestLibraryRepository(tmp_path))
