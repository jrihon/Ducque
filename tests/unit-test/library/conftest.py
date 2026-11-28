from pathlib import Path

import pytest

from ducque.library import Library
from ducque.library.library_handler import LibraryRepository


@pytest.fixture
def snapshot_main_directory() -> Path:
    return Path("tests/test-data/snapshots")


@pytest.fixture
def library_fixture(tmp_path) -> Library:
    """Create a mock database location for the LibraryRepository, then pass it to Library()."""

    repository = LibraryRepository(
        ducque_db=Path(tmp_path) / "ducque_library.yml",
        user_db=Path(tmp_path) / "user_library.yml",
    )
    # Initialise the databases
    repository.ducque_db.touch()
    repository.user_db.touch()

    return Library(repository)
