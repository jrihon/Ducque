from pathlib import Path

import pytest

from ducque.library import Library
from ducque.library.library_handler import LibraryRepository


@pytest.fixture
def snapshot_main_directory() -> Path:
    return Path("tests/test-data/snapshots")


@pytest.fixture
def library_handler_fixture(tmp_path) -> LibraryRepository:
    return LibraryRepository(
        ducque_db=Path(tmp_path) / "ducque_db.yml",
        user_db=Path(tmp_path) / "user_db.yml",
    )


@pytest.fixture
def library_fixture(library_handler_fixture: LibraryRepository) -> Library:
    return Library(library_handler_fixture)
