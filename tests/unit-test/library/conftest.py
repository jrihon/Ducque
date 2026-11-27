from pathlib import Path
from typing import override

import pytest

from ducque.library import Library
from ducque.library.library_handler import AbstractLibraryRepository


class TestLibraryRepository(AbstractLibraryRepository):
    """Mock library that will work with an internal dictionary."""

    def __init__(self):
        pass

    @override
    def fetch_entries(): ...

    @override
    def add_entry(): ...

    @override
    def remove_entry(): ...


@pytest.fixture
def snapshot_main_directory() -> Path:
    return Path("tests/test-data/snapshots")


@pytest.fixture
def library_handler_fixture() -> TestLibraryRepository:
    return TestLibraryRepository()


@pytest.fixture
def library_fixture(library_handler_fixture: TestLibraryRepository) -> Library:
    return Library(library_handler_fixture)
