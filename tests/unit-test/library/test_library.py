from ducque.library import Library


def test_where(library_fixture: Library):

    ducque_lib, user_lib = library_fixture.where()
    assert isinstance(ducque_lib, str)
    assert isinstance(user_lib, str)


def test_library_addition(library_fixture: Library):
    """Test if molecule can be added to the library."""


def test_library_removal(library_fixture: Library):
    """Test if molecule can be removed from the library."""
