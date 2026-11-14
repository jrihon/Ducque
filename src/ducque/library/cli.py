import typer

from ducque.library.library import Library

library_app = typer.Typer(help="Manage the nucleotide library.")


@library_app.command()
def list():
    """List all entries currently existing in the library."""
    a = Library().list_entries()
    print("listing all entries ...")
    print(a)


@library_app.command()
def add():
    """Add a molecule to the library."""
    Library()
    print("molecule added")


@library_app.command()
def recreate():
    """Recreate a molecule, pdb-formatted, from the library."""
    Library()
    print("molecule recreated from : ")


@library_app.command()
def where() -> None:
    """Tell the user where the ducque libraries are installed."""

    void = "⭕"
    correct = "✅"
    ducque_lib, user_lib = Library().where()

    print(f"{correct} directory exists at '{ducque_lib}'")

    if user_lib is not None:
        print(f"{correct} directory created at '{user_lib}'")
    else:
        print(f"{void} directory not yet created.")
