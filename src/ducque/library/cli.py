from pathlib import Path
from typing import Annotated

import typer

from ducque.library.library import Library

library_app = typer.Typer(help="Manage the nucleotide library.")


@library_app.command()
def list():
    """List all entries currently existing in the library."""
    print("listing all entries ...")
    entries = Library().list_entries()
    print(entries)


@library_app.command()
def add(
    yaml_file: Annotated[
        Path,
        typer.Argument(
            exists=True, help="Input file to add molecule to database, YAML format."
        ),
    ],
    override: Annotated[
        bool, typer.Option(help="Override existing keys in the database.")
    ] = False,
):
    """Add a molecule to the library."""

    print(f"Adding {yaml_file}")
    if override:
        print("Override activated!")

    with yaml_file.open("r") as yml:
        Library().add(yml.read(), override)


@library_app.command()
def remove(yaml_file: str):
    """Remove a molecule from the library."""
    print(f"Removing {yaml_file}")
    Library().remove(yaml_file)


@library_app.command()
def recreate(entry: str):
    """Recreate a molecule, pdb-formatted, from the library."""
    Library().recreate(entry)
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
