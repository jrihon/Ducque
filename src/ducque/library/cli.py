import typer

from ducque.library.library import Library

library_app = typer.Typer(help="Manage the nucleotide library")


@library_app.command()
def list_entries():
    Library()
    print("list_entries")


@library_app.command()
def add():
    Library()
    print("molecule added")


@library_app.command()
def where() -> None:
    Library().where()
