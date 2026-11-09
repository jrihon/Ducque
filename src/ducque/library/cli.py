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

    void = "⭕"
    correct = "✅"
    ducque_lib, user_lib = Library().where()

    print(f"{correct} directory exists at '{ducque_lib}'")

    if user_lib is not None:
        print(f"{correct} directory created at '{user_lib}'")
    else:
        print(f"{void} directory not yet created.")
