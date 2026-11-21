import typer

from ducque.template.template import Template

template_app = typer.Typer(help="Create template input files for ducque.")


@template_app.command()
def get(template_type: str):
    """Return the queried template inputfile to stdout."""
    Template(template_type).print_yaml()
