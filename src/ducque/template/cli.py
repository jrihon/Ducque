import typer

from ducque.template.template import Template

template_app = typer.Typer(help="Create template input files for ducque.")


@template_app.command()
def template(template_type: str):
    """Return the queried template inputfile to stdout."""
    a = Template(template_type)
    print("return the given template ...")
    print(a)
