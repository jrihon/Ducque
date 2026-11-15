from pathlib import Path
from typing import Any

import typer
import yaml
from pydantic import ValidationError

from ducque.library.errors import format_errors
from ducque.models.library_models import (
    LinkerValidator,
    MoietyType,
    NucleobaseValidator,
    SugarValidator,
    moiety_type_field,
)

type ModelValidator = SugarValidator | LinkerValidator | NucleobaseValidator


class LibaryValidator:
    """Factory class to validate the input yaml that add gets passed to the library."""

    def __init__(self, yaml_file: Path):

        self.py_obj: dict[str, Any] = yaml.safe_load(yaml_file)
        validator_cls = self._assert_molecule_type()

        try:
            self.model = validator_cls.model_validate(self.py_obj)
        except ValidationError as err:
            typer.secho(format_errors(err), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)

    def _assert_molecule_type(self) -> type[ModelValidator]:
        """Use the `moiety_type_field` to validate the input molecule type.

        This pydantic.TypeAdapter is used to parse the given dict.
        This only checks for the `moiety_type` field, to do an initial check.
        """
        moiety_type = moiety_type_field.validate_python(self.py_obj)

        match moiety_type:
            case MoietyType.SUGAR:
                return SugarValidator
            case MoietyType.LINKER:
                return LinkerValidator
            case MoietyType.NUCLEOBASE:
                return NucleobaseValidator
