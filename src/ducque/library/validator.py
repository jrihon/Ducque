from typing import Any

import typer
import yaml
from pydantic import ValidationError

from ducque.library.errors import format_errors
from ducque.models.library_models import (
    LinkerValidator,
    MoietyType,
    MoleculeValidator,
    NucleobaseValidator,
    SugarValidator,
    moiety_type_field,
)


class MoietyTypeValidator:
    """Factory class to validate the type of input, based on the MoietyType.

    A MoleculeValidator is returned.
    """

    def __init__(self, yaml_content: str):

        self.py_obj: dict[str, Any] = yaml.safe_load(yaml_content)
        validator_cls = self._assert_molecule_type()

        try:
            self.model: MoleculeValidator = validator_cls.model_validate(self.py_obj)
        except ValidationError as err:
            typer.secho(format_errors(err), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)

    def _assert_molecule_type(self) -> type[MoleculeValidator]:
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
