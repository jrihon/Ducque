from enum import StrEnum
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator, TypeAdapter


def _lower(string: str):
    """Transform string before validation occurs."""
    return string.strip().lower()


def _path_exists(path: Path):
    """Transform string before validation occurs."""
    if not path.exists():
        raise FileNotFoundError(path)


##
# The moiety type defines the model to validate with
# validation beforehand to optimise.
##
class MoietyType(StrEnum):
    SUGAR = "sugar"
    LINKER = "linker"
    NUCLEOBASE = "nucleobase"


MoietyTypeField = Annotated[MoietyType, BeforeValidator(_lower)]

moiety_type_field = TypeAdapter(MoietyTypeField)


##
# Moiety BaseModels
##
class MoietyValidatorMixin(BaseModel):
    """Class to use as base to inherit from, for the validator."""

    file_path: Annotated[Path, BeforeValidator(_path_exists)]
    moiety_type: MoietyTypeField


class SugarValidator(MoietyValidatorMixin): ...


class LinkerValidator(MoietyValidatorMixin): ...


class NucleobaseValidator(MoietyValidatorMixin): ...
