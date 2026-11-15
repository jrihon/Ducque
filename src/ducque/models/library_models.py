from enum import StrEnum
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, BeforeValidator, TypeAdapter

from ducque.models.utils import enum_lower, path_exists


##
# The moiety type defines the model to validate with
# validation beforehand to optimise.
##
class MoietyType(StrEnum):
    SUGAR = "sugar"
    LINKER = "linker"
    NUCLEOBASE = "nucleobase"


MoietyTypeField = Annotated[MoietyType, BeforeValidator(enum_lower)]

moiety_type_field = TypeAdapter(MoietyTypeField)


##
# Moiety BaseModels
##
class MoietyValidatorMixin(BaseModel):
    """Class to use as base to inherit from, for the validator."""

    file_path: Annotated[Path, BeforeValidator(path_exists)]
    moiety_type: MoietyTypeField


class SugarValidator(MoietyValidatorMixin): ...


class LinkerValidator(MoietyValidatorMixin): ...


class NucleobaseValidator(MoietyValidatorMixin): ...
