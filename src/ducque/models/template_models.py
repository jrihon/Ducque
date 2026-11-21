from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel

from ducque.models.construct_models import ConstructValidator
from ducque.models.library_models import (
    LinkerValidator,
    MoietyType,
    NucleobaseValidator,
    SugarValidator,
)

_base_path_test = Path(".")


class TemplateType(StrEnum):
    """List of possible templates to produce.

    We use the `.model_construct()` method to bypass validation, since the source is trusted.
    """

    ADD_LINKER = "add-linker"
    ADD_NUCLEOBASE = "add-nucleobase"
    ADD_SUGAR = "add-sugar"
    CONSTRUCT = "construct"

    @classmethod
    def _fill_sugar(cls) -> SugarValidator:

        return SugarValidator.model_construct(
            file_path=_base_path_test, moiety_type=MoietyType.SUGAR
        )

    @classmethod
    def _fill_linker(cls) -> LinkerValidator:

        return LinkerValidator.model_construct(
            file_path=_base_path_test, moiety_type=MoietyType.LINKER
        )

    @classmethod
    def _fill_nucleobase(cls) -> NucleobaseValidator:

        return NucleobaseValidator.model_construct(
            file_path=_base_path_test, moiety_type=MoietyType.NUCLEOBASE
        )

    @classmethod
    def _fill_construct(cls) -> ConstructValidator:

        return ConstructValidator.model_construct(
            sequence=["foo"], complement=["bar"], output_filename="sugar"
        )

    def get_model(self) -> BaseModel:

        match self:
            case TemplateType.ADD_SUGAR:
                return TemplateType._fill_sugar()
            case TemplateType.ADD_NUCLEOBASE:
                return TemplateType._fill_nucleobase()
            case TemplateType.ADD_LINKER:
                return TemplateType._fill_linker()
            case TemplateType.CONSTRUCT:
                return TemplateType._fill_construct()
