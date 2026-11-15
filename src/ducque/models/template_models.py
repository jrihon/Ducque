from enum import StrEnum

from pydantic import BaseModel

from ducque.models.construct_models import ConstructValidator
from ducque.models.library_models import (
    LinkerValidator,
    NucleobaseValidator,
    SugarValidator,
)


class TemplateType(StrEnum):
    """List of possible templates to produce."""

    ADD_LINKER = "add-linker"
    ADD_NUCLEOBASE = "add-nucleobase"
    ADD_SUGAR = "add-sugar"
    CONSTRUCT = "construct"

    @classmethod
    def _fill_sugar(cls) -> SugarValidator:

        return SugarValidator(file_path="foo", moiety_type="sugar")

    @classmethod
    def _fill_linker(cls) -> LinkerValidator:

        return LinkerValidator(file_path="foo", moiety_type="sugar")

    @classmethod
    def _fill_nucleobase(cls) -> NucleobaseValidator:

        return NucleobaseValidator(file_path="foo", moiety_type="sugar")

    @classmethod
    def _fill_construct(cls) -> ConstructValidator:

        return ConstructValidator(
            sequence=["foo"], complement=["bar"], output_filename="sugar"
        )

    @classmethod
    def get_model(cls) -> BaseModel:

        match cls:
            case cls.ADD_SUGAR:
                return cls._fill_sugar()
            case cls.ADD_NUCLEOBASE:
                return cls._fill_linker()
            case cls.ADD_LINKER:
                return cls._fill_nucleobase()
            case cls.CONSTRUCT:
                return cls._fill_construct()
            case _:
                raise TypeError
