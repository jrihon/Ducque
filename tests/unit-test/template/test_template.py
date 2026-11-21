import pytest

from ducque.models.construct_models import ConstructValidator
from ducque.models.library_models import (
    LinkerValidator,
    MoietyType,
    MoleculeValidator,
    NucleobaseValidator,
    SugarValidator,
)
from ducque.template import Template

# TODO: need to add snapshot testing when *Validators are complete.


@pytest.mark.parametrize(
    "template_type,expected_model,expected_enum",
    [
        ("add-sugar", SugarValidator, MoietyType.SUGAR),
        ("add-linker", LinkerValidator, MoietyType.LINKER),
        ("add-nucleobase", NucleobaseValidator, MoietyType.NUCLEOBASE),
    ],
)
def test_template_factory_moietytype(
    template_type: str,
    expected_model: type[MoleculeValidator],
    expected_enum: MoietyType,
) -> None:

    template = Template(template_type)

    assert isinstance(template.template, expected_model)
    assert template.template.moiety_type is expected_enum


@pytest.mark.parametrize("template_type", [("construct")])
def test_template_factory_construct(template_type: str) -> None:

    template = Template(template_type)

    assert isinstance(template.template, ConstructValidator)
