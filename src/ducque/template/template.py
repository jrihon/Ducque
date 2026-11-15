from typing import Annotated

from pydantic import BeforeValidator, TypeAdapter

from ducque.models.template_models import TemplateType
from ducque.models.utils import enum_lower

TemplateTypeArgument = Annotated[TemplateType, BeforeValidator(enum_lower)]

template_type_argument = TypeAdapter(TemplateTypeArgument)


class Template:
    """Factory class to decide which template to return"""

    def __init__(self, template_type: str) -> None:
        self.template_type = template_type_argument.validate_strings(template_type)

    def build_template(self) -> dict[str, str]:
        """Build the template with the enum variants."""

    def print(py_obj: dict[str, str]) -> None:
        """The output is yaml-formatted."""
