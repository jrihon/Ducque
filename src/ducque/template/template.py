from typing import Annotated

from pydantic import BeforeValidator, TypeAdapter
from yaml import safe_dump

from ducque.models.construct_models import ConstructValidator
from ducque.models.library_models import MoleculeValidator
from ducque.models.template_models import TemplateType
from ducque.models.utils import enum_lower

TemplateTypeArgument = Annotated[TemplateType, BeforeValidator(enum_lower)]

template_type_argument = TypeAdapter(TemplateTypeArgument)


class Template:
    """Factory class to decide which template to return"""

    def __init__(self, template_type: str) -> None:
        self.template_type = template_type_argument.validate_strings(template_type)
        self.template = self._build_template()

    def _build_template(self) -> MoleculeValidator | ConstructValidator:
        """Build the template with the enum variants."""
        return self.template_type.get_model()

    def print_yaml(self) -> None:
        """The output is yaml-formatted."""

        # Difference between
        # model_dump(mode="json") output a py_dict for json
        # model_dump_json() outputs a JSON string
        obj = self.template.model_dump(mode="json")

        yaml_string = safe_dump(obj, sort_keys=False)
        print(yaml_string, end="")
