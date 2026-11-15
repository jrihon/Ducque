from pathlib import Path

from pydantic import BaseModel


class ConstructValidator(BaseModel):
    """Class to use as base to inherit from, for the validator."""

    sequence: list[str]
    complement: str | list[str]
    output_filename: Path
