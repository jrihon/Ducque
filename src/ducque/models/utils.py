from pathlib import Path


def enum_lower(string: str):
    """Transform string before validation occurs."""
    return string.strip().lower()


def path_exists(path: Path):
    """Transform string before validation occurs."""
    if not path.exists():
        raise FileNotFoundError(path)
    return path
