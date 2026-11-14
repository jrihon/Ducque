from pydantic import ValidationError


def format_errors(err: ValidationError):
    """
    for pydantic.ValidationError:

    e.error_count()   # int, the number of errors
    e.errors()        # list[dict], one dict per error
    e.title           # "SugarValidator", the model name
    str(e)            # Pydantic's own multi-line summary
    e.json()          # the errors as a JSON string

    What each error dict contains

    {
        "type":  "enum",                                # machine-readable error code
        "loc":   ("moiety_type",),                      # path to the field (a tuple)
        "msg":   "Input should be 'sugar', 'linker' or 'nucleobase'",
        "input": "sugr",                                # the value that failed
        "ctx":   {"expected": "'sugar', 'linker' or 'nucleobase'"},   # only on some errors
        "url":   "https://errors.pydantic.dev/2.x/v/enum",
    }
    """


class DucqueLibraryNotFound(FileNotFoundError):
    """More specific error to essentially the same problem."""
