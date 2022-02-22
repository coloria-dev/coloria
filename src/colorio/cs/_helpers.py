import re

_string_to_cs = {}


def _normalize_string(string: str) -> str:
    # remove all alphanumerical characters
    string = re.sub(r"[^a-zA-Z0-9]", "", string)
    string = string.lower()
    return string


def register(string, cs):
    _string_to_cs[_normalize_string(string)] = cs


def string_to_cs(string):
    try:
        return _string_to_cs[_normalize_string(string)]
    except KeyError:
        keys = ", ".join(_string_to_cs.keys())
        msg = f"Could not find color space `{string}`. Try one of {keys}."
        raise KeyError(msg)
