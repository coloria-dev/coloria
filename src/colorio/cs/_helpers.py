import re

string_to_cs = {}


def _normalize_string(string: str) -> str:
    # remove all alphanumerical characters
    string = re.sub(r"[^a-zA-Z0-9]", "", string)
    string = string.lower()
    return string


def register(string, cs):
    string_to_cs[_normalize_string(string)] = cs
