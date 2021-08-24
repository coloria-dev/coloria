import inspect
from typing import Type

from numpy.typing import ArrayLike

from ..cs import ColorSpace


def create_cs_class_instance(
    cs_class: Type[ColorSpace], whitepoint: ArrayLike, c: float, Y_b: float, L_A: float
):
    spec = inspect.getfullargspec(cs_class)
    kwargs = {}
    if "whitepoint" in spec.args:
        kwargs["whitepoint"] = whitepoint

    if "c" in spec.args:
        kwargs["c"] = c

    if "Y_b" in spec.args:
        kwargs["Y_b"] = Y_b

    if "L_A" in spec.args:
        kwargs["L_A"] = L_A

    return cs_class(**kwargs)
