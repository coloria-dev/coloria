from __future__ import annotations

import inspect
from typing import Type

import numpy as np
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


def stress_absolute(
    target: ArrayLike, actual: ArrayLike, weights: float | ArrayLike = 1.0
) -> float:
    target = np.asarray(target)
    actual = np.asarray(actual)
    weights = np.asarray(weights)
    wtarget = weights * target
    alpha = np.dot(wtarget, actual) / np.dot(wtarget, target)
    diff = alpha * target - actual
    val = np.dot(weights * diff, diff) / np.dot(weights * actual, actual)
    return 100 * np.sqrt(val)


def stress_relative(
    target: ArrayLike, actual: ArrayLike, weights: float | ArrayLike = 1.0
) -> float:
    target = np.asarray(target)
    actual = np.asarray(actual)
    weights = np.asarray(weights)
    wtarget = weights * target
    alpha = np.sum(wtarget) / np.sum(wtarget * target / actual)
    diff = alpha * target - actual
    val = np.sum(weights * diff**2 / actual) / np.sum(weights * actual)
    return 100 * np.sqrt(val)
