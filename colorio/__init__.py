from . import cs, data, diff, illuminants, observers
from .__about__ import __version__
from ._exceptions import ColorioError
from ._rgb_gamut import (
    plot_rgb_slice,
    save_rgb_gamut,
    save_rgb_slice,
    show_rgb_gamut,
    show_rgb_slice,
)
from ._surface_gamut import save_surface_gamut, show_surface_gamut
from ._tools import (
    show_primary_srgb_gradients,
    show_srgb255_gradient,
    show_xy_gamut,
    xy_gamut_mesh,
)
from ._visible_gamut import (
    plot_visible_slice,
    save_visible_gamut,
    show_visible_gamut,
    show_visible_slice,
)

__all__ = [
    "__version__",
    "data",
    "cs",
    "diff",
    "illuminants",
    "observers",
    "ColorioError",
    #
    "show_xy_gamut",
    "show_visible_gamut",
    "xy_gamut_mesh",
    "show_visible_slice",
    "plot_visible_slice",
    "show_rgb_slice",
    "plot_rgb_slice",
    "save_rgb_slice",
    "save_rgb_gamut",
    "show_rgb_gamut",
    "save_cone_gamut",
    "save_visible_gamut",
    "show_srgb255_gradient",
    "show_primary_srgb_gradients",
    "save_surface_gamut",
    "show_surface_gamut",
]
