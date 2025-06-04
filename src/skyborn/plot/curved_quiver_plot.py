"""
Functions for curved quiver plots.
"""

from __future__ import annotations

import warnings
from collections.abc import Hashable
from typing import TYPE_CHECKING

import matplotlib.font_manager
from .modplot import CurvedQuiverplotSet
from typing import Literal
from matplotlib.patches import FancyArrowPatch

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = ["curved_quiver", "add_curved_quiverkey"]


def curved_quiver(
    ds: xr.Dataset,
    x: Hashable,
    y: Hashable,
    u: Hashable,
    v: Hashable,
    ax: Axes | None = None,
    density=1,
    linewidth=None,
    color=None,
    cmap=None,
    norm=None,
    arrowsize=1,
    arrowstyle="-|>",
    transform=None,
    zorder=None,
    start_points=None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
) -> CurvedQuiverplotSet:
    """
    Plot streamlines of a vector flow.

    .. warning::

        This function is experimental and the API is subject to change. Please use with caution.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset`.
        Wind dataset.
    x : Hashable or None, optional.
        Variable name for x-axis.
    y : Hashable or None, optional.
        Variable name for y-axis.
    u : Hashable or None, optional.
        Variable name for the u velocity (in `x` direction).
    v : Hashable or None, optional.
        Variable name for the v velocity (in `y` direction).
    ax : :py:class:`matplotlib.axes.Axes`, optional.
        Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
    density : float or (float, float)
        Controls the closeness of streamlines. When ``density = 1``, the domain
        is divided into a 30x30 grid. *density* linearly scales this grid.
        Each cell in the grid can have, at most, one traversing streamline.
        For different densities in each direction, use a tuple
        (density_x, density_y).
    linewidth : float or 2D array
        The width of the streamlines. With a 2D array the line width can be
        varied across the grid. The array must have the same shape as *u*
        and *v*.
    color : color or 2D array
        The streamline color. If given an array, its values are converted to
        colors using *cmap* and *norm*.  The array must have the same shape
        as *u* and *v*.
    cmap, norm
        Data normalization and colormapping parameters for *color*; only used
        if *color* is an array of floats. See `~.Axes.imshow` for a detailed
        description.
    arrowsize : float
        Scaling factor for the arrow size.
    arrowstyle : str
        Arrow style specification.
        See `~matplotlib.patches.FancyArrowPatch`.
    start_points : (N, 2) array
        Coordinates of starting points for the streamlines in data coordinates
        (the same coordinates as the *x* and *y* arrays).
    zorder : float
        The zorder of the streamlines and arrows.
        Artists with lower zorder values are drawn first.
    integration_direction : {'forward', 'backward', 'both'}, default: 'both'
        Integrate the streamline in forward, backward or both directions.
    broken_streamlines : boolean, default: True
        If False, forces streamlines to continue until they
        leave the plot domain.  If True, they may be terminated if they
        come too close to another streamline.

    Returns
    -------
    CurvedQuiverplotSet
        Container object with attributes

        - ``lines``: `.LineCollection` of streamlines

        - ``arrows``: `.PatchCollection` containing `.FancyArrowPatch`
          objects representing the arrows half-way along streamlines.

            This container will probably change in the future to allow changes
            to the colormap, alpha, etc. for both lines and arrows, but these
            changes should be backward compatible.

    .. seealso::
        - https://github.com/matplotlib/matplotlib/issues/20038
        - https://github.com/kieranmrhunt/curved-quivers
        - https://github.com/Deltares/dfm_tools/issues/483
        - https://github.com/NCAR/geocat-viz/issues/4
        - https://docs.xarray.dev/en/stable/generated/xarray.Dataset.plot.streamplot.html#xarray.Dataset.plot.streamplot
    """
    from .modplot import velovect

    ds = ds.sortby(y)
    x = ds[x].data
    y = ds[y].data
    u = ds[u].data
    v = ds[v].data

    # https://scitools.org.uk/cartopy/docs/latest/gallery/miscellanea/logo.html#sphx-glr-gallery-miscellanea-logo-py
    if ax == None:
        ax = plt.gca()
    if type(transform).__name__ == "PlateCarree":
        transform = transform._as_mpl_transform(ax)

    # https://github.com/Deltares/dfm_tools/issues/294
    # https://github.com/Deltares/dfm_tools/blob/main/dfm_tools/modplot.py
    obj = velovect(
        ax,
        x,
        y,
        u,
        v,
        density=density,
        linewidth=linewidth,
        color=color,
        cmap=cmap,
        norm=norm,
        arrowsize=arrowsize,
        arrowstyle=arrowstyle,
        transform=transform,
        zorder=zorder,
        start_points=start_points,
        integration_direction=integration_direction,
        grains=grains,
        broken_streamlines=broken_streamlines,
    )
    return obj


def add_curved_quiverkey(
    curved_quiver: CurvedQuiverplotSet,
    X: float,
    Y: float,
    U: float,
    label: str,
    ax: matplotlib.axes.Axes = None,
    color: str = "black",
    angle: float = 0.0,
    labelpos: Literal["N", "S", "E", "W"] = "N",
    labelsep: float = 0.02,
    labelcolor: str = None,
    fontproperties: matplotlib.font_manager.FontProperties = None,
    zorder: float = None,
):
    """
    Add a key to a quiver plot.

    The positioning of the key depends on X, Y, coordinates, and labelpos.
    If labelpos is 'N' or 'S', X, Y give the position of the middle of the key arrow.
    If labelpos is 'E', X, Y positions the head, and if labelpos is 'W', X, Y positions the tail;
    in either of these two cases, X, Y is somewhere in the middle of the arrow+label key object.

    .. warning::

        This function is experimental and the API is subject to change. Please use with caution.

    Parameters
    ----------
    Q : :py:class:`easyclimate.modplot.CurvedQuiverplotSet`
        A `.CurvedQuiverplotSet` object as returned by a call to `curved_quiver()`.
    X, Y : float
        The location of the key.
    U : float
        The length of the key.
    label : str
        The key label (e.g., length and units of the key).
    ax : :py:class:`matplotlib.axes.Axes`, optional.
        Axes on which to plot.
    angle : float, default: `0.0`.
        The angle of the key arrow, in degrees anti-clockwise from the
        horizontal axis.
    labelpos : {'N', 'S', 'E', 'W'}, default: `N`.
        Position the label above, below, to the right, to the left of the
        arrow, respectively.
    labelsep : float, default: `0.02`.
        Distance in inches between the arrow and the label.
    labelcolor : str.
            Label color.
    fontproperties : dict, optional
        A dictionary with keyword arguments accepted by the
        `~matplotlib.font_manager.FontProperties` initializer:
        *family*, *style*, *variant*, *size*, *weight*.
    zorder : float
            The zorder of the key.
    """
    if ax == None:
        ax = plt.gca()

    pos = (X, Y)

    # Calculate arrow length in data coordinates
    resolution = curved_quiver.resolution
    magnitude = curved_quiver.magnitude
    arrow_length = U * resolution / np.mean(magnitude)

    # Convert angle to radians
    angle_rad = np.radians(angle)

    # Calculate the components of the arrow based on the angle
    dx = arrow_length * np.cos(angle_rad)
    dy = arrow_length * np.sin(angle_rad)

    # Define arrow properties
    arrow_props = dict(linewidth=0.02)

    # Draw the arrow with the given angle
    # ax.annotate(
    #     "",
    #     xy=(pos[0] + dx, pos[1] + dy),  # Arrow endpoint
    #     xytext=pos,  # Arrow start
    #     arrowprops=arrow_props,
    #     xycoords="axes fraction",
    #     zorder=zorder if zorder is not None else 2,
    # )
    arrow = FancyArrowPatch(
        (pos[0], pos[1]),
        (pos[0] + arrow_length, pos[1]),
        mutation_scale=15,
        arrowstyle="->",
        **arrow_props,
    )
    ax.add_patch(arrow)

    # Calculate the label position based on `labelpos`
    if labelpos == "N":  # Label above the arrow
        label_x = (
            pos[0] + dx / 2 + labelsep * np.cos(angle_rad + np.pi / 2)
        )  # Label at the center, above
        label_y = (
            pos[1] + dy / 2 + labelsep * np.sin(angle_rad + np.pi / 2)
        )  # Label at the center, above
    elif labelpos == "S":  # Label below the arrow (centered)
        label_x = (
            pos[0] + dx / 2 + labelsep * np.cos(angle_rad + np.pi / 2)
        )  # Move in opposite direction
        label_y = (
            pos[1] + dy / 2 - labelsep * np.sin(angle_rad + np.pi / 2)
        )  # Move in opposite direction
    elif labelpos == "E":  # Label to the right of the arrow
        label_x = pos[0] + dx + labelsep * np.cos(angle_rad)
        label_y = pos[1] + dy + labelsep * np.sin(angle_rad)
    elif labelpos == "W":  # Label to the left of the arrow
        label_x = pos[0] + dx - labelsep * np.cos(angle_rad)
        label_y = pos[1] + dy - labelsep * np.sin(angle_rad)

    # Add the label text
    ax.text(
        label_x,
        label_y,
        label,
        color=labelcolor if labelcolor is not None else color,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        fontproperties=fontproperties,
        zorder=zorder if zorder is not None else 2,
    )
