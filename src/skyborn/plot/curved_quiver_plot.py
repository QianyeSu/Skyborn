"""
Functions for curved quiver plots.
"""

from __future__ import annotations
from collections.abc import Hashable
from typing import TYPE_CHECKING

# import matplotlib.font_manager
from .modplot import CurvedQuiverplotSet
from typing import Literal
from matplotlib.patches import FancyArrowPatch

import xarray as xr
import matplotlib.pyplot as plt

# import numpy as np
# import matplotlib
from matplotlib.artist import Artist, allow_rasterization
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.figure import Figure
from matplotlib.backend_bases import RendererBase

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
    if ax is None:
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


class CurvedQuiverLegend(Artist):
    """Curved quiver legend with white background box, arrows scaled according to actual wind speed"""

    def __init__(
        self,
        ax,
        curved_quiver_set,  # Pass the original curved quiver object
        U: float,
        units: str = "m/s",
        width: float = 0.15,
        height: float = 0.08,
        loc: Literal[
            "lower left", "lower right", "upper left", "upper right"
        ] = "lower right",
        # Added labelpos parameter
        labelpos: Literal["N", "S", "E", "W"] = "E",
        max_arrow_length: float = 0.08,  # Maximum arrow length
        arrow_props: dict = None,
        patch_props: dict = None,
        text_props: dict = None,
        padding: float = 0.01,
        margin: float = 0.02,
        reference_speed: float = 2.0,
        # If None, automatically determined based on whether units is empty
        center_label: bool = None,
    ) -> None:
        """
        Initialize curved quiver legend with arrows proportional to wind speed

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to add the legend to
        curved_quiver_set : CurvedQuiverplotSet
            Original curved quiver object, used to get scaling information
        U : float
            Wind speed value represented by the arrow
        units : str
            Unit string, if empty string(""), units won't be displayed and label will be auto-centered
        labelpos : {'N', 'S', 'E', 'W'}, default: 'E'
            Label position relative to arrow:
            'N' - Label above the arrow
            'S' - Label below the arrow
            'E' - Label to the right of the arrow
            'W' - Label to the left of the arrow
        center_label : bool or None
            Whether to center the label. If None, auto-center when units is empty
        max_arrow_length : float
            Maximum arrow length (relative to legend box)
        """
        super().__init__()
        self.reference_speed = reference_speed
        self.margin = margin
        self.ax = ax
        self.curved_quiver_set = curved_quiver_set
        self.U = U
        self.units = units
        self.labelpos = labelpos

        # Automatically determine whether to center label based on units
        self.show_units = units != ""
        if center_label is None:
            self.center_label = not self.show_units  # Auto-center if no units
        else:
            self.center_label = center_label

        self.width = width
        self.height = height
        self.loc = loc
        self.max_arrow_length = max_arrow_length
        self.padding = padding

        # Set default properties
        self.arrow_props = self._setup_arrow_props(arrow_props)
        self.patch_props = self._setup_patch_props(patch_props)
        self.text_props = self._setup_text_props(text_props)

        # Calculate actual arrow length (based on wind speed scale)
        self.arrow_length = self._calculate_arrow_length()

        # Create text content
        if self.center_label:
            self.text_content = f"{U}"
        else:
            self.text_content = f"{U}" if not self.show_units else f"{U} {units}"

        # Create temporary text object to measure size
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10, 10))  # Create temporary figure
        renderer = fig.canvas.get_renderer()

        temp_text = Text(0, 0, self.text_content, **self.text_props)
        temp_text.set_figure(fig)
        bbox = temp_text.get_window_extent(renderer)
        text_width = bbox.width / fig.dpi / 10  # 10 is the width of temporary figure
        text_height = bbox.height / fig.dpi / 10

        plt.close(fig)  # Close temporary figure

        # Adjust box size to ensure it can contain text and arrow
        # Adjust for different labelpos
        if labelpos in ["E", "W"]:  # Horizontal layout
            if self.center_label:
                # Center mode
                min_width = self.arrow_length + 3 * self.padding + text_width
                total_content_width = self.arrow_length + self.padding + text_width
            else:
                # Adjust based on labelpos
                if labelpos == "E":  # Text to the right of the arrow
                    min_width = (
                        2 * self.padding
                        + self.arrow_length
                        + self.padding
                        + text_width
                        + self.padding
                    )
                else:  # labelpos == 'W', text to the left of the arrow
                    min_width = (
                        2 * self.padding
                        + text_width
                        + self.padding
                        + self.arrow_length
                        + self.padding
                    )
                total_content_width = min_width - 2 * self.padding

            # Ensure box is wide enough
            if min_width > self.width:
                self.width = min_width

            # Adjust height to fit text and arrow
            min_height = max(text_height, 0.03) + 2 * self.padding

        else:  # Vertical layout (labelpos in ['N', 'S'])
            # In vertical direction, text and arrow are stacked
            min_height = (
                text_height + 2 * self.padding + 0.03
            )  # 0.03 is an estimate for arrow height

            # Ensure box is wide enough to accommodate the maximum width of arrow and text
            min_width = max(self.arrow_length, text_width) + 2 * self.padding

            if min_width > self.width:
                self.width = min_width

        # Ensure box height is sufficient
        if min_height > self.height:
            self.height = min_height

        # Calculate bottom-left coordinates of legend box based on position
        self._calculate_position()

        # Create background box
        self.patch = Rectangle(
            xy=(self.x, self.y),
            width=self.width,
            height=self.height,
            transform=ax.transAxes,
            **self.patch_props,
        )

        # Set arrow and text positions based on labelpos
        self._position_arrow_and_text(text_width, text_height)

        # Set z-order
        self.set_zorder(10)
        self.patch.set_zorder(9)
        self.arrow.set_zorder(11)
        self.text.set_zorder(11)

        # Add to axes
        ax.add_artist(self.patch)
        ax.add_artist(self.arrow)
        ax.add_artist(self.text)
        ax.add_artist(self)

    def _position_arrow_and_text(self, text_width, text_height):
        """Position arrow and text based on labelpos"""
        box_center_x = self.x + self.width / 2
        box_center_y = self.y + self.height / 2

        # Set default text alignment
        self.text_props.update(
            {"verticalalignment": "center", "horizontalalignment": "center"}
        )

        if self.center_label:
            # Center mode - entire combination (arrow+text) is centered
            if self.labelpos in ["E", "W"]:  # Horizontal layout
                # Calculate total content width
                total_content_width = self.arrow_length + self.padding + text_width

                # Group start position to center the entire element
                group_start_x = box_center_x - (total_content_width / 2)

                if self.labelpos == "E":  # Text to the right of arrow
                    # Set arrow position
                    arrow_start_x = group_start_x
                    arrow_end_x = arrow_start_x + self.arrow_length

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, box_center_y),
                        (arrow_end_x, box_center_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Set text position
                    text_x = arrow_end_x + self.padding
                    self.text_props["horizontalalignment"] = "left"

                else:  # self.labelpos == 'W', text to the left of arrow
                    # Place text first
                    text_x = group_start_x
                    self.text_props["horizontalalignment"] = "right"

                    # Set arrow position
                    arrow_start_x = text_x + text_width + self.padding
                    arrow_end_x = arrow_start_x + self.arrow_length

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, box_center_y),
                        (arrow_end_x, box_center_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                # Create text
                self.text = Text(
                    text_x,
                    box_center_y,
                    self.text_content,
                    transform=self.ax.transAxes,
                    **self.text_props,
                )

            else:  # Vertical layout (self.labelpos in ['N', 'S'])
                if self.labelpos == "N":  # Text above the arrow
                    # Set arrow position (horizontally centered)
                    arrow_start_x = box_center_x - self.arrow_length / 2
                    arrow_end_x = box_center_x + self.arrow_length / 2

                    # Vertical position (lower half)
                    arrow_y = box_center_y - text_height / 2 - self.padding / 2

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, arrow_y),
                        (arrow_end_x, arrow_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Set text position (upper half)
                    text_y = box_center_y + self.padding / 2
                    self.text_props["verticalalignment"] = "bottom"

                else:  # self.labelpos == 'S', text below the arrow
                    # Set arrow position (horizontally centered)
                    arrow_start_x = box_center_x - self.arrow_length / 2
                    arrow_end_x = box_center_x + self.arrow_length / 2

                    # Vertical position (upper half)
                    arrow_y = box_center_y + text_height / 2 + self.padding / 2

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, arrow_y),
                        (arrow_end_x, arrow_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Set text position (lower half)
                    text_y = box_center_y - self.padding / 2
                    self.text_props["verticalalignment"] = "top"

                # Create text (horizontally centered)
                self.text = Text(
                    box_center_x,
                    text_y,
                    self.text_content,
                    transform=self.ax.transAxes,
                    **self.text_props,
                )

        else:  # Non-centered mode
            if self.labelpos in ["E", "W"]:  # Horizontal layout
                if self.labelpos == "E":  # Text to the right of arrow
                    # Layout starting from left
                    arrow_start_x = self.x + self.padding
                    arrow_end_x = arrow_start_x + self.arrow_length

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, box_center_y),
                        (arrow_end_x, box_center_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Create text label
                    text_x = arrow_end_x + self.padding
                    self.text_props["horizontalalignment"] = "left"

                else:  # self.labelpos == 'W', text to the left of arrow
                    # Layout starting from right
                    arrow_end_x = self.x + self.width - self.padding
                    arrow_start_x = arrow_end_x - self.arrow_length

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, box_center_y),
                        (arrow_end_x, box_center_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Create text label
                    text_x = arrow_start_x - self.padding
                    self.text_props["horizontalalignment"] = "right"

                self.text = Text(
                    text_x,
                    box_center_y,
                    self.text_content,
                    transform=self.ax.transAxes,
                    **self.text_props,
                )

            else:  # Vertical layout (self.labelpos in ['N', 'S'])
                # Horizontally center the arrow
                arrow_start_x = self.x + (self.width - self.arrow_length) / 2
                arrow_end_x = arrow_start_x + self.arrow_length

                if self.labelpos == "N":  # Text above the arrow
                    # Arrow at bottom
                    arrow_y = self.y + self.padding + 0.015  # Give the arrow some space

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, arrow_y),
                        (arrow_end_x, arrow_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Text at top
                    text_y = self.y + self.height - self.padding
                    self.text_props["verticalalignment"] = "top"

                else:  # self.labelpos == 'S', text below the arrow
                    # Arrow at top
                    arrow_y = self.y + self.height - self.padding - 0.015

                    # Create arrow
                    self.arrow = FancyArrowPatch(
                        (arrow_start_x, arrow_y),
                        (arrow_end_x, arrow_y),
                        transform=self.ax.transAxes,
                        **self.arrow_props,
                    )

                    # Text at bottom
                    text_y = self.y + self.padding
                    self.text_props["verticalalignment"] = "bottom"

                # Create text (horizontally centered)
                self.text = Text(
                    box_center_x,
                    text_y,
                    self.text_content,
                    transform=self.ax.transAxes,
                    **self.text_props,
                )

    def _calculate_arrow_length(self):
        """Calculate arrow length based on actual wind speed and original data"""
        try:
            # Use the scale information from the original curved_quiver
            if hasattr(self.curved_quiver_set, "resolution") and hasattr(
                self.curved_quiver_set, "magnitude"
            ):
                # resolution = self.curved_quiver_set.resolution
                # magnitude = self.curved_quiver_set.magnitude

                # Use reference speed to scale the arrow
                reference_speed = getattr(
                    self, "reference_speed", 2.0
                )  # Default reference speed is 2.0

                # Calculate the scale factor
                scale_factor = self.U / reference_speed

                # Ensure arrow length is within reasonable range
                arrow_length = min(
                    # Allow arrows up to 4x maximum
                    scale_factor * self.max_arrow_length,
                    self.max_arrow_length * 4,
                )
                arrow_length = max(
                    arrow_length, self.max_arrow_length * 0.2
                )  # Minimum length
                return arrow_length

            # Method 2: If scale information is unavailable, use simple linear scaling
            # Assume common wind speed range is 0-20 m/s
            scale_factor = min(self.U / 20.0, 1.0)
            arrow_length = max(
                scale_factor * self.max_arrow_length, self.max_arrow_length * 0.2
            )
            return arrow_length

        except Exception as e:
            print(f"Warning: Could not calculate proportional arrow length: {e}")
            # Fall back to fixed length
            return self.max_arrow_length * 0.6

    def _calculate_position(self):
        """Calculate legend box position based on loc parameter"""
        margin = getattr(self, "margin", 0.02)  # Use class attribute or default

        if self.loc == "lower left":
            self.x = margin
            self.y = margin
        elif self.loc == "lower right":
            self.x = 1 - self.width - margin
            self.y = margin
        elif self.loc == "upper left":
            self.x = margin
            self.y = 1 - self.height - margin
        elif self.loc == "upper right":
            self.x = 1 - self.width - margin
            self.y = 1 - self.height - margin
        else:
            raise ValueError(
                f"loc must be one of ['lower left', 'lower right', 'upper left', 'upper right'], got {self.loc}"
            )

    def _setup_arrow_props(self, arrow_props):
        """Set default arrow properties"""
        defaults = {
            "arrowstyle": "->",
            "mutation_scale": 20,
            "linewidth": 1.5,
            "color": "black",
        }
        if arrow_props:
            defaults.update(arrow_props)
        return defaults

    def _setup_patch_props(self, patch_props):
        """Set default background box properties"""
        defaults = {
            "linewidth": 1,
            "edgecolor": "black",
            "facecolor": "white",
            "alpha": 0.9,
        }
        if patch_props:
            defaults.update(patch_props)
        return defaults

    def _setup_text_props(self, text_props):
        """Set default text properties"""
        defaults = {
            "fontsize": 10,
            "verticalalignment": "center",
            "horizontalalignment": "left",
            "color": "black",
        }
        if text_props:
            defaults.update(text_props)
        return defaults

    def set_figure(self, fig: Figure) -> None:
        """Set figure object"""
        super().set_figure(fig)
        self.patch.set_figure(fig)
        self.arrow.set_figure(fig)
        self.text.set_figure(fig)

    @allow_rasterization
    def draw(self, renderer: RendererBase) -> None:
        """Draw the legend"""
        if self.get_visible():
            # Ensure all components are drawn
            self.patch.draw(renderer)
            self.arrow.draw(renderer)
            self.text.draw(renderer)
            self.stale = False


def add_curved_quiverkey(
    ax,
    curved_quiver_set,  # Pass the curved_quiver return object
    U: float = 2.0,
    units: str = "m/s",
    loc: str = "lower right",
    labelpos: str = "E",  # label position
    **kwargs,
) -> CurvedQuiverLegend:
    """
    Convenience function: Add proportionally scaled curved quiver legend to axes

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object
    curved_quiver_set : CurvedQuiverplotSet
        Object returned by curved_quiver function
    U : float
        Wind speed value represented by the arrow
    units : str
        Unit
    loc : str
        Position
    labelpos : {'N', 'S', 'E', 'W'}, default: 'E'
        Label position relative to arrow:
        'N' - Label above the arrow
        'S' - Label below the arrow
        'E' - Label to the right of the arrow
        'W' - Label to the left of the arrow
    **kwargs
        Other parameters passed to CurvedQuiverLegend

    Returns
    -------
    CurvedQuiverLegend
        Legend object
    """
    return CurvedQuiverLegend(
        ax, curved_quiver_set, U, units=units, loc=loc, labelpos=labelpos, **kwargs
    )
