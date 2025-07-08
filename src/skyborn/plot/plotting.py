import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np


def add_equal_axes(ax, loc, pad, width):
    """
    Add a new Axes with equal height or width next to the original Axes and return that object.

    Parameters
    ----------
    ax : Axes or array_like of Axes
        The original Axes, or can be an array of Axes.

    loc : {'left', 'right', 'bottom', 'top'}
        Position of the new Axes relative to the old Axes.

    pad : float
        Spacing between the new Axes and the old Axes.

    width: float
        When loc='left' or 'right', width represents the width of the new Axes.
        When loc='bottom' or 'top', width represents the height of the new Axes.

    Returns
    -------
    ax_new : Axes
        New Axes object.
    """
    # Whether ax is a single Axes or a group of Axes, get the size and position of ax.
    axes = np.atleast_1d(ax).ravel()
    bbox = mtransforms.Bbox.union([ax.get_position() for ax in axes])

    # Determine the size and position of the new Axes.
    if loc == "left":
        x0_new = bbox.x0 - pad - width
        x1_new = x0_new + width
        y0_new = bbox.y0
        y1_new = bbox.y1
    elif loc == "right":
        x0_new = bbox.x1 + pad
        x1_new = x0_new + width
        y0_new = bbox.y0
        y1_new = bbox.y1
    elif loc == "bottom":
        x0_new = bbox.x0
        x1_new = bbox.x1
        y0_new = bbox.y0 - pad - width
        y1_new = y0_new + width
    elif loc == "top":
        x0_new = bbox.x0
        x1_new = bbox.x1
        y0_new = bbox.y1 + pad
        y1_new = y0_new + width

    # Create new Axes.
    fig = axes[0].get_figure()
    bbox_new = mtransforms.Bbox.from_extents(x0_new, y0_new, x1_new, y1_new)
    ax_new = fig.add_axes(bbox_new)

    return ax_new


def createFigure(figsize=(12, 8), dpi=300, subplotAdj=None, **kwargs):
    figsize = figsize
    figure = plt.figure(figsize=figsize, dpi=dpi, **kwargs)
    if subplotAdj is not None:
        plt.subplots_adjust(**subplotAdj)
    return figure
