import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np

def add_equal_axes(ax, loc, pad, width):
    '''
    在原有的Axes旁新添一个等高或等宽的Axes并返回该对象.

    Parameters
    ----------
    ax : Axes or array_like of Axes
        原有的Axes,也可以为一组Axes构成的数组.

    loc : {'left', 'right', 'bottom', 'top'}
        新Axes相对于旧Axes的位置.

    pad : float
        新Axes与旧Axes的间距.

    width: float
        当loc='left'或'right'时,width表示新Axes的宽度.
        当loc='bottom'或'top'时,width表示新Axes的高度.

    Returns
    -------
    ax_new : Axes
        新Axes对象.
    '''
    # 无论ax是单个还是一组Axes,获取ax的大小位置.
    axes = np.atleast_1d(ax).ravel()
    bbox = mtransforms.Bbox.union([ax.get_position() for ax in axes])

    # 决定新Axes的大小位置.
    if loc == 'left':
        x0_new = bbox.x0 - pad - width
        x1_new = x0_new + width
        y0_new = bbox.y0
        y1_new = bbox.y1
    elif loc == 'right':
        x0_new = bbox.x1 + pad
        x1_new = x0_new + width
        y0_new = bbox.y0
        y1_new = bbox.y1
    elif loc == 'bottom':
        x0_new = bbox.x0
        x1_new = bbox.x1
        y0_new = bbox.y0 - pad - width
        y1_new = y0_new + width
    elif loc == 'top':
        x0_new = bbox.x0
        x1_new = bbox.x1
        y0_new = bbox.y1 + pad
        y1_new = y0_new + width

    # 创建新Axes.
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
