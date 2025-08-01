"""
@author: Kieran Hunt
https://github.com/kieranmrhunt/curved-quivers/blob/master/modplot.py
https://stackoverflow.com/questions/51843313/flow-visualisation-in-python-using-curved-path-following-vectors
aligned with matplotlib.streamplot and improved by Jelmer Veenstra (30-03-2023), https://github.com/veenstrajelmer
matplotlib.streamplot available at https://raw.githubusercontent.com/matplotlib/matplotlib/main/lib/matplotlib/streamplot.py

Streamline plotting for 2D vector fields.

"""

import numpy as np
import matplotlib as mpl
from matplotlib import cm, patches
import matplotlib.colors as mcolors
import matplotlib.collections as mcollections
import matplotlib.lines as mlines

__all__ = ["velovect"]


def velovect(
    axes,
    x,
    y,
    u,
    v,
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
):
    """
    Draw streamlines of a vector flow.

    Parameters
    ----------
    x, y : 1D/2D arrays
        Evenly spaced strictly increasing arrays to make a grid.  If 2D, all
        rows of *x* must be equal and all columns of *y* must be equal; i.e.,
        they must be as if generated by ``np.meshgrid(x_1d, y_1d)``.
    u, v : 2D arrays
        *x* and *y*-velocities. The number of rows and columns must match
        the length of *y* and *x*, respectively.
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
    minlength : float
        Minimum length of streamline in axes coordinates.
    start_points : (N, 2) array
        Coordinates of starting points for the streamlines in data coordinates
        (the same coordinates as the *x* and *y* arrays).
    zorder : float
        The zorder of the streamlines and arrows.
        Artists with lower zorder values are drawn first.
    maxlength : float
        Maximum length of streamline in axes coordinates.
    integration_direction : {'forward', 'backward', 'both'}, default: 'both'
        Integrate the streamline in forward, backward or both directions.
    data : indexable object, optional
        DATA_PARAMETER_PLACEHOLDER
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
    """

    grid = Grid(x, y)
    mask = StreamMask(10)
    dmap = DomainMap(grid, mask)

    if zorder is None:
        zorder = mlines.Line2D.zorder

    # default to data coordinates
    if transform is None:
        transform = axes.transData

    if color is None:
        color = axes._get_lines.get_next_color()

    if linewidth is None:
        linewidth = mpl.rcParams["lines.linewidth"]

    line_kw = {}
    arrow_kw = dict(arrowstyle=arrowstyle, mutation_scale=10 * arrowsize)

    mpl._api.check_in_list(
        ["both", "forward", "backward"], integration_direction=integration_direction
    )

    # if integration_direction == 'both': #done for magnitude
    #    maxlength /= 2.

    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        if color.shape != grid.shape:
            raise ValueError(
                "If 'color' is given, it must match the shape of " "the (x, y) grid"
            )
        line_colors = [[]]  # Empty entry allows concatenation of zero arrays.
        color = np.ma.masked_invalid(color)
    else:
        line_kw["color"] = color
        arrow_kw["color"] = color

    if isinstance(linewidth, np.ndarray):
        if linewidth.shape != grid.shape:
            raise ValueError(
                "If 'linewidth' is given, it must match the " "shape of the (x, y) grid"
            )
        line_kw["linewidth"] = []
    else:
        line_kw["linewidth"] = linewidth
        arrow_kw["linewidth"] = linewidth

    line_kw["zorder"] = zorder
    arrow_kw["zorder"] = zorder

    # Sanity checks.
    if u.shape != grid.shape or v.shape != grid.shape:
        raise ValueError("'u' and 'v' must match the shape of the (x, y) grid")

    u = np.ma.masked_invalid(u)
    v = np.ma.masked_invalid(v)
    magnitude = np.sqrt(u**2 + v**2)
    magnitude /= np.max(magnitude)
    if integration_direction == "both":
        magnitude /= 2.0

    resolution = density / np.max(grains)
    integrate = _get_integrator(
        u, v, dmap, resolution, magnitude, integration_direction
    )

    trajectories = []

    if start_points is None:
        start_points = _gen_starting_points(x, y, grains)

    sp2 = np.asanyarray(start_points, dtype=float).copy()

    # Check if start_points are outside the data boundaries
    for xs, ys in sp2:
        if not (
            grid.x_origin <= xs <= grid.x_origin + grid.width
            and grid.y_origin <= ys <= grid.y_origin + grid.height
        ):
            raise ValueError(
                f"Starting point ({xs}, {ys}) outside of " "data boundaries"
            )

    # Convert start_points from data to array coords
    # Shift the seed points from the bottom left of the data so that
    # data2grid works properly.
    sp2[:, 0] -= grid.x_origin
    sp2[:, 1] -= grid.y_origin

    for xs, ys in sp2:
        xg, yg = dmap.data2grid(xs, ys)
        # Floating point issues can cause xg, yg to be slightly out of
        # bounds for xs, ys on the upper boundaries. Because we have
        # already checked that the starting points are within the original
        # grid, clip the xg, yg to the grid to work around this issue
        xg = np.clip(xg, 0, grid.nx - 1)
        yg = np.clip(yg, 0, grid.ny - 1)

        t = integrate(xg, yg, broken_streamlines)
        if t is not None:
            trajectories.append(t)

    if use_multicolor_lines:
        if norm is None:
            norm = mcolors.Normalize(color.min(), color.max())
        cmap = cm._ensure_cmap(cmap)

    streamlines = []
    arrows = []
    for t in trajectories:
        tgx, tgy = t.T
        # Rescale from grid-coordinates to data-coordinates.
        tx, ty = dmap.grid2data(tgx, tgy)
        tx += grid.x_origin
        ty += grid.y_origin

        # Create multiple tiny segments if varying width or color is given
        if isinstance(linewidth, np.ndarray) or use_multicolor_lines:
            points = np.transpose([tx, ty]).reshape(-1, 1, 2)
            streamlines.extend(np.hstack([points[:-1], points[1:]]))
        else:
            points = np.transpose([tx, ty])
            streamlines.append(points)

        # Add arrows halfway along each trajectory.
        s = np.cumsum(np.hypot(np.diff(tx), np.diff(ty)))
        n = np.searchsorted(s, s[-1])
        arrow_tail = (tx[n], ty[n])
        arrow_head = (np.mean(tx[n : n + 2]), np.mean(ty[n : n + 2]))

        if isinstance(linewidth, np.ndarray):
            line_widths = interpgrid(linewidth, tgx, tgy)[:-1]
            line_kw["linewidth"].extend(line_widths)
            arrow_kw["linewidth"] = line_widths[n]

        if use_multicolor_lines:
            color_values = interpgrid(color, tgx, tgy)[:-1]
            line_colors.append(color_values)
            arrow_kw["color"] = cmap(norm(color_values[n]))

        p = patches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform, **arrow_kw
        )
        arrows.append(p)

    lc = mcollections.LineCollection(streamlines, transform=transform, **line_kw)
    lc.sticky_edges.x[:] = [grid.x_origin, grid.x_origin + grid.width]
    lc.sticky_edges.y[:] = [grid.y_origin, grid.y_origin + grid.height]
    if use_multicolor_lines:
        lc.set_array(np.ma.hstack(line_colors))
        lc.set_cmap(cmap)
        lc.set_norm(norm)
    axes.add_collection(lc)

    ac = mcollections.PatchCollection(arrows)
    # Adding the collection itself is broken; see #2341.
    for p in arrows:
        axes.add_patch(p)

    axes.autoscale_view()
    stream_container = CurvedQuiverplotSet(
        lc,
        ac,
        resolution,
        magnitude,
        zorder,
        transform,
        axes,
        linewidth,
        color,
        cmap,
        arrowsize,
        arrowstyle,
        start_points,
        integration_direction,
        grains,
        broken_streamlines,
    )
    return stream_container


class CurvedQuiverplotSet:
    def __init__(
        self,
        lines,
        arrows,
        resolution,
        magnitude,
        zorder,
        transform,
        axes,
        linewidth,
        color,
        cmap,
        arrowsize,
        arrowstyle,
        start_points,
        integration_direction,
        grains,
        broken_streamlines,
    ):
        self.lines = lines
        self.arrows = arrows
        self.resolution = resolution
        self.magnitude = magnitude
        self.zorder = zorder
        self.transform = transform
        self.axes = axes
        self.linewidth = linewidth
        self.color = color
        self.cmap = cmap
        self.arrowsize = arrowsize
        self.arrowstyle = arrowstyle
        self.start_points = start_points
        self.integration_direction = integration_direction
        self.grains = grains
        self.broken_streamlines = broken_streamlines

        # New attribute added for handling scaling
        self.max_magnitude = np.max(magnitude) if magnitude is not None else None

        # Set scaling factor based on integration_direction
        self.scale_factor = 2.0 if integration_direction == "both" else 1.0

    def get_scale_factor(self):
        """Get the current scaling factor"""
        return self.scale_factor

    def scale_value(self, value):
        """Convert physical value to display value (for conversion from actual wind speed to display scale)"""
        return value / self.scale_factor

    def unscale_value(self, value):
        """Convert display value back to physical value (for conversion from display scale to actual wind speed)"""
        return value * self.scale_factor


# Coordinate definitions
# ========================


class DomainMap:
    """
    Map representing different coordinate systems.

    Coordinate definitions:

    * axes-coordinates goes from 0 to 1 in the domain.
    * data-coordinates are specified by the input x-y coordinates.
    * grid-coordinates goes from 0 to N and 0 to M for an N x M grid,
      where N and M match the shape of the input data.
    * mask-coordinates goes from 0 to N and 0 to M for an N x M mask,
      where N and M are user-specified to control the density of streamlines.

    This class also has methods for adding trajectories to the StreamMask.
    Before adding a trajectory, run `start_trajectory` to keep track of regions
    crossed by a given trajectory. Later, if you decide the trajectory is bad
    (e.g., if the trajectory is very short) just call `undo_trajectory`.
    """

    def __init__(self, grid, mask):
        self.grid = grid
        self.mask = mask
        # Constants for conversion between grid- and mask-coordinates
        self.x_grid2mask = (mask.nx - 1) / (grid.nx - 1)
        self.y_grid2mask = (mask.ny - 1) / (grid.ny - 1)

        self.x_mask2grid = 1.0 / self.x_grid2mask
        self.y_mask2grid = 1.0 / self.y_grid2mask

        self.x_data2grid = 1.0 / grid.dx
        self.y_data2grid = 1.0 / grid.dy

    def grid2mask(self, xi, yi):
        """Return nearest space in mask-coords from given grid-coords."""
        return round(xi * self.x_grid2mask), round(yi * self.y_grid2mask)

    def mask2grid(self, xm, ym):
        return xm * self.x_mask2grid, ym * self.y_mask2grid

    def data2grid(self, xd, yd):
        return xd * self.x_data2grid, yd * self.y_data2grid

    def grid2data(self, xg, yg):
        return xg / self.x_data2grid, yg / self.y_data2grid

    def start_trajectory(self, xg, yg, broken_streamlines=True):
        xm, ym = self.grid2mask(xg, yg)
        self.mask._start_trajectory(xm, ym, broken_streamlines)

    def reset_start_point(self, xg, yg):
        xm, ym = self.grid2mask(xg, yg)
        self.mask._current_xy = (xm, ym)

    def update_trajectory(self, xg, yg, broken_streamlines=True):
        if not self.grid.within_grid(xg, yg):
            raise InvalidIndexError
        xm, ym = self.grid2mask(xg, yg)
        self.mask._update_trajectory(xm, ym, broken_streamlines)

    def undo_trajectory(self):
        self.mask._undo_trajectory()


class Grid:
    """Grid of data."""

    def __init__(self, x, y):

        if np.ndim(x) == 1:
            pass
        elif np.ndim(x) == 2:
            x_row = x[0]
            if not np.allclose(x_row, x):
                raise ValueError("The rows of 'x' must be equal")
            x = x_row
        else:
            raise ValueError("'x' can have at maximum 2 dimensions")

        if np.ndim(y) == 1:
            pass
        elif np.ndim(y) == 2:
            yt = np.transpose(y)  # Also works for nested lists.
            y_col = yt[0]
            if not np.allclose(y_col, yt):
                raise ValueError("The columns of 'y' must be equal")
            y = y_col
        else:
            raise ValueError("'y' can have at maximum 2 dimensions")

        if not (np.diff(x) > 0).all():
            raise ValueError("'x' must be strictly increasing")
        if not (np.diff(y) > 0).all():
            raise ValueError("'y' must be strictly increasing")

        self.nx = len(x)
        self.ny = len(y)

        self.dx = x[1] - x[0]
        self.dy = y[1] - y[0]

        self.x_origin = x[0]
        self.y_origin = y[0]

        self.width = x[-1] - x[0]
        self.height = y[-1] - y[0]

        if not np.allclose(np.diff(x), self.width / (self.nx - 1)):
            raise ValueError("'x' values must be equally spaced")
        if not np.allclose(np.diff(y), self.height / (self.ny - 1)):
            raise ValueError("'y' values must be equally spaced")

    @property
    def shape(self):
        return self.ny, self.nx

    def within_grid(self, xi, yi):
        """Return whether (*xi*, *yi*) is a valid index of the grid."""
        # Note that xi/yi can be floats; so, for example, we can't simply check
        # `xi < self.nx` since *xi* can be `self.nx - 1 < xi < self.nx`
        return 0 <= xi <= self.nx - 1 and 0 <= yi <= self.ny - 1


class StreamMask:
    """
    Mask to keep track of discrete regions crossed by streamlines.

    The resolution of this grid determines the approximate spacing between
    trajectories. Streamlines are only allowed to pass through zeroed cells:
    When a streamline enters a cell, that cell is set to 1, and no new
    streamlines are allowed to enter.
    """

    def __init__(self, density):
        try:
            self.nx, self.ny = (30 * np.broadcast_to(density, 2)).astype(int)
        except ValueError as err:
            raise ValueError("'density' must be a scalar or be of length " "2") from err
        if self.nx < 0 or self.ny < 0:
            raise ValueError("'density' must be positive")
        self._mask = np.zeros((self.ny, self.nx))
        self.shape = self._mask.shape

        self._current_xy = None

    def __getitem__(self, args):
        return self._mask[args]

    def _start_trajectory(self, xm, ym, broken_streamlines=True):
        """Start recording streamline trajectory"""
        self._traj = []
        self._update_trajectory(xm, ym, broken_streamlines)

    def _undo_trajectory(self):
        """Remove current trajectory from mask"""
        for t in self._traj:
            self._mask[t] = 0

    def _update_trajectory(self, xm, ym, broken_streamlines=True):
        """
        Update current trajectory position in mask.

        If the new position has already been filled, raise `InvalidIndexError`.
        """
        if self._current_xy != (xm, ym):
            if self[ym, xm] == 0:
                self._traj.append((ym, xm))
                self._mask[ym, xm] = 1
                self._current_xy = (xm, ym)
            else:
                if broken_streamlines:
                    raise InvalidIndexError
                else:
                    pass


class InvalidIndexError(Exception):
    pass


class TerminateTrajectory(Exception):
    pass


# Integrator definitions
# =======================


def _get_integrator(u, v, dmap, resolution, magnitude, integration_direction):

    # rescale velocity onto grid-coordinates for integrations.
    u, v = dmap.data2grid(u, v)

    # speed (path length) will be in axes-coordinates
    u_ax = u / (dmap.grid.nx - 1)
    v_ax = v / (dmap.grid.ny - 1)
    speed = np.ma.sqrt(u_ax**2 + v_ax**2)

    def forward_time(xi, yi):
        if not dmap.grid.within_grid(xi, yi):
            raise OutOfBounds
        ds_dt = interpgrid(speed, xi, yi)
        if ds_dt == 0:
            raise TerminateTrajectory()
        dt_ds = 1.0 / ds_dt
        ui = interpgrid(u, xi, yi)
        vi = interpgrid(v, xi, yi)
        return ui * dt_ds, vi * dt_ds

    def backward_time(xi, yi):
        dxi, dyi = forward_time(xi, yi)
        return -dxi, -dyi

    def integrate(x0, y0, broken_streamlines=True):
        """
        Return x, y grid-coordinates of trajectory based on starting point.

        Integrate both forward and backward in time from starting point in
        grid coordinates.

        Integration is terminated when a trajectory reaches a domain boundary
        or when it crosses into an already occupied cell in the StreamMask. The
        resulting trajectory is None if it is shorter than `minlength`.
        """

        stotal, xy_traj = 0.0, []

        try:
            dmap.start_trajectory(x0, y0, broken_streamlines)
        except InvalidIndexError:
            return None
        if integration_direction in ["both", "backward"]:
            s, xyt = _integrate_rk12(
                x0, y0, dmap, backward_time, resolution, magnitude, broken_streamlines
            )
            stotal += s
            xy_traj += xyt[::-1]

        if integration_direction in ["both", "forward"]:
            dmap.reset_start_point(x0, y0)
            s, xyt = _integrate_rk12(
                x0, y0, dmap, forward_time, resolution, magnitude, broken_streamlines
            )
            stotal += s
            xy_traj += xyt[1:]

        if len(xy_traj) > 1:
            return np.broadcast_arrays(xy_traj, np.empty((1, 2)))[0]
        else:  # reject short trajectories
            dmap.undo_trajectory()
            return None

    return integrate


class OutOfBounds(IndexError):
    pass


def _integrate_rk12(x0, y0, dmap, f, resolution, magnitude, broken_streamlines=True):
    """
    2nd-order Runge-Kutta algorithm with adaptive step size.

    This method is also referred to as the improved Euler's method, or Heun's
    method. This method is favored over higher-order methods because:

    1. To get decent looking trajectories and to sample every mask cell
       on the trajectory we need a small timestep, so a lower order
       solver doesn't hurt us unless the data is *very* high resolution.
       In fact, for cases where the user inputs
       data smaller or of similar grid size to the mask grid, the higher
       order corrections are negligible because of the very fast linear
       interpolation used in `interpgrid`.

    2. For high resolution input data (i.e. beyond the mask
       resolution), we must reduce the timestep. Therefore, an adaptive
       timestep is more suited to the problem as this would be very hard
       to judge automatically otherwise.

    This integrator is about 1.5 - 2x as fast as RK4 and RK45 solvers (using
    similar Python implementations) in most setups.
    """
    # This error is below that needed to match the RK4 integrator. It
    # is set for visual reasons -- too low and corners start
    # appearing ugly and jagged. Can be tuned.
    maxerror = 0.003

    # This limit is important (for all integrators) to avoid the
    # trajectory skipping some mask cells. We could relax this
    # condition if we use the code which is commented out below to
    # increment the location gradually. However, due to the efficient
    # nature of the interpolation, this doesn't boost speed by much
    # for quite a bit of complexity.
    maxds = min(1.0 / dmap.mask.nx, 1.0 / dmap.mask.ny, 0.1)

    ds = maxds
    stotal = 0
    xi = x0
    yi = y0
    xyf_traj = []
    m_total = []

    while True:
        try:
            if dmap.grid.within_grid(xi, yi):
                xyf_traj.append((xi, yi))
                m_total.append(interpgrid(magnitude, xi, yi))
                maxlength = resolution * np.mean(m_total)
            else:
                raise OutOfBounds

            # Compute the two intermediate gradients.
            # f should raise OutOfBounds if the locations given are
            # outside the grid.
            k1x, k1y = f(xi, yi)
            k2x, k2y = f(xi + ds * k1x, yi + ds * k1y)

        except OutOfBounds:
            # Out of the domain during this step.
            # Take an Euler step to the boundary to improve neatness
            # unless the trajectory is currently empty.
            if xyf_traj:
                ds, xyf_traj = _euler_step(xyf_traj, dmap, f)
                stotal += ds
            break
        except TerminateTrajectory:
            break

        dx1 = ds * k1x
        dy1 = ds * k1y
        dx2 = ds * 0.5 * (k1x + k2x)
        dy2 = ds * 0.5 * (k1y + k2y)

        ny, nx = dmap.grid.shape
        # Error is normalized to the axes coordinates
        error = np.hypot((dx2 - dx1) / (nx - 1), (dy2 - dy1) / (ny - 1))

        # Only save step if within error tolerance
        if error < maxerror:
            xi += dx2
            yi += dy2
            try:
                dmap.update_trajectory(xi, yi, broken_streamlines)
            except InvalidIndexError:
                break
            if stotal + ds > maxlength:
                break
            stotal += ds

        # recalculate stepsize based on step error
        if error == 0:
            ds = maxds
        else:
            ds = min(maxds, 0.85 * ds * (maxerror / error) ** 0.5)

    return stotal, xyf_traj


def _euler_step(xyf_traj, dmap, f):
    """Simple Euler integration step that extends streamline to boundary."""
    ny, nx = dmap.grid.shape
    xi, yi = xyf_traj[-1]
    cx, cy = f(xi, yi)
    if cx == 0:
        dsx = np.inf
    elif cx < 0:
        dsx = xi / -cx
    else:
        dsx = (nx - 1 - xi) / cx
    if cy == 0:
        dsy = np.inf
    elif cy < 0:
        dsy = yi / -cy
    else:
        dsy = (ny - 1 - yi) / cy
    ds = min(dsx, dsy)
    xyf_traj.append((xi + cx * ds, yi + cy * ds))
    return ds, xyf_traj


# Utility functions
# ========================


def interpgrid(a, xi, yi):
    """Fast 2D, linear interpolation on an integer grid"""

    Ny, Nx = np.shape(a)
    if isinstance(xi, np.ndarray):
        x = xi.astype(int)
        y = yi.astype(int)
        # Check that xn, yn don't exceed max index
        xn = np.clip(x + 1, 0, Nx - 1)
        yn = np.clip(y + 1, 0, Ny - 1)
    else:
        x = int(xi)
        y = int(yi)
        # conditional is faster than clipping for integers
        if x == (Nx - 1):
            xn = x
        else:
            xn = x + 1
        if y == (Ny - 1):
            yn = y
        else:
            yn = y + 1

    a00 = a[y, x]
    a01 = a[y, xn]
    a10 = a[yn, x]
    a11 = a[yn, xn]
    xt = xi - x
    yt = yi - y
    a0 = a00 * (1 - xt) + a01 * xt
    a1 = a10 * (1 - xt) + a11 * xt
    ai = a0 * (1 - yt) + a1 * yt

    if not isinstance(xi, np.ndarray):
        if np.ma.is_masked(ai):
            raise TerminateTrajectory

    return ai


def _gen_starting_points(x, y, grains):
    if isinstance(grains, tuple):
        nx, ny = grains
    elif isinstance(grains, int):
        nx = ny = grains

    eps = np.finfo(np.float32).eps

    tmp_x = np.linspace(x.min() + eps, x.max() - eps, nx)
    tmp_y = np.linspace(y.min() + eps, y.max() - eps, ny)

    xs = np.tile(tmp_x, ny)
    ys = np.repeat(tmp_y, nx)

    seed_points = np.array([xs, ys])

    return seed_points.T
