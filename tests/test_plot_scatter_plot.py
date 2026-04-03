"""Tests for the display-space-thinned scatter plotting helpers."""

from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
from matplotlib.collections import PathCollection

from skyborn.plot import scatter as public_scatter
from skyborn.plot.scatter import _array_scatter as array_scatter


class TestScatterPlot:
    """Tests for the low-level scatter plotting core."""

    def test_array_scatter_grid_mask_uses_default_display_thinning(self):
        x = np.linspace(0.0, 359.0, 72)
        y = np.linspace(-90.0, 90.0, 37)
        mask = np.ones((37, 72), dtype=bool)

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, x, y, where=mask, s=4.0, c="k")

        assert isinstance(result, PathCollection)
        assert 0 < len(result.get_offsets()) < mask.sum()

        plt.close(fig)

    def test_array_scatter_paired_points_preserves_all_points_by_default(self):
        x = np.linspace(0.0, 10.0, 25)
        y = np.sin(x)

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, x, y, s=9.0, c="tab:blue")

        assert isinstance(result, PathCollection)
        assert len(result.get_offsets()) == len(x)

        plt.close(fig)

    def test_array_scatter_paired_points_can_be_thinned_explicitly(self):
        x = np.linspace(0.0, 1.0, 120)
        y = np.zeros_like(x) + 0.5

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, x, y, density=0.2, s=6.0, c="0.2")

        assert isinstance(result, PathCollection)
        assert 0 < len(result.get_offsets()) < len(x)

        plt.close(fig)

    def test_array_scatter_accepts_2d_meshgrid_coordinates(self):
        x1d = np.linspace(100.0, 120.0, 8)
        y1d = np.linspace(10.0, 30.0, 6)
        x2d, y2d = np.meshgrid(x1d, y1d, indexing="xy")
        mask = (x2d + y2d) > 125.0

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, x2d, y2d, where=mask, distance=1e-6)

        assert isinstance(result, PathCollection)
        assert len(result.get_offsets()) == int(mask.sum())

        plt.close(fig)

    def test_array_scatter_subsets_grid_shaped_size_and_color_fields(self):
        x = np.linspace(0.0, 4.0, 5)
        y = np.linspace(10.0, 13.0, 4)
        mask = np.array(
            [
                [True, False, True, False, True],
                [False, True, False, True, False],
                [True, True, False, False, True],
                [False, False, True, True, False],
            ],
            dtype=bool,
        )
        sizes = np.arange(mask.size, dtype=float).reshape(mask.shape) + 5.0
        colors = np.arange(mask.size, dtype=float).reshape(mask.shape)
        linewidths = sizes / 10.0

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(
            ax,
            x,
            y,
            where=mask,
            s=sizes,
            c=colors,
            linewidths=linewidths,
            cmap="viridis",
            distance=1e-6,
        )

        np.testing.assert_allclose(result.get_sizes(), sizes[mask])
        np.testing.assert_allclose(result.get_array(), colors[mask])
        np.testing.assert_allclose(result.get_linewidths(), linewidths[mask])

        plt.close(fig)

    def test_array_scatter_transposes_dataarray_mask_to_xy_dims(self):
        lon = xr.DataArray(
            np.linspace(100.0, 110.0, 4),
            dims=("lon",),
            coords={"lon": np.linspace(100.0, 110.0, 4)},
        )
        lat = xr.DataArray(
            np.linspace(20.0, 26.0, 3),
            dims=("lat",),
            coords={"lat": np.linspace(20.0, 26.0, 3)},
        )
        mask = xr.DataArray(
            np.array(
                [
                    [False, True, False],
                    [True, False, False],
                    [False, False, True],
                    [False, True, False],
                ],
                dtype=bool,
            ),
            dims=("lon", "lat"),
            coords={"lon": lon, "lat": lat},
        )

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, lon, lat, where=mask, distance=1e-6)

        expected_lon, expected_lat = np.meshgrid(lon.data, lat.data, indexing="xy")
        expected = np.column_stack(
            [
                expected_lon[mask.transpose("lat", "lon")],
                expected_lat[mask.transpose("lat", "lon")],
            ]
        )
        offsets = np.asarray(result.get_offsets(), dtype=float)

        np.testing.assert_allclose(offsets, expected)

        plt.close(fig)

    def test_array_scatter_rejects_where_and_mask_together(self):
        x = np.linspace(0.0, 1.0, 4)
        y = np.linspace(0.0, 1.0, 3)
        mask = np.ones((3, 4), dtype=bool)

        fig, ax = plt.subplots(figsize=(6, 4))
        with pytest.raises(ValueError, match="Use only one of 'where' or 'mask'"):
            array_scatter(ax, x, y, where=mask, mask=mask)
        plt.close(fig)

    def test_array_scatter_rejects_selection_shape_mismatch(self):
        x = np.linspace(0.0, 4.0, 5)
        y = np.linspace(0.0, 3.0, 4)
        bad_mask = np.ones((5, 4), dtype=bool)

        fig, ax = plt.subplots(figsize=(6, 4))
        with pytest.raises(ValueError, match="must match candidate shape"):
            array_scatter(ax, x, y, where=bad_mask)
        plt.close(fig)

    def test_array_scatter_accepts_min_distance_as_backward_compatible_alias(self):
        x = np.linspace(0.0, 1.0, 20)
        y = np.linspace(0.0, 1.0, 10)
        mask = np.ones((10, 20), dtype=bool)

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, x, y, where=mask, min_distance=1e-6)

        assert isinstance(result, PathCollection)
        assert len(result.get_offsets()) == int(mask.sum())

        plt.close(fig)

    def test_array_scatter_rejects_distance_and_min_distance_together(self):
        x = np.linspace(0.0, 1.0, 4)
        y = np.linspace(0.0, 1.0, 3)
        mask = np.ones((3, 4), dtype=bool)

        fig, ax = plt.subplots(figsize=(6, 4))
        with pytest.raises(
            TypeError, match="Use only one of 'distance' or 'min_distance'"
        ):
            array_scatter(ax, x, y, where=mask, distance=0.02, min_distance=0.01)
        plt.close(fig)

    def test_array_scatter_returns_empty_collection_for_empty_selection(self):
        x = np.linspace(0.0, 4.0, 5)
        y = np.linspace(0.0, 3.0, 4)
        mask = np.zeros((4, 5), dtype=bool)

        fig, ax = plt.subplots(figsize=(6, 4))
        result = array_scatter(ax, x, y, where=mask, c="k")

        assert isinstance(result, PathCollection)
        assert len(result.get_offsets()) == 0

        plt.close(fig)


class TestScatterWrapper:
    """Tests for the public scatter wrapper."""

    def test_public_scatter_is_exported(self):
        assert callable(public_scatter)

    def test_public_scatter_uses_current_axes_when_ax_is_omitted(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        x = np.linspace(0.0, 1.0, 5)
        y = np.linspace(1.0, 2.0, 5)

        result = public_scatter(x, y, c="tab:red")

        assert isinstance(result, PathCollection)
        assert result.axes is ax

        plt.close(fig)

    @patch("skyborn.plot.scatter._array_scatter")
    def test_public_scatter_converts_crs_like_transform(self, mock_array_scatter):
        fig, ax = plt.subplots(figsize=(6, 4))
        mock_result = Mock(spec=PathCollection)
        mock_array_scatter.return_value = mock_result
        converted_transform = Mock()

        class DummyCRS:
            def _as_mpl_transform(self, target_ax):
                assert target_ax is ax
                return converted_transform

        public_scatter(
            ax,
            np.linspace(0.0, 1.0, 4),
            np.linspace(0.0, 1.0, 4),
            transform=DummyCRS(),
        )

        assert mock_array_scatter.call_args.kwargs["transform"] is converted_transform
        plt.close(fig)

    @patch("skyborn.plot.scatter._array_scatter")
    def test_public_scatter_rejects_duplicate_axes(self, mock_array_scatter):
        fig, ax = plt.subplots(figsize=(6, 4))

        with pytest.raises(TypeError, match="received Axes both positionally"):
            public_scatter(ax, [0, 1], [0, 1], ax=ax)

        mock_array_scatter.assert_not_called()
        plt.close(fig)
