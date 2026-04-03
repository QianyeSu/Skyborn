"""Public Skyborn API for curly-vector plotting.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-03-01 14:58:56
"""

from .ncl_vector import CurlyVectorKey, curly_vector
from .plotting import add_equal_axes, createFigure
from .scatter import scatter
from .streamline import streamline
from .vector_key import curly_vector_key
from .vector_plot import CurlyVectorPlotSet

__all__ = [
    "CurlyVectorKey",
    "CurlyVectorPlotSet",
    "add_equal_axes",
    "createFigure",
    "curly_vector",
    "curly_vector_key",
    "scatter",
    "streamline",
]
