"""Skyborn Tropical Cyclone Ventilation Index Module.

This module implements the ventilated potential intensity (vPI) framework
and the associated genesis potential index (GPIv) from Chavas, Camargo &
Tippett (2025, J. Climate).

The module provides:

- Vertical wind shear (VWS) between 200 and 850 hPa
- Entropy deficit (Chi) at 600 hPa
- Ventilation index (VI = VWS * Chi / PI)
- Ventilated potential intensity (vPI) via analytic cubic solution
- Clipped 850-hPa absolute vorticity (eta_c)
- Ventilated genesis potential index (GPIv)

References
----------
Chavas, D. R., S. J. Camargo, and M. K. Tippett, 2025: Tropical Cyclone
Genesis Potential Using a Ventilated Potential Intensity. J. Climate, 38,
1667–1689, https://doi.org/10.1175/JCLI-D-24-0186.1

Tang, B., and K. Emanuel, 2012: A ventilation index for tropical cyclones.
Bull. Amer. Meteor. Soc., 93, 1901–1912.

Komacek, T. D., D. R. Chavas, and D. S. Abbot, 2020: Hurricane genesis is
favorable on terrestrial exoplanets orbiting late-type M dwarf stars.
Astrophys. J., 898, 115.
"""

from .ventilation import (
    absolute_vorticity_850,
    entropy_deficit,
    genesis_potential_index,
    ventilated_pi,
    ventilation_index,
    vertical_wind_shear,
)

__version__ = "1.0.0"

__all__ = [
    "absolute_vorticity_850",
    "entropy_deficit",
    "genesis_potential_index",
    "ventilated_pi",
    "ventilation_index",
    "vertical_wind_shear",
]
