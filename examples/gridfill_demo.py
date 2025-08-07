#!/usr/bin/env python
"""
GridFill Integration Demo
========================

This example demonstrates the integrated gridfill functionality in skyborn.
"""

from skyborn.gridfill import fill as gridfill_fill
import numpy as np
import matplotlib.pyplot as plt
import skyborn

# Create some sample data with missing values
np.random.seed(42)
nlat, nlon = 50, 100
lat = np.linspace(-90, 90, nlat)
lon = np.linspace(0, 360, nlon)

# Create a synthetic field with some structure
LON, LAT = np.meshgrid(lon, lat)
data = np.sin(np.radians(3 * LON)) * np.cos(np.radians(2 * LAT))
data = data + 0.2 * np.random.randn(nlat, nlon)

# Create missing values (simulate satellite gaps)
missing_mask = np.zeros_like(data, dtype=bool)
missing_mask[20:30, 40:60] = True  # Large gap
missing_mask[10:15, 80:90] = True  # Small gap

# Create masked array
data_with_gaps = np.ma.array(data, mask=missing_mask)

print("Original data shape:", data.shape)
print("Number of missing values:", np.sum(missing_mask))

# Use skyborn's integrated gridfill to fill the gaps
print("Filling missing values using gridfill...")
filled_data, converged = skyborn.fill(data_with_gaps, xdim=1, ydim=0, eps=1e-4)

print("Fill converged:", converged[0])
print("Filled data shape:", filled_data.shape)

# Calculate error (where we have original data)
valid_mask = ~missing_mask
error = np.abs(filled_data - data)[valid_mask]
print(f"RMS error in valid regions: {np.sqrt(np.mean(error**2)):.6f}")

# Also test the direct module access

filled_data2, converged2 = gridfill_fill(data_with_gaps, xdim=1, ydim=0, eps=1e-4)
print("Direct gridfill access also works:", converged2[0])

print("\nGridfill integration successful!")
print("You can now use both:")
print("  - skyborn.fill() for direct access")
print("  - skyborn.gridfill.fill() for module access")
