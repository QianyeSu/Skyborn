# CDO remap vendor subset

Imported from `cdo`.

Included upstream scope:
- `remap_conserv*` for conservative remapping
- `remap_bilinear*` for linear / bilinear remapping
- `remap_knn*` for nearest-neighbor style remapping
- direct shared remap/search helpers required by those three methods
- required `mpim_grid`, `lib/healpix`, and `lib/yac/src` support files
- selected direct helper files such as `knndata.cc`, `interpol.h`, and the
  `lib/yac/src/compute_weights*` / `yac_lapack_interface.h` pieces required by
  the remap subset

Local narrowing notes:
- `src/rect_grid_search.cc` and `src/rect_grid_search.h` are a small Skyborn
  extraction of the upstream `src/interpol.cc` regular-grid search helpers.
- This keeps the vendor subtree self-contained for regular-grid point/cell
  lookup without pulling in the broader CDO runtime from `interpol.cc`.

This is intentionally not a full CDO mirror.
