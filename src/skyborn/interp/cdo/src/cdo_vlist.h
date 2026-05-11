#ifndef CDO_VLIST_H
#define CDO_VLIST_H

#include <string>
#include <cdi.h>

#include "cdo_options.h"
#include "varray.h"
#include "cdo_varlist.h"

void vlist_define_timestep_type(int vlistID, int operfunc);

int vlist_compare_x(int vlistID1, int vlistID2, int cmpFlag);
bool vlist_is_szipped(int vlistID);

size_t vlist_check_gridsize(int vlistID);
Varray<double> vlist_read_vct(int vlistID, int &zaxisID_ML, int &numHybridLevels, int &numFullLevels, int &numHalfLevels);
void vlist_change_hybrid_zaxis(int vlistID1, int vlistID2, int zaxisID1, int zaxisID2);

void cdo_compare_grids(int gridID1, int gridID2);  // TODO: Check if this belongs here or if it should be in griddes.*

int vlist_get_first_gaussian_grid(int vlistID);
int vlist_get_first_spectral_grid(int vlistID);
int vlist_get_first_fourier_grid(int vlistID);

void cdo_check_missval(double missval, std::string const &varname);

void vlist_unpack(int vlistID);

int zaxis_check_levels(int zaxisID1, int zaxisID2);

#endif
