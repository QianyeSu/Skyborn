#ifndef CDO_CDI_WRAPPER_H
#define CDO_CDI_WRAPPER_H

#include <string>

#include "mpim_grid/grid_healpix.h"

namespace cdo
{
// convert a CDI filetype to cstring
const char *filetype_to_cstr(int filetype);

// convert a CDI datatype to cstring
const char *datatype_to_cstr(int datatype);

int str_to_datatype(std::string const &datatypeStr);

const char *get_steptype_name(int tsteptype);

}  // namespace cdo

// create Taxis(cdi) with added check/setting of taxisType
int cdo_taxis_create(int taxisType);
void cdo_taxis_copy_timestep(int taxisIDdes, int taxisIDsrc);

// cdi
void grid_gen_xvals(int xsize, double xfirst, double xlast, double xinc, double *xvals);
void grid_gen_yvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals);

void cdo_def_table_id(int tableID);

size_t cdo_vlist_gridsizemax(int vlistID);

namespace cdo
{

int inq_att_int(int cdiID, int varID, std::string const &attname);
std::string inq_att_string(int cdiID, int varID, std::string const &attname);

int inq_key_int(int cdiID, int varID, int key);
std::string inq_key_string(int cdiID, int varID, int key);

std::string inq_var_name(int vlistID, int varID);
std::string inq_var_longname(int vlistID, int varID);
std::string inq_var_units(int vlistID, int varID);

HpParams get_healpix_params(int gridID);

std::string get_chunkspec_string(int vlistID, int varID);
// void delete_chunks_dimT(int vlistID);

}  // namespace cdo

#endif
