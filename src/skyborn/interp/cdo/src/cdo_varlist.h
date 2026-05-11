#ifndef CDO_VARLIST_H
#define CDO_VARLIST_H

#include <cassert>
#include <string>
#include <vector>
#include <map>

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"

enum struct MapFlag
{
  Undefined = 0,
  Left = 1,
  Right = 2,
  Intersect = 3
};

namespace CmpVarList
{
constexpr int Name = 1;
constexpr int Grid = 2;
constexpr int NumLevels = 4;
constexpr int GridSize = 8;
constexpr int Dim = GridSize | NumLevels | Grid;
constexpr int All = Name | Dim;
};  // namespace CmpVarList

struct CdoVar
{
  size_t counter{ 0 };
  std::string name{};
  std::string longname{};
  std::string units{};
  std::string stdname{};
  MemType memType{ MemType::Native };
  int gridID{ -1 };
  int zaxisID{ -1 };
  int gridType{ -1 };
  int zaxisType{ -1 };
  int timeType{ -1 };
  int stepType{ -1 };
  size_t gridsize{ 0 };
  int nlevels{ 0 };
  int dataType{ -1 };
  double missval{ 0 };
  double addOffset{ 0.0 };
  double scaleFactor{ 1.0 };
  int code{ 0 };
  int param{ 0 };
  int nwpv{ 1 };  // number of words per value; real:1  complex:2
  bool isConstant{ false };
  bool isPacked{ false };
  int ID{ 0 };
};

using CdoVars = std::vector<CdoVar>;

void cdoVars_init(CdoVars &cdoVars, int vlistID);

class VarList
{
public:
  CdoVars vars{};
  int vlistID{ CDI_UNDEFID };

  VarList() {}
  explicit VarList(int _vlistID);

  // clang-format off
  void isInit() const { assert(vlistID != CDI_UNDEFID); }
  int numVars() const noexcept { isInit(); return static_cast<int>(vars.size()); }
  int maxFields() const noexcept { isInit(); return m_maxFields; }
  int numSteps() const noexcept { isInit(); return m_numSteps; }
  int numZaxes() const noexcept { isInit(); return m_numZaxes; }
  int numGrids() const noexcept { isInit(); return m_numGrids; }
  int numConstVars() const noexcept { isInit(); return m_numConstVars; }
  int numVaryingVars() const noexcept { isInit(); return m_numVaryingVars; }
  size_t gridsizeMax() const noexcept { isInit(); return m_gridsizeMax; }
  // clang-format on

private:
  // clang-format off
  int m_maxFields{ 0 };
  int m_numSteps{ 0 };
  int m_numZaxes{ 0 };
  int m_numGrids{ 0 };
  int m_numConstVars{ 0 };
  int m_numVaryingVars{ 0 };
  size_t m_gridsizeMax{ 0 };
  void set_num_const_vars(CdoVars const &cdoVars);
  void set_num_varying_vars(CdoVars const &cdoVars);
  // clang-format on
};

struct VarIDs
{
  int sgeopotID{ CDI_UNDEFID };
  int geopotID{ CDI_UNDEFID };
  int taID{ CDI_UNDEFID };
  int psID{ CDI_UNDEFID };
  int lnpsID{ CDI_UNDEFID };
  int lnpsID2{ CDI_UNDEFID };
  int gheightID{ CDI_UNDEFID };
  int husID{ CDI_UNDEFID };
  int clwcID{ CDI_UNDEFID };
  int ciwcID{ CDI_UNDEFID };
};

VarIDs varList_search_varIDs(VarList const &varList, int numFullLevels);

void vlist_compare(int vlistID1, int vlistID2, int cmpFlag);

void varList_compare(VarList const &varList1, VarList const &varList2, int cmpFlag = CmpVarList::All);
void varList_map(VarList const &varList1, VarList const &varList2, MapFlag mapFlag, std::map<int, int> &mapOfVarIDs);

void varList_set_memtype(VarList &varList, MemType memType);
void varList_set_unique_memtype(VarList &varList);
int varList_get_psvarid(VarList const &varList, int zaxisID);

void print_found_variables(VarIDs const &varIDs, CdoVars const &vars);

#endif
