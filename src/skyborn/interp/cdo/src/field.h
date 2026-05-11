/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_H
#define FIELD_H

#include <cstdio>
#include <utility>
#include "varray.h"
#include "cdo_options.h"
#include "cdo_varlist.h"
#include "cdo_vlist.h"

// clang-format off
const auto memtype_is_float_float   = [](auto a, auto b) noexcept { return (a == MemType::Float  && b == MemType::Float); };
const auto memtype_is_double_float  = [](auto a, auto b) noexcept { return (a == MemType::Double && b == MemType::Float); };
const auto memtype_is_float_double  = [](auto a, auto b) noexcept { return (a == MemType::Float  && b == MemType::Double); };
const auto memtype_is_double_double = [](auto a, auto b) noexcept { return (a == MemType::Double && b == MemType::Double); };
// clang-format on

// clang-format off
template <typename FUNC, typename FIELD>
inline auto
field_operation(FUNC func, FIELD &field)
-> decltype(func(field.vec_f))
{
  if      (field.memType == MemType::Float)  return func(field.vec_f);
  else if (field.memType == MemType::Double) return func(field.vec_d);
  else throw std::runtime_error("Type of field unsupported!");
}

template <typename FUNC, typename FIELD, typename... ARGS>
inline auto
field_operation(FUNC func, FIELD &field, ARGS &...args)
-> decltype(func(field.vec_f, args...))
{
  if      (field.memType == MemType::Float)  return func(field.vec_f, args...);
  else if (field.memType == MemType::Double) return func(field.vec_d, args...);
  else throw std::runtime_error("Type of field unsupported!");
}

template<typename FUNC, typename FIELD1, typename FIELD2>
inline auto
field_operation2(FUNC func, FIELD1& field1, FIELD2& field2)
-> decltype(func(field1.vec_f, field2.vec_d))
{
  if      (memtype_is_float_float(field1.memType, field2.memType))   return func(field1.vec_f, field2.vec_f);
  else if (memtype_is_float_double(field1.memType, field2.memType))  return func(field1.vec_f, field2.vec_d);
  else if (memtype_is_double_float(field1.memType, field2.memType))  return func(field1.vec_d, field2.vec_f);
  else if (memtype_is_double_double(field1.memType, field2.memType)) return func(field1.vec_d, field2.vec_d);
  else throw std::runtime_error("Type of fields unsupported!");
}

template<typename FUNC, typename FIELD1, typename FIELD2, typename... ARGS>
inline auto
field_operation2(FUNC func, FIELD1& field1, FIELD2& field2, ARGS&... args)
-> decltype(func(field1.vec_f, field2.vec_d, args...))
{
  if      (memtype_is_float_float(field1.memType, field2.memType))   return func(field1.vec_f, field2.vec_f, args...);
  else if (memtype_is_float_double(field1.memType, field2.memType))  return func(field1.vec_f, field2.vec_d, args...);
  else if (memtype_is_double_float(field1.memType, field2.memType))  return func(field1.vec_d, field2.vec_f, args...);
  else if (memtype_is_double_double(field1.memType, field2.memType)) return func(field1.vec_d, field2.vec_d, args...);
  else throw std::runtime_error("Type of fields unsupported!");
}
// clang-format on

enum field_flag
{
  FIELD_VEC = 2,   // allocated memory
  FIELD_FLT = 4,   // 32-bit float
  FIELD_DBL = 8,   // 64-bit float
  FIELD_NAT = 16,  // native: 32-bit float for 32-bit float data, otherweise 64-bit float
};

// clang-format off
class  // Field
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
Field
// clang-format on
{
public:
  int fpeRaised = 0;
  int nwpv = 1;  // number of words per value; real:1  complex:2
  int grid = -1;
  MemType memType = MemType::Native;  // MemType::Float or MemType::Double

  size_t gridsize = 0;
  size_t size = 0;
  size_t nsamp = 0;

  size_t numMissVals = 0;
  double missval = 0;

  Varray<float> vec_f;
  Varray<double> vec_d;
  Varray<double> weightv;

  Field() {}
  void init(const CdoVar &var);
  void resize(size_t count);
  void resize(size_t count, double value);
  void resizef(size_t count);
  void resizef(size_t count, float value);
  bool empty() const;
  void check_gridsize() const;

  double
  operator[](size_t pos) const
  {
    if (pos < size)
    {
      if (memType == MemType::Float)
        return vec_f[pos];
      else if (memType == MemType::Double)
        return vec_d[pos];
      else
        throw std::runtime_error("Type of field unsupported!");
    }

    return double();
  }

  bool
  hasData() const
  {
    return (memType == MemType::Float) ? !vec_f.empty() : !vec_d.empty();
  }

private:
  size_t m_count = 0;
};

// clang-format off
class  // Field3D
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
Field3D : public Field
// clang-format on
{
public:
  size_t nlevels = 0;

  Field3D() {}
  void init(const CdoVar &var);
};

struct FieldInfo
{
  int varID = 0;
  int levelID = 0;
  void
  set(int _varID, int _levelID)
  {
    this->varID = _varID;
    this->levelID = _levelID;
  }
  std::pair<int, int>
  get() const
  {
    return std::make_pair(varID, levelID);
  }
};

using FieldVector = std::vector<Field>;
using FieldVector2D = std::vector<FieldVector>;
using FieldVector3D = std::vector<FieldVector2D>;

using Field3DVector = std::vector<Field3D>;

void field_fill(Field &field, double value);
void field_ncopy(size_t n, Field const &fieldIn, Field &fieldOut);
void field_copy(Field const &fieldIn, Field &fieldOut);
void field_copy(const Field3D &fieldIn, Field3D &fieldOut);
void field_copy(const Field3D &fieldIn, int levelID, Field &fieldOut);
void field_add(Field &field1, const Field3D &field2, int levelID);
size_t field_num_NANs(Field const &field);
size_t field_num_mv(Field &field);
MinMax field_min_max(Field const &field);

template <class UnaryOperation>
void
field_transform(Field const &fieldIn, Field &fieldOut, UnaryOperation unary_op)
{
  if (fieldIn.memType == MemType::Float && fieldOut.memType == MemType::Float)
    varray_transform(fieldIn.vec_f, fieldOut.vec_f, unary_op);
  else
    varray_transform(fieldIn.vec_d, fieldOut.vec_d, unary_op);
}

#endif /* FIELD_H */
