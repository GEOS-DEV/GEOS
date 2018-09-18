/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


#ifndef __TABLE__
#define __TABLE__

#ifdef GEOSX_USE_ATK
#include "slic/slic.hpp"
#endif

#include "codingUtilities/StringUtilities.hpp"
#include "codingUtilities/Utilities.hpp"
#include <vector>

// If we wanted to be super fancy we could move this inside the class and then specialize on the dimension of the table class.
// Seems like too much of a headache though - a separate interpolation class might be a better solution.
namespace TableInterpolation{
  enum Order{
	   	zeroth,
	 	linear
  };
}

template<unsigned int T_dim, class T>
class Table
{
public:

private:
 bool m_set;
 bool m_zeroGradient;
 TableInterpolation::Order m_interpolation;
 typename std::vector<T>::size_type m_size;

 std::vector<T> m_p;
 std::vector<unsigned int> m_mult;
 std::vector<std::vector<realT> > m_x;

public:
  ///Default constructor
  Table() : m_set(false), m_zeroGradient(false),m_interpolation(TableInterpolation::linear), m_size(0), m_p(), m_mult(T_dim), m_x(T_dim)
  {
  }

  template<class ARRAY>
  void SetGrid(const ARRAY& x)
  {
    m_size = 1;
    for (unsigned int i = 0; i < T_dim; i++)
    {
      m_x[i] = x[i];
      //validate user input
      {
        realT xlast = -std::numeric_limits<realT>::max();
        for( int j = 0; j < x[i].size(); j++)
        {
          if(x[i][j] <= xlast) {
#ifdef GEOSX_USE_ATK
            SLIC_ERROR("Table:SetGrid - ticks for axis " + std::to_string(i) +" must be a monotonic increasing vector of values");
#endif
          }
          xlast = x[i][j];
        }
      }

      m_mult[i] = 1;
      m_size *= x[i].size();
      m_mult[i] = i == 0 ? 1 : m_mult[i-1] * x[i-1].size();
    }

    m_p.resize(m_size);
  }

  template<class VECTOR>
  void SetAxisValues(const unsigned dim, const VECTOR& x){
      m_x[dim] = x;

      m_size = 1;
      for (unsigned int i = 0; i < T_dim; i++)
      {
          m_size *= m_x[i].size();
          m_mult[i] = i == 0 ? 1 : m_mult[i-1] * m_x[i-1].size();
      }
      m_p.resize(m_size);

  }

  template<class ARRAY>
  void SetValues(const ARRAY& p)
  {
    m_p = p;
    if(m_size == 0 || m_p.size() != m_size)
#ifdef GEOSX_USE_ATK
      SLIC_ERROR("Table:SetValues - Must set axes before attempting to set data");
#endif
    m_set = true;
  }

  bool SetZeroGradient(const bool zeroGrad = true)
  {
    m_zeroGradient = zeroGrad;
    return m_zeroGradient;
  }

  void SetInterpolation(TableInterpolation::Order interpolate ){m_interpolation =  interpolate;}

  bool ZeroGradient() const
  {
    return m_zeroGradient;
  }

  unsigned int Dimensions() const { return T_dim; }

  unsigned int Dimension(const unsigned dim) const
  {
    if(dim >= T_dim) {
#ifdef GEOSX_USE_ATK
      SLIC_ERROR("Table dimension out of range");
#endif
    }
    return this->m_x[dim].size();
  }

  const std::vector<realT>& AxisValues(const unsigned int dim) const
  {
    if(dim >= T_dim) {
#ifdef GEOSX_USE_ATK
      SLIC_ERROR("Table dimension out of range");
#endif
    }
    return m_x[dim];
  }

  unsigned int ValuesPerBlock(const unsigned int dim) const
  {
    if(dim >= T_dim) {
#ifdef GEOSX_USE_ATK
      SLIC_ERROR("Table dimension out of range");
#endif
    }
    return m_mult[dim];
  }

  const std::vector<realT>& Values() const
  {
    return m_p;
  }

//  template<class ARRAY>
//  static unsigned int LookupIndex(const ARRAY& indices,
//                                  const ARRAY& dims)
//  {
//    unsigned int ii = 0;
//    for(unsigned int i = 0; i < dims.size(); i++)
//    {
//      unsigned int icurr = indices[i];
//      for(unsigned int j = 0; j < i; j++)
//        icurr *= dims[j];
//      ii += icurr;
//    }
//    return ii;
//  }

  realT LookupIndices(const unsigned int dim, const realT x, unsigned int& i0, unsigned int& i1) const
  {
    i0 = 0;
    i1 = 0;
    realT weight0 = 1.0;

    const std::vector<realT>& xcurr = m_x[dim];
    const unsigned int ilast = xcurr.size() - 1;

    if (x < *xcurr.begin() || ilast == 0)
    {
      i0 = 0;
      i1 = 0;
    }
    else if (x >= xcurr.back())
    {
      i0 = ilast;
      i1 = i0;
    }
    else
    {
      // We are within the table
      for (i1 = 1; i1 <= ilast; i1++)
      {
        if (xcurr[i1] >= x)
          break;
      }
      i0 = i1 - 1;
      weight0 = 1.0 - (x - xcurr[i0]) / (xcurr[i1] - xcurr[i0]);
    }

    return weight0;
  }

private:

  template<class ARRAY>
  T Lookup(const ARRAY& xx, const unsigned int dim, unsigned int& i, const TableInterpolation::Order interpolate = TableInterpolation::linear) const
  {
    //get the weighting and dimensional indices (i0 & i1)
    unsigned int i0, i1;
    const realT weightCurrent = LookupIndices(dim, xx[dim], i0, i1);

    //get the components
    unsigned int icurr = i + m_mult[dim] * i0;
    T ret = dim > 0 ? Lookup(xx, dim-1, icurr, interpolate) : m_p[icurr];

    if(interpolate == TableInterpolation::linear)
    {
      ret *= weightCurrent;
      icurr = i + m_mult[dim] * i1;
      T ret1 = dim > 0 ? Lookup(xx, dim-1, icurr, interpolate) : m_p[icurr];
      ret1 *= (1.0 - weightCurrent);
      ret += ret1;
    }

    //return the result
    return ret;
  }

  template<class ARRAY>
  T Gradient(const ARRAY& xx, const unsigned int dim, const unsigned int gdim, unsigned int& i) const
  {
    //get the weighting and dimensional indices (i0 & i1)
    unsigned int i0, i1;
    const realT weightCurrent = LookupIndices(dim, xx[dim], i0, i1);

    //get the components
    unsigned int icurr = i + m_mult[dim] * i0;
    T ret0 = dim > 0 ? Gradient(xx, dim-1, gdim, icurr) : m_p[icurr];
    icurr = i + m_mult[dim] * i1;
    T ret = dim > 0 ? Gradient(xx, dim-1, gdim, icurr) : m_p[icurr];
    if(dim != gdim)
    {
      ret0 *= weightCurrent;
      ret *= (1.0 - weightCurrent);
      ret += ret0;
    }
    else
    {
      ret -= ret0;
      realT fct = m_x[dim][i1] - m_x[dim][i0];
      if( !isZero(fct) )
        fct = 1.0 / fct;
      ret *= fct;
    }

    //return the result
    return ret;
  }

 public:

  ///Look up a value in the list
  /**
   \param xx N-dimensional point to lookup
   \param interpolate A flag to determine how to interpolate the point on the grid
   \return Value corresponding to N-dimensional point
   */
  template<class ARRAY>
  T Lookup(const ARRAY& xx, const TableInterpolation::Order interpolation) const
  {
    T ret = static_cast<T>(0);
    if(T_dim < 1)
      return ret;

    unsigned int i = 0, dim = T_dim - 1;
    ret = Lookup(xx, dim, i, interpolation);
    return ret;
  }

  //In most cases this function should be used as a table should know how it wants to be interpolated
  template<class ARRAY>
  T Lookup(const ARRAY& xx) const{
	 return Lookup(xx, m_interpolation);
  }

  ///Look up a gradient in the list
  /**
   If the gradient requested lies directly on a C1 discontinuity
   then the left segment slope is used; if the gradient requested
   lies beyond the boundaries, a 0 gradient is returned
   \param xx Point to lookup
   \param interpolate A flag to determine whether to linearly interpolate the point on the grid
   \return Value corresponding to txyz point
   */
  template<class ARRAY>
  T Gradient(const ARRAY& xx, const unsigned int gradientDimension = 0) const
  {
    T ret = 0;
    if(m_zeroGradient || !m_set || T_dim < 1)
      return ret;
    if(gradientDimension >= T_dim) {
#ifdef GEOSX_USE_ATK
      SLIC_ERROR("Table:Gradient - cannot have gradient dimension >= dimension");
#endif
    }

    unsigned int i = 0, dim = T_dim - 1;
    return Gradient(xx, dim, gradientDimension, i);
  }
};
#endif //__TABLE__
