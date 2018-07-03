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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file typedefs.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include "intrinsic_typedefs.h"
#include "legacy/ArrayT/array.h"
#include "legacy/ArrayT/Array2dT.h"
#include "math/TensorT/TensorT.h"
#include "common/SortedArray.hpp"

/// set number of degrees of freedom
const unsigned int nsdof = 3;

/// define the templated R2TensorT<> with nsdof as its template argument
typedef R2TensorT<nsdof>    R2Tensor;

/// define the templated R2SymTensorT<> with nsdof as its template argument
typedef R2SymTensorT<nsdof> R2SymTensor;

/// define the templated R1TensorT<> with nsdof as its template argument
typedef R1TensorT<nsdof>    R1Tensor;


typedef array<int> array<integer>;
typedef array<realT> array<real64>;
typedef array<std::string> array<string>;
typedef array<array<integer>*> iArrayPtrs;
typedef array<array<real64>*> rArrayPtrs;

typedef Array2dT<int> iArray2d;
typedef Array2dT<localIndex> lArray2d;
typedef Array2dT<globalIndex> gArray2d;
typedef Array2dT<realT> rArray2d;
typedef Array2dT<std::pair<int,localIndex> > pArray2d;

typedef array<localIndex> lArray1d;
typedef array<globalIndex> gArray1d;
typedef array<std::pair<int,localIndex> > pArray1d;


template< typename T >
using set = SortedArray<T>;

typedef SortedArray<int> iSet;
typedef SortedArray<localIndex> lSet;
typedef SortedArray<globalIndex> gSet;
typedef SortedArray<std::pair<int,localIndex> > pSet;
typedef SortedArray<std::string> sSet;

template<typename T> struct type_name
{
  static const char* name();  // { static_assert(false, "You are missing a
                              // DECL_TYPE_NAME"); }
};

template<> struct type_name<int>                { static const char* name() {return "int";}
};
template<> struct type_name<unsigned int>       { static const char* name() {return "unsigned int";}
};
template<> struct type_name<long long>          { static const char* name() {return "long long";}
};
template<> struct type_name<unsigned long long> { static const char* name() {return "unsigned long long";}
};
template<> struct type_name<float>              { static const char* name() {return "float";}
};
template<> struct type_name<double>             { static const char* name() {return "double";}
};

/// stores pointer to elemenRegion and localIndex of element
//class ElementRegionT;
//typedef std::pair< ElementRegionT*, localIndex > ElementIdPair;

template <class Type>
struct locallyIndexedType
{
  localIndex m_index;
  Type m_x;

public:
  locallyIndexedType(localIndex indx, Type x): m_index(indx),m_x(x){ /** empty
                                                                     **/}
  virtual localIndex GetIndex(){return m_index;}
  virtual Type GetValue(){return m_x;}
};

typedef locallyIndexedType<realT> locallyIndexedReal;


#endif /* TYPEDEFS_H_ */
