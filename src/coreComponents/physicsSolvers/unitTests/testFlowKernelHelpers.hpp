/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#ifndef GEOSX_TESTFLOWKERNELHELPERS_HPP
#define GEOSX_TESTFLOWKERNELHELPERS_HPP

#include "common/DataTypes.hpp"

//using namespace geosx;

template<typename T, int NDIM>
using Array = LvArray::Array<T, NDIM, localIndex>;

template<typename T, int NDIM>
using ArrayView = LvArray::ArrayView<T, NDIM, localIndex>;

namespace detail
{

template<typename T, int NDIM>
void setArrayElement( ArrayView<T, NDIM> const & arr,
                      localIndex const dstIndex,
                      localIndex const srcIndex,
                      T const * const data )
{
  localIndex const elemSize = arr.size() / arr.size(0); // assume singleParameterResizeIndex == 0
  T * const elemPtr = arr.data(dstIndex);
  T const * const dataPtr = data + srcIndex * elemSize;

  for (localIndex i = 0; i < elemSize; ++i)
  {
    elemPtr[i] = dataPtr[i];
  }
}

}

/**
 * @brief Specializations of this struct provide a way to create data structures
 *        that flow kernels want (element/material accessors)
 * @tparam FULL whether to create "full" accessors (i.e. nested arrays of Array)
 *              or for a single region (i.e. just an Array)
 */
template<bool FULL>
struct AccessorHelper { };

template<>
struct AccessorHelper<false>
{
  template<int NDIM, typename T>
  using ElementAccessor = Array<T, NDIM>;

  template<int NDIM, typename T>
  using MaterialAccessor = Array<T, NDIM>;

  template<int NDIM, typename T, typename... DIMS >
  static ElementAccessor<NDIM, T>
  makeElementAccessor( T const * const data,
                       localIndex const stencilSize,
                       arraySlice1d<localIndex const> const & stencilRegIndices,
                       arraySlice1d<localIndex const> const & stencilSubRegIndices,
                       arraySlice1d<localIndex const> const & stencilElemIndices,
                       DIMS... otherDims )
  {
    localIndex numElems = 0;
    for (int i = 0; i < stencilSize; ++i)
    {
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    ElementAccessor<NDIM, T> acc( numElems, otherDims... );

    for (int i = 0; i < stencilSize; ++i)
    {
      detail::setArrayElement( acc, stencilElemIndices[i], i, data );
    }

    return acc;
  }

  template<int NDIM, typename T, typename... DIMS >
  static MaterialAccessor<NDIM, T>
  makeMaterialAccessor( T const * const data,
                        localIndex const stencilSize,
                        arraySlice1d<localIndex const> const & stencilRegIndices,
                        arraySlice1d<localIndex const> const & stencilSubRegIndices,
                        arraySlice1d<localIndex const> const & stencilElemIndices,
                        localIndex matIndex, DIMS... otherDims )
  {
    localIndex numElems = 0;
    for (int i = 0; i < stencilSize; ++i)
    {
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    MaterialAccessor<NDIM, T> acc( numElems, 1, otherDims... );

    for (int i = 0; i < stencilSize; ++i)
    {
      detail::setArrayElement( acc, stencilElemIndices[i], i, data );
    }

    return acc;
  }
};

template<>
struct AccessorHelper<true>
{
  template<int NDIM, typename T>
  using ElementAccessor = ElementRegionManager::ElementViewAccessor<Array<T, NDIM>>;

  template<int NDIM, typename T>
  using MaterialAccessor = ElementRegionManager::MaterialViewAccessor<Array<T, NDIM>>;

  template<int NDIM, typename T, typename... DIMS >
  static ElementAccessor<NDIM, T>
  makeElementAccessor( T const * const data,
                       localIndex const stencilSize,
                       arraySlice1d<localIndex const> const & stencilRegIndices,
                       arraySlice1d<localIndex const> const & stencilSubRegIndices,
                       arraySlice1d<localIndex const> const & stencilElemIndices,
                       DIMS... otherDims )
  {
    localIndex numRegions = 0, numSubRegions = 0, numElems = 0;
    for (int i = 0; i < stencilSize; ++i)
    {
      numRegions = std::max( numRegions, stencilRegIndices[i] + 1 );
      numSubRegions = std::max( numSubRegions, stencilSubRegIndices[i] + 1 );
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    ElementAccessor<NDIM, T> acc;
    acc.resize( numRegions );
    for (localIndex kr = 0; kr < numRegions; ++kr)
    {
      acc[kr].resize( numSubRegions );
      for (localIndex ksr = 0; ksr < numSubRegions; ++ksr)
      {
        acc[kr][ksr].resize( numElems, otherDims... );
      }
    }

    for (int i = 0; i < stencilSize; ++i)
    {
      detail::setArrayElement( acc[stencilRegIndices[i]][stencilSubRegIndices[i]], stencilElemIndices[i], i, data );
    }

    return acc;
  }

  template<int NDIM, typename T, typename... DIMS>
  static MaterialAccessor<NDIM, T>
  makeMaterialAccessor( T const * const data,
                        localIndex const stencilSize,
                        arraySlice1d<localIndex const> const & stencilRegIndices,
                        arraySlice1d<localIndex const> const & stencilSubRegIndices,
                        arraySlice1d<localIndex const> const & stencilElemIndices,
                        localIndex matIndex,
                        DIMS... otherDims )
  {
    localIndex numRegions = 0, numSubRegions = 0, numElems = 0;
    for (int i = 0; i < stencilSize; ++i)
    {
      numRegions = std::max( numRegions, stencilRegIndices[i] + 1 );
      numSubRegions = std::max( numSubRegions, stencilSubRegIndices[i] + 1 );
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    MaterialAccessor<NDIM, T> acc;
    acc.resize( numRegions );
    for (localIndex kr = 0; kr < numRegions; ++kr)
    {
      acc[kr].resize( numSubRegions );
      for (localIndex ksr = 0; ksr < numSubRegions; ++ksr)
      {
        acc[kr][ksr].resize( matIndex + 1 );
        acc[kr][ksr][matIndex].resize( numElems, 1, otherDims... );
      }
    }

    for (int i = 0; i < stencilSize; ++i)
    {
      detail::setArrayElement( acc[stencilRegIndices[i]][stencilSubRegIndices[i]][matIndex],
                               stencilElemIndices[i],
                               i,
                               data );
    }

    return acc;
  }
};

#endif //GEOSX_TESTFLOWKERNELHELPERS_HPP
