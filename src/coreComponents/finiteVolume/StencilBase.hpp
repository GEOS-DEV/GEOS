/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CellElementStencilTPFA.hpp
 */

#ifndef GEOSX_FINITEVOLUME_STENCILBASE_HPP_
#define GEOSX_FINITEVOLUME_STENCILBASE_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{

template< typename LEAFCLASSTRAITS >
class StencilWrapperBase
{
public:
  StencilWrapperBase( typename LEAFCLASSTRAITS::IndexContainerType & elementRegionIndices,
                      typename LEAFCLASSTRAITS::IndexContainerType & elementSubRegionIndices,
                      typename LEAFCLASSTRAITS::IndexContainerType & elementIndices,
                      typename LEAFCLASSTRAITS::WeightContainerType & weights ):
    m_elementRegionIndices( elementRegionIndices.toViewConst() ),
    m_elementSubRegionIndices( elementSubRegionIndices.toViewConst() ),
    m_elementIndices( elementIndices.toViewConst() ),
    m_weights( weights.toViewConst() )
  {};

  /// Default copy constructor
  StencilWrapperBase( StencilWrapperBase const & ) = default;

  /// Default move constructor
  StencilWrapperBase( StencilWrapperBase && ) = default;

  /// Deleted copy assignment operator
  StencilWrapperBase & operator=( StencilWrapperBase const & ) = delete;

  /// Deleted move assignment operator
  StencilWrapperBase & operator=( StencilWrapperBase && ) = delete;

  /**
   * @brief Const access to the element regions indices.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType getElementRegionIndices() const { return m_elementRegionIndices; }

  /**
   * @brief Const access to the element subregions indices.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType getElementSubRegionIndices() const { return m_elementSubRegionIndices; }

  /**
   * @brief Const access to the element indices.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType getElementIndices() const { return m_elementIndices; }

  /**
   * @brief Const access to the stencil weights.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::WeightContainerViewConstType getWeights() const { return m_weights; }

//  template< typename PERMTYPE >
//  void computeTransmissibility( localIndex iconn,
//                                PERMTYPE permeability,
//                                real64 (& transmissibility)[2] );
//
//  template< typename PERMTYPE >
//  void dTrans_dPressure( localIndex iconn,
//                         PERMTYPE dPerm_dPressure,
//                         real64 (&dTrans_dPressure )[2] );


  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const = 0;

protected:
  typename LEAFCLASSTRAITS::IndexContainerViewConstType m_elementRegionIndices;

  typename LEAFCLASSTRAITS::IndexContainerViewConstType m_elementSubRegionIndices;

  typename LEAFCLASSTRAITS::IndexContainerViewConstType m_elementIndices;

  typename LEAFCLASSTRAITS::WeightContainerViewConstType m_weights;

};


/**
 * @class StencilBase
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
class StencilBase
{
public:

  StencilBase():
    m_elementRegionIndices(),
    m_elementSubRegionIndices(),
    m_elementIndices(),
    m_weights(),
    m_connectorIndices()
  {}

  /**
   * @brief Destructor.
   */
  virtual ~StencilBase() = default;

  /**
   * @brief Constructor.
   */
  StencilBase( StencilBase const & ) = default;

  /**
   * @brief Move constructor.
   */
  StencilBase( StencilBase && ) = default;

  /**
   * @brief Reserve the size of the stencil.
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size );

  /**
   * @brief Move the data arrays associated with the stencil to a specified
   *   memory space.
   * @param space The target memory space.
   *
   * @note The existence of this function indicates we need to redesign the
   * stencil classes.
   */
  virtual void move( LvArray::MemorySpace const space );


  /**
   * @brief Add an entry to the stencil.
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] elementRegionIndices The element region indices for each point in the stencil entry
   * @param[in] elementSubRegionIndices The element sub-region indices for each point in the stencil entry
   * @param[in] elementIndices The element indices for each point in the stencil entry
   * @param[in] weights The weights each point in the stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) = 0;

  /**
   * @brief Zero weights for a stencil entry.
   * @param[in] connectorIndex The index of the connector element that the stencil acts across for which the weights are
   *                           to be zero.
   * @return True if a valid connectorIndex was found, and had its corresponding weights set to zero.
   */
  virtual bool zero( localIndex const connectorIndex );

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const = 0;

  /**
   * @brief Set the name used in data movement logging callbacks.
   * @param name the name prefix for the stencil's data arrays
   */
  void setName( string const & name );

  /**
   * @brief Const access to the element regions indices.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType getElementRegionIndices() const { return m_elementRegionIndices.toViewConst(); }

  /**
   * @brief Const access to the element subregions indices.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType getElementSubRegionIndices() const { return m_elementSubRegionIndices.toViewConst(); }

  /**
   * @brief Const access to the element indices.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType getElementIndices() const { return m_elementIndices.toViewConst(); }

  /**
   * @brief Const access to the stencil weights.
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::WeightContainerViewConstType getWeights() const { return m_weights.toViewConst(); }

protected:
  /// The container for the element region indices for each point in each stencil
  typename LEAFCLASSTRAITS::IndexContainerType m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  typename LEAFCLASSTRAITS::IndexContainerType m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  typename LEAFCLASSTRAITS::IndexContainerType m_elementIndices;

  /// The container for the weights for each point in each stencil
  typename LEAFCLASSTRAITS::WeightContainerType m_weights;

  /// The map that provides the stencil index given the index of the underlying connector object.
  map< localIndex, localIndex > m_connectorIndices;

};



template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 2 );
  m_elementSubRegionIndices.reserve( size * 2 );
  m_elementIndices.reserve( size * 2 );
  m_weights.reserve( size * 2 );
}


template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
bool StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::zero( localIndex const connectorIndex )
{
  return
    executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for( localIndex i = 0; i < static_cast< LEAFCLASS * >(this)->stencilSize( connectorIndex ); ++i )
    {
      m_weights[connectionListIndex][i] = 0;
    }
  } );
}

template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::setName( string const & name )
{
  m_elementRegionIndices.setName( name + "/elementRegionIndices" );
  m_elementSubRegionIndices.setName( name + "/elementSubRegionIndices" );
  m_elementIndices.setName( name + "/elementIndices" );
  m_weights.setName( name + "/weights" );
}

template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::move( LvArray::MemorySpace const space )
{
  m_elementRegionIndices.move( space, true );
  m_elementSubRegionIndices.move( space, true );
  m_elementIndices.move( space, true );
  m_weights.move( space, true );
}


//template< typename LEAFCLASSTRAITS >
//template< typename PERMTYPE >
//void StencilWrapperBase<LEAFCLASSTRAITS>::computeTransmissibility( localIndex iconn,
//                                                                   PERMTYPE permeability,
//                                                                   real64 (& transmissibility)[2] )
//{
//  localIndex const er0  =  m_elementRegionIndices[iconn][0];
//  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
//  localIndex const ei0  =  m_elementIndices[iconn][0];
//
//  localIndex const er1  =  m_elementRegionIndices[iconn][1];
//  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
//  localIndex const ei1  =  m_elementIndices[iconn][1];
//
//  real64 const t0 = m_weights[iconn][0] * permeability[er0][esr0][ei0];
//  real64 const t1 = m_weights[iconn][1] * permeability[er1][esr1][ei1];
//
//  real64 const harmonicWeight   = t0*t1 / (t0+t1);
//  real64 const arithmeticWeight = (t0+t1)/2;
//
//  real64 const meanPermCoeff = 1.0; //TODO make it a member
//
//  transmissibility[0] = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
//  transmissibility[1] = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
//}
//
//template< typename LEAFCLASSTRAITS >
//template< typename PERMTYPE >
//void StencilWrapperBase<LEAFCLASSTRAITS>::dTrans_dPressure( localIndex iconn,
//                                                            PERMTYPE dPerm_dPressure,
//                                                            real64 (&dTrans_dPressure )[2] )
//{
//  localIndex const er0  =  m_elementRegionIndices[iconn][0];
//  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
//  localIndex const ei0  =  m_elementIndices[iconn][0];
//
//  localIndex const er1  =  m_elementRegionIndices[iconn][1];
//  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
//  localIndex const ei1  =  m_elementIndices[iconn][1];
//
//  real64 const dt0 = m_weights[iconn][0] * dPerm_dPressure[er0][esr0][ei0];
//  real64 const dt1 = m_weights[iconn][1] * dPerm_dPressure[er1][esr1][ei1];
//
//  // TODO fix this with proper derivative calculation.
//  dTrans_dPressure[0] = dt0;
//  dTrans_dPressure[1] = dt1;
//}



} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_STENCILBASE_HPP_ */
