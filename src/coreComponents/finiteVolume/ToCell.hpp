//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Fri 04 2024.
//////////////////////////////////////////////////////////////////////////////////

#ifndef GEOSX_TOCELL_HPP
#define GEOSX_TOCELL_HPP

#include "StencilBase.hpp"


namespace geos
{

class ToCellWrapperBase : public StencilWrapperBase< TwoPointStencilTraits >
{

public:

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  ToCellWrapperBase( const IndexContainerType & elementRegionIndices, const IndexContainerType & elementSubRegionIndices,
                     const IndexContainerType & elementIndices, const WeightContainerType & weights ):
    StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ){};

  localIndex size() const
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const override
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  using StencilWrapperBase< TwoPointStencilTraits >::computeWeights;

  /**
   * @brief Compute the stabilization weights
   * @param[in] iconn connection index
   * @param[out] stabilizationWeight view weights
   */
  GEOS_HOST_DEVICE
  void computeStabilizationWeights( localIndex iconn,
                                    real64 ( & stabilizationWeight )[1][2] ) const
  {
    GEOS_UNUSED_VAR( iconn, stabilizationWeight );
  }

  /**
   * @brief Remove the contribution of the aperture from the weight in the stencil (done before aperture update)
   *
   * @param iconn connection index
   * @param hydraulicAperture hydraulic apertures of the fractures
   */
  GEOS_HOST_DEVICE
  void removeHydraulicApertureContribution( localIndex const iconn,
                                            ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const;

  /**
   * @brief Add the contribution of the aperture to the weight in the stencil (done after aperture update)
   *
   * @param iconn connection index
   * @param hydraulicAperture hydraulic apertures of the fractures
   */
  GEOS_HOST_DEVICE
  void addHydraulicApertureContribution( localIndex const iconn,
                                         ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const;



};

template< typename LEAFCLASS >
class ToCellBase : public StencilBase< TwoPointStencilTraits, LEAFCLASS >
{
public:
//  /**
//   * @brief Move the data arrays associated with the stencil to a specified
//   *   memory space.
//   * @param space The target memory space.
//   *
//   * @note The existence of this function indicates we need to redesign the
//   * stencil classes.
//   */
//  virtual void move( LvArray::MemorySpace const space ) override;
  using StencilBase< TwoPointStencilTraits, LEAFCLASS >::move;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override { return StencilBase< TwoPointStencilTraits, LEAFCLASS >::m_elementRegionIndices.size( 0 ); }


  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const
  {
    GEOS_UNUSED_VAR( index );
    return TwoPointStencilTraits::maxStencilSize;
  }

};
GEOS_HOST_DEVICE
inline void
ToCellWrapperBase::
  removeHydraulicApertureContribution( localIndex const iconn,
                                       ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  // only the fracture side is modified, k=1
  localIndex constexpr k = 1;
  localIndex const er = m_elementRegionIndices[iconn][k];
  localIndex const esr = m_elementSubRegionIndices[iconn][k];
  localIndex const ei = m_elementIndices[iconn][k];

  m_weights[iconn][k] = m_weights[iconn][k] * hydraulicAperture[er][esr][ei];
}

GEOS_HOST_DEVICE
inline void
ToCellWrapperBase::
  addHydraulicApertureContribution( localIndex const iconn,
                                    ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  // only the fracture side is modified, k=1
  localIndex constexpr k = 1;
  localIndex const er = m_elementRegionIndices[iconn][k];
  localIndex const esr = m_elementSubRegionIndices[iconn][k];
  localIndex const ei = m_elementIndices[iconn][k];

  m_weights[iconn][k] = m_weights[iconn][k] / hydraulicAperture[er][esr][ei];
}


GEOS_HOST_DEVICE
template< typename LEAF >
void ToCellBase< LEAF >::add( localIndex const numPts,
                              localIndex const * const elementRegionIndices,
                              localIndex const * const elementSubRegionIndices,
                              localIndex const * const elementIndices,
                              real64 const * const weights,
                              localIndex const connectorIndex )
{
  GEOS_ERROR_IF_NE_MSG( numPts, 2, "Number of cells in faceToCell stencil should be 2" );

  localIndex const oldSize = StencilBase< TwoPointStencilTraits, LEAF >::m_elementRegionIndices.size( 0 );
  localIndex const newSize = oldSize + 1;
  StencilBase< TwoPointStencilTraits, LEAF >::m_elementRegionIndices.resize( newSize, numPts );
  StencilBase< TwoPointStencilTraits, LEAF >::m_elementSubRegionIndices.resize( newSize, numPts );
  StencilBase< TwoPointStencilTraits, LEAF >::m_elementIndices.resize( newSize, numPts );
  StencilBase< TwoPointStencilTraits, LEAF >::m_weights.resize( newSize, numPts );

  for( localIndex a=0; a<numPts; ++a )
  {
    StencilBase< TwoPointStencilTraits, LEAF >::m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    StencilBase< TwoPointStencilTraits, LEAF >::m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    StencilBase< TwoPointStencilTraits, LEAF >::m_elementIndices( oldSize, a ) = elementIndices[a];
    StencilBase< TwoPointStencilTraits, LEAF >::m_weights( oldSize, a ) = weights[a];
  }
  StencilBase< TwoPointStencilTraits, LEAF >::m_connectorIndices[connectorIndex] = oldSize;
}
} // geos

#endif //GEOSX_TOCELL_HPP
