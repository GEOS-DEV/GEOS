/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPointFluxApproximation.hpp
 */

#ifndef GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
#define GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_

#include "finiteVolume/FluxApproximationBase.hpp"

namespace geosx
{

/**
 * @class TwoPointFluxApproximation
 *
 * Provides management of the interior stencil points when using a two-point flux approximation.
 */
class TwoPointFluxApproximation : public FluxApproximationBase
{
public:

  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static string catalogName() { return "TwoPointFluxApproximation"; }

  TwoPointFluxApproximation() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the TwoPointFluxApproximation in the data repository
   * @param parent the parent group of this group.
   */
  TwoPointFluxApproximation( string const & name, dataRepository::Group * const parent );

  /**
   * @brief View keys.
   */
  struct viewKeyStruct : FluxApproximationBase::viewKeyStruct
  {
    /// @return The key for edfStencil
    static constexpr char const * edfmStencilString() { return "edfmStencil"; }
    /// @return The key for faceElementToCellStencil
    static constexpr char const * faceToCellStencilString() { return "faceElementToCellStencil"; }
    /// @return The key for the meanPermCoefficient
    static constexpr char const * meanPermCoefficientString() { return "meanPermCoefficient"; }
    /// @return The key for the usePEDFM flag
    static constexpr char const * usePEDFMString() { return "usePEDFM"; }
  };

protected:
  virtual void registerCellStencil( Group & stencilGroup ) const override;

  virtual void computeCellStencil( MeshLevel & mesh ) const override;

  virtual void registerFractureStencil( Group & stencilGroup ) const override;

  virtual void addToFractureStencil( MeshLevel & mesh,
                                     string const & faceElementRegionName,
                                     bool const initFlag ) const override;

  virtual void registerBoundaryStencil( Group & stencilGroup,
                                        string const & setName ) const override;

  virtual void computeBoundaryStencil( MeshLevel & mesh,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & faceSet ) const override;

  virtual void registerAquiferStencil( Group & stencilGroup,
                                       string const & setName ) const override;

  virtual void computeAquiferStencil( DomainPartition & domain,
                                      MeshLevel & mesh ) const override;

  virtual void addEmbeddedFracturesToStencils( MeshLevel & mesh,
                                               string const & embeddedSurfaceRegionName ) const override;

  /**
   * @brief adds fracture-fracture connections to the surfaceElement stencil
   * @param mesh the mesh object
   * @param embeddedSurfaceRegionName name of the fracture region
   */
  void addFractureFractureConnections( MeshLevel & mesh,
                                       string const & embeddedSurfaceRegionName ) const;

  /**
   * @brief adds fracture-matrix connections to the edfm stencil
   * @param mesh the mesh object
   * @param embeddedSurfaceRegionName name of the fracture region
   */
  void addFractureMatrixConnections( MeshLevel & mesh,
                                     string const & embeddedSurfaceRegionName ) const;
private:

  /// mean permeability coefficient
  real64 m_meanPermCoefficient;
  /// flag to determine whether or not to use projection EDFM
  integer m_useProjectionEmbeddedFractureMethod;
};

}


#endif //GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
