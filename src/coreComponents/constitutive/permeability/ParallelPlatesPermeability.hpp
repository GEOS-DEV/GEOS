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
 * @file ParallelPlatesPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class ParallelPlatesPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ParallelPlatesPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dPressure,
                                    arrayView4d< real64 > const & dPerm_dDispJump )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dDispJump( dPerm_dDispJump )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice2d< real64 > const & dPerm_dDispJump ) const
  {
    // TODO: maybe move to this computation or have the possibility of choosing.
//    real64 const perm  = newHydraulicAperture*newHydraulicAperture*newHydraulicAperture / 12.0;
//    real64 const dPerm = newHydraulicAperture*newHydraulicAperture / 4.0;

    real64 const perm = 0.25 * ( oldHydraulicAperture*oldHydraulicAperture*oldHydraulicAperture +
                                 oldHydraulicAperture*oldHydraulicAperture*newHydraulicAperture +
                                 oldHydraulicAperture*newHydraulicAperture*newHydraulicAperture +
                                 newHydraulicAperture*newHydraulicAperture*newHydraulicAperture ) / 12;

    real64 const dPerm  = 0.25 * ( oldHydraulicAperture*oldHydraulicAperture +
                                   2*oldHydraulicAperture*newHydraulicAperture +
                                   3*newHydraulicAperture*newHydraulicAperture ) / 12;


    for( int dim=0; dim < 3; dim++ )
    {
      permeability[dim]     = perm;
      dPerm_dDispJump[dim][0]  = dPerm;
      dPerm_dDispJump[dim][1]  = 0.0;
      dPerm_dDispJump[dim][2]  = 0.0;
    }
  }

  GEOSX_HOST_DEVICE
  virtual void updateFromAperture( localIndex const k,
                                   localIndex const q,
                                   real64 const & oldHydraulicAperture,
                                   real64 const & newHydraulicAperture ) const override final
  {
    GEOSX_UNUSED_VAR( q );

    compute( oldHydraulicAperture,
             newHydraulicAperture,
             m_permeability[k][0],
             m_dPerm_dDispJump[k][0] );
  }

  GEOSX_HOST_DEVICE
  virtual void updateFromApertureAndShearDisplacement( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const & pressure,
                                                       real64 const ( &dispJump )[3],
                                                       real64 const ( &traction )[3] ) const override final
  {
    GEOSX_UNUSED_VAR( dispJump, traction, pressure );

    updateFromAperture( k, q, oldHydraulicAperture, newHydraulicAperture );
  }

private:

  arrayView4d< real64 > m_dPerm_dDispJump;

};


class ParallelPlatesPermeability : public PermeabilityBase
{
public:

  ParallelPlatesPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ParallelPlatesPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ParallelPlatesPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dDispJump );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {} viewKeys;

private:

  array3d< real64 > m_dPerm_dAperture;

  /// Derivative of fracture permeability w.r.t. displacement jump
  array4d< real64 > m_dPerm_dDispJump;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
