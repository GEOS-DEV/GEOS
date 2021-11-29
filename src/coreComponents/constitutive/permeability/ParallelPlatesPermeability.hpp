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
                                    arrayView3d< real64 > const & dPerm_dAper )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dAperture( dPerm_dAper )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dAperture ) const
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
      dPerm_dAperture[dim]  = dPerm;
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
             m_dPerm_dAperture[k][0] );
  }

private:

  arrayView3d< real64 > m_dPerm_dAperture;

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
                          m_dPerm_dAperture );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {} viewKeys;

private:

  array3d< real64 > m_dPerm_dAperture;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
