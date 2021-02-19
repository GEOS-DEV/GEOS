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
 * @file FracturePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_

#include "constitutive/permeability/FracturePermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class ParallelPlatesPermeabilityUpdate : public FracturePermeabilityBaseUpdate
{
public:

  ParallelPlatesPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dAperture )
    : FracturePermeabilityBaseUpdate( permeability, dPerm_dAperture )
  {}

  /// Default copy constructor
  ParallelPlatesPermeabilityUpdate( ParallelPlatesPermeabilityUpdate const & ) = default;

  /// Default move constructor
  ParallelPlatesPermeabilityUpdate( ParallelPlatesPermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  ParallelPlatesPermeabilityUpdate & operator=( ParallelPlatesPermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  ParallelPlatesPermeabilityUpdate & operator=( ParallelPlatesPermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const & effectiveAperture,
                        arraySlice1d< real64 > const & permeability,
                        arraySlice1d< real64 > const & dPerm_dAperture ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & effectiveAperture ) const
  {
    compute( effectiveAperture,
             m_permeability[k][q],
             m_dPerm_dAperture[k][q] );
  }

};


class ParallelPlatesPermeability : public FracturePermeabilityBase
{
public:
  ParallelPlatesPermeability( string const & name, Group * const parent );

  virtual ~ParallelPlatesPermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ParallelPlatesPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ParallelPlatesPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dAperture );
  }


  struct viewKeyStruct : public FracturePermeabilityBase::viewKeyStruct
  {} viewKeys;

protected:

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ParallelPlatesPermeabilityUpdate::compute( real64 const & effectiveAperture,
                                                arraySlice1d< real64 > const & permeability,
                                                arraySlice1d< real64 > const & dPerm_dAperture ) const
{
  permeability[0] = effectiveAperture*effectiveAperture*effectiveAperture /12;
  dPerm_dAperture[0]  = effectiveAperture*effectiveAperture / 12;
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
