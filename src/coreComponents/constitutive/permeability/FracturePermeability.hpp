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
 * @file PoreVolumeCompressibleSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class FracturePermeabilityUpdate
{
public:

  FracturePermeabilityUpdate( arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & dPerm_dAperture )
    : m_permeability( permeability ),
    m_dPerm_dAperture( dPerm_dAperture )
  {}

  /// Default copy constructor
  FracturePermeabilityUpdate( FracturePermeabilityUpdate const & ) = default;

  /// Default move constructor
  FracturePermeabilityUpdate( FracturePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  FracturePermeabilityUpdate & operator=( FracturePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  FracturePermeabilityUpdate & operator=( FracturePermeabilityUpdate && ) = delete;

protected:

  arrayView3d< real64 > m_permeability;
  arrayView3d< real64 > m_dPerm_dAperture;

private:

  void compute( real64 const & effectiveAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dAperture );

  void update(  );
};


class FracturePermeability : public PermeabilityBase
{
public:
  FracturePermeability( string const & name, Group * const parent );

  virtual ~FracturePermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "FracturePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FracturePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dAperture );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {} viewKeys;

protected:
  virtual void postProcessInput() override;

private:

  array3d< real64 > m_dPerm_dAperture;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
