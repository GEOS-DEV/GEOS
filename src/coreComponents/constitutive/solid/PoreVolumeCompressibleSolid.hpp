/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PoreVolumeCompressibleSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POREVOLUMECOMPRESSIBLESOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POREVOLUMECOMPRESSIBLESOLID_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{


class PoreVolumeCompressibleSolid : public ConstitutiveBase
{
public:
  PoreVolumeCompressibleSolid( std::string const & name, Group * const parent );

  virtual ~PoreVolumeCompressibleSolid() override;

  void DeliverClone( string const & name,
                     Group * const parent,
                     std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  static std::string CatalogName() { return "PoreVolumeCompressibleSolid"; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdatePointPressure( real64 const & pres,
                                         localIndex const k,
                                         localIndex const q ) override final;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    dataRepository::ViewKey compressibility   = { "compressibility"   };
    dataRepository::ViewKey referencePressure = { "referencePressure" };
  } viewKeys;

protected:
  virtual void PostProcessInput() override;

private:

  /// scalar compressibility parameter
  real64 m_compressibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  array2d< real64 > m_poreVolumeMultiplier;
  array2d< real64 > m_dPVMult_dPressure;

  ExponentialRelation< real64, ExponentApproximationType::Linear > m_poreVolumeRelation;
};

inline void PoreVolumeCompressibleSolid::StateUpdatePointPressure( real64 const & pres,
                                                                   localIndex const k,
                                                                   localIndex const q )
{
  m_poreVolumeRelation.Compute( pres, m_poreVolumeMultiplier[k][q], m_dPVMult_dPressure[k][q] );
}

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_SOLID_POREVOLUMECOMPRESSIBLESOLID_HPP_
