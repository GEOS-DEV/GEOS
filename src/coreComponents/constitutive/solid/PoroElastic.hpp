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
 * @file PoroElastic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_

#include "SolidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

template< typename BASE >
class PoroElastic : public BASE
{
public:
  PoroElastic( string const & name, dataRepository::Group * const parent );
  virtual ~PoroElastic() override;


  static std::string CatalogName() { return string( "Poro" ) + BASE::m_catalogNameString; }
  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void PostProcessInput() override;

  virtual void
  DeliverClone( string const & name,
                dataRepository::Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  inline virtual void
  StateUpdatePointPressure( real64 const & pres,
                            localIndex const k,
                            localIndex const q ) override
  {
    m_poreVolumeRelation.Compute( pres, m_poreVolumeMultiplier[k][q], m_dPVMult_dPressure[k][q] );
  }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {

    static constexpr auto compressibilityString =  "compressibility";
    static constexpr auto referencePressureString =  "referencePressure";
    static constexpr auto biotCoefficientString =  "BiotCoefficient";
  };


protected:
  /// scalar compressibility parameter
  real64 m_compressibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// scalar Biot's coefficient
  real64 m_biotCoefficient;

  array2d< real64 > m_poreVolumeMultiplier;
  array2d< real64 > m_dPVMult_dPressure;

  ExponentialRelation< real64, ExponentApproximationType::Linear > m_poreVolumeRelation;

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_ */
