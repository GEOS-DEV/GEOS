/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file CoupledSolidBase.cpp
 */

#include "CoupledSolidBase.hpp"

#include "constitutive/solid/SolidFields.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CoupledSolidBase::CoupledSolidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model." );

  registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the porosity model." );

  registerWrapper( viewKeyStruct::permeabilityModelNameString(), &m_permeabilityModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the permeability model." );

  registerWrapper( viewKeyStruct::solidInternalEnergyModelNameString(), &m_solidInternalEnergyModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the solid internal energy model." );

  registerWrapper( viewKeyStruct::solidThermalExpansionModelNameString(), &m_solidThermalExpansionModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the solid thermal expansion model. When omitted, the thermal expansion coefficient "
                    "of the solid model is used (constant in temperature)." );
}

void CoupledSolidBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  if( m_solidThermalExpansionModelName.empty() )
  {
    return;
  }

  GEOS_THROW_IF( !supportsSolidThermalExpansionModel(),
                 GEOS_FMT( "{} is not supported by {}",
                           viewKeyStruct::solidThermalExpansionModelNameString(), getCatalogName() ),
                 InputError, getWrapperDataContext( viewKeyStruct::solidThermalExpansionModelNameString() ) );

  GEOS_THROW_IF( !getParent().hasGroup< SolidThermalExpansion >( m_solidThermalExpansionModelName ),
                 GEOS_FMT( "The SolidThermalExpansion model '{}' could not be found", m_solidThermalExpansionModelName ),
                 InputError, getWrapperDataContext( viewKeyStruct::solidThermalExpansionModelNameString() ) );

  // The coefficient must come from one place only.
  bool const solidHasCoefficient =
    getBaseSolidModel().getWrapperBase( SolidBase::viewKeyStruct::defaultThermalExpansionCoefficientString() ).getSuccessfulReadFromInput();
  GEOS_THROW_IF( solidHasCoefficient,
                 GEOS_FMT( "The thermal expansion coefficient is given both by {} on the solid model '{}' and by the "
                           "SolidThermalExpansion model '{}'; remove one of them",
                           SolidBase::viewKeyStruct::defaultThermalExpansionCoefficientString(),
                           m_solidModelName, m_solidThermalExpansionModelName ),
                 InputError, getWrapperDataContext( viewKeyStruct::solidThermalExpansionModelNameString() ),
                 getBaseSolidModel().getWrapperDataContext( SolidBase::viewKeyStruct::defaultThermalExpansionCoefficientString() ) );

  // With the model, the porosity uses its coefficient too (Biot-Coussy theory), so the porosity coefficient
  // must not be given either.
  string const porosityCoefficientKey = BiotPorosity::viewKeyStruct::defaultThermalExpansionCoefficientString();
  PorosityBase const & porosityModel = getBasePorosityModel();
  bool const porosityHasCoefficient = porosityModel.hasWrapper( porosityCoefficientKey ) &&
                                      porosityModel.getWrapperBase( porosityCoefficientKey ).getSuccessfulReadFromInput();
  GEOS_THROW_IF( porosityHasCoefficient,
                 GEOS_FMT( "{} cannot be given on the porosity model '{}' together with the SolidThermalExpansion model '{}': "
                           "the porosity then uses the coefficient of that model; remove {}",
                           porosityCoefficientKey, m_porosityModelName, m_solidThermalExpansionModelName, porosityCoefficientKey ),
                 InputError, getWrapperDataContext( viewKeyStruct::solidThermalExpansionModelNameString() ),
                 porosityModel.getWrapperDataContext( porosityCoefficientKey ) );
}

SolidThermalExpansionUpdates CoupledSolidBase::createThermalExpansionUpdates() const
{
  if( hasSolidThermalExpansionModel() )
  {
    return getSolidThermalExpansionModel().createKernelUpdates();
  }
  return createThermalExpansionUpdates( getBaseSolidModel() );
}

SolidThermalExpansionUpdates CoupledSolidBase::createThermalExpansionUpdates( SolidBase const & solidModel )
{
  return SolidThermalExpansionUpdates( solidModel.getReference< array1d< real64 > >( fields::solid::thermalExpansionCoefficient::key() ).toViewConst() );
}

} /* namespace constitutive */

} /* namespace geos */
