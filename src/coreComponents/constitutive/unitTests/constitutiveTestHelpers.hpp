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

#pragma once

// Source includes
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/relativePermeability/relativePermeabilitySelector.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

namespace geosx
{
namespace testing
{

template< typename MODEL, typename VAR, typename D_VAR_D_SAT >
void testNumericalDerivatives( dataRepository::Group & parent,
                               MODEL & model,
                               arraySlice1d< real64 const > const saturation,
                               real64 const perturbParameter,
                               real64 const relTol,
                               std::string const & varName,
                               VAR && varAccessor,
                               D_VAR_D_SAT && dVar_dSat_accessor )
{
  localIndex const NP = model.numFluidPhases();
  auto const & phases = model.phaseNames();

  // create a clone of the rel perm to run updates on
  std::unique_ptr< constitutive::ConstitutiveBase > modelCopyPtr;
  model.DeliverClone( "fluidCopy", &parent, modelCopyPtr );
  MODEL & modelCopy = *modelCopyPtr->group_cast< MODEL * >();

  model.AllocateConstitutiveData( model.getParent(), 1 );
  modelCopy.AllocateConstitutiveData( model.getParent(), 1 );

  arraySlice1d< real64 const > const var = varAccessor( model );
  arraySlice2d< real64 const > const dPhaseRelPerm_dSat = dVar_dSat_accessor( model );
  arraySlice1d< real64 const > const varCopy = varAccessor( modelCopy );

  // set the fluid state to current
  constitutive::constitutiveUpdatePassThru( model, [&] ( auto & castedModel )
  {
    typename TYPEOFREF( castedModel ) ::KernelWrapper relPermWrapper = castedModel.createKernelWrapper();
    relPermWrapper.Update( 0, 0, saturation );
  } );

  // update saturation and check derivatives
  auto dPhaseRelPerm_dS = testing::invertLayout( dPhaseRelPerm_dSat, NP, NP );

  array1d< real64 > satNew( NP );
  for( localIndex jp = 0; jp < NP; ++jp )
  {
    real64 const dS = perturbParameter * (saturation[jp] + perturbParameter);
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      satNew[ip] = saturation[ip];
    }
    satNew[jp] += dS;

    constitutive::constitutiveUpdatePassThru( modelCopy, [&] ( auto & castedRelPerm )
    {
      typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();
      relPermWrapper.Update( 0, 0, satNew );
    } );

    checkDerivative( varCopy.toSliceConst(),
                     var.toSliceConst(),
                     dPhaseRelPerm_dS[jp].toSliceConst(),
                     dS,
                     relTol,
                     varName,
                     "phaseVolFrac[" + phases[jp] + "]",
                     phases );
  }
}

template< typename BASE >
class ConstitutiveTestBase : public ::testing::Test
{
public:
  ConstitutiveTestBase():
    m_node(),
    m_parent( "parent", m_node )
  {
    m_parent.resize( 1 );
  }

  void initialize( BASE * model )
  {
    m_model = model;
    m_parent.Initialize( &m_parent );
    m_parent.InitializePostInitialConditions( &m_parent );
  }

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  BASE * m_model;
};

} // namespace testing
} // namespace geosx
