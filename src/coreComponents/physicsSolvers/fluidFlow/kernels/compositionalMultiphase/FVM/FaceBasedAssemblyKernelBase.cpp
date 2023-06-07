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
 * @file FaceBasedAssemblyKernelBase.cpp
 */

#include "FaceBasedAssemblyKernelBase.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** FaceBasedAssemblyKernel ********************************/

FaceBasedAssemblyKernelBase::FaceBasedAssemblyKernelBase( integer const numPhases,
                                                          globalIndex const rankOffset,
                                                          integer const hasCapPressure,
                                                          DofNumberAccessor const & dofNumberAccessor,
                                                          CompFlowAccessors const & compFlowAccessors,
                                                          MultiFluidAccessors const & multiFluidAccessors,
                                                          CapPressureAccessors const & capPressureAccessors,
                                                          PermeabilityAccessors const & permeabilityAccessors,
                                                          real64 const & dt,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
  : m_numPhases( numPhases ),
  m_rankOffset( rankOffset ),
  m_hasCapPressure( hasCapPressure ),
  m_dt( dt ),
  m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
  m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
  m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
  m_ghostRank( compFlowAccessors.get( fields::ghostRank {} ) ),
  m_gravCoef( compFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
  m_pres( compFlowAccessors.get( fields::flow::pressure {} ) ),
  m_dCompFrac_dCompDens( compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity {} ) ),
  m_dPhaseVolFrac( compFlowAccessors.get( fields::flow::dPhaseVolumeFraction {} ) ),
  m_phaseMob( compFlowAccessors.get( fields::flow::phaseMobility {} ) ),
  m_dPhaseMob( compFlowAccessors.get( fields::flow::dPhaseMobility {} ) ),
  m_phaseMassDens( multiFluidAccessors.get( fields::multifluid::phaseMassDensity {} ) ),
  m_dPhaseMassDens( multiFluidAccessors.get( fields::multifluid::dPhaseMassDensity {} ) ),
  m_phaseCompFrac( multiFluidAccessors.get( fields::multifluid::phaseCompFraction {} ) ),
  m_dPhaseCompFrac( multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction {} ) ),
  m_phaseCapPressure( capPressureAccessors.get( fields::cappres::phaseCapPressure {} ) ),
  m_dPhaseCapPressure_dPhaseVolFrac( capPressureAccessors.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
  m_localMatrix( localMatrix ),
  m_localRhs( localRhs )
{ }

}

}
