/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file ThermalCompressibleSolid.cpp
 */

#include "ThermalCompressibleSolid.hpp"
#include "porosity/PressurePorosity.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"
#include "constitutive/permeability/ExponentialDecayPermeability.hpp"
#include "constitutive/permeability/ParallelPlatesPermeability.hpp"
#include "constitutive/permeability/PressurePermeability.hpp"
#include "constitutive/permeability/SlipDependentPermeability.hpp"
#include "constitutive/permeability/WillisRichardsPermeability.hpp"
#include "constitutive/thermalConductivity/SinglePhaseThermalConductivity.hpp"
#include "constitutive/thermalConductivity/MultiPhaseConstantThermalConductivity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename THERMAL_COND_TYPE >
ThermalCompressibleSolid< PORO_TYPE, PERM_TYPE, THERMAL_COND_TYPE >::ThermalCompressibleSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >( name, parent )
{}

template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename THERMAL_COND_TYPE >
ThermalCompressibleSolid< PORO_TYPE, PERM_TYPE, THERMAL_COND_TYPE >::~ThermalCompressibleSolid() = default;

// Register all ThermalCompressibleSolid model types with single phase flow.
typedef ThermalCompressibleSolid< PressurePorosity, ConstantPermeability, SinglePhaseThermalConductivity > CompressibleRockConstantSinglePhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, CarmanKozenyPermeability, SinglePhaseThermalConductivity > CompressibleRockCKSinglePhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, PressurePermeability, SinglePhaseThermalConductivity > CompressibleRockPressurePermSinglePhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, ExponentialDecayPermeability, SinglePhaseThermalConductivity > FaultEDSinglePhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, ParallelPlatesPermeability, SinglePhaseThermalConductivity > FractureRockSinglePhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, SlipDependentPermeability, SinglePhaseThermalConductivity > FaultSinglePhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, WillisRichardsPermeability, SinglePhaseThermalConductivity > FaultWRSinglePhaseThermalCond;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockConstantSinglePhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockCKSinglePhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockPressurePermSinglePhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FractureRockSinglePhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultEDSinglePhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultSinglePhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultWRSinglePhaseThermalCond, string const &, Group * const )

// Register all ThermalCompressibleSolid model types with multi phase flow.
typedef ThermalCompressibleSolid< PressurePorosity, ConstantPermeability, MultiPhaseConstantThermalConductivity > CompressibleRockConstantMultiPhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, CarmanKozenyPermeability, MultiPhaseConstantThermalConductivity > CompressibleRockCKMultiPhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, PressurePermeability, MultiPhaseConstantThermalConductivity > CompressibleRockPressurePermMultiPhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, ExponentialDecayPermeability, MultiPhaseConstantThermalConductivity > FaultEDMultiPhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, ParallelPlatesPermeability, MultiPhaseConstantThermalConductivity > FractureRockMultiPhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, SlipDependentPermeability, MultiPhaseConstantThermalConductivity > FaultMultiPhaseThermalCond;
typedef ThermalCompressibleSolid< PressurePorosity, WillisRichardsPermeability, MultiPhaseConstantThermalConductivity > FaultWRMultiPhaseThermalCond;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockConstantMultiPhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockCKMultiPhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockPressurePermMultiPhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FractureRockMultiPhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultEDMultiPhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultMultiPhaseThermalCond, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultWRMultiPhaseThermalCond, string const &, Group * const )
}
} /* namespace geos */
