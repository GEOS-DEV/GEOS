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
 * @file CompositionalMultiphaseBaseComputeHydrostaticEquilibrium.cpp
 */

#include "CompositionalMultiphaseBase.hpp"

#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace constitutive;

void CompositionalMultiphaseBase::computeHydrostaticEquilibrium()
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  integer const numComps = m_numComponents;
  integer const numPhases = m_numPhases;

  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  // Step 1: count individual equilibriums (there may be multiple ones)

  std::map< string, localIndex > equilNameToEquilId;
  localIndex equilCounter = 0;

  fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
  {

    // collect all the equilibrium names to idx
    equilNameToEquilId[bc.getName()] = equilCounter;
    equilCounter++;

    // check that the gravity vector is aligned with the z-axis
    GEOS_THROW_IF( !isZero( gravVector[0] ) || !isZero( gravVector[1] ),
                   getCatalogName() << " " << getDataContext() <<
                   ": the gravity vector specified in this simulation (" << gravVector[0] << " " << gravVector[1] << " " << gravVector[2] <<
                   ") is not aligned with the z-axis. \n"
                   "This is incompatible with the " << EquilibriumInitialCondition::catalogName() << " " << bc.getDataContext() <<
                   "used in this simulation. To proceed, you can either: \n" <<
                   "   - Use a gravityVector aligned with the z-axis, such as (0.0,0.0,-9.81)\n" <<
                   "   - Remove the hydrostatic equilibrium initial condition from the XML file",
                   InputError );

  } );

  if( equilCounter == 0 )
  {
    return;
  }

  // Step 2: find the min elevation and the max elevation in the targetSets

  array1d< real64 > globalMaxElevation( equilNameToEquilId.size() );
  array1d< real64 > globalMinElevation( equilNameToEquilId.size() );
  findMinMaxElevationInEquilibriumTarget( domain,
                                          equilNameToEquilId,
                                          globalMaxElevation,
                                          globalMinElevation );

  // Step 3: for each equil, compute a fine table with hydrostatic pressure vs elevation if the region is a target region

  // first compute the region filter
  std::set< string > regionFilter;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    for( string const & regionName : regionNames )
    {
      regionFilter.insert( regionName );
    }

    fsManager.apply< ElementSubRegionBase,
                     EquilibriumInitialCondition >( 0.0,
                                                    mesh,
                                                    EquilibriumInitialCondition::catalogName(),
                                                    [&] ( EquilibriumInitialCondition const & fs,
                                                          string const &,
                                                          SortedArrayView< localIndex const > const & targetSet,
                                                          ElementSubRegionBase & subRegion,
                                                          string const & )
    {
      // Step 3.1: retrieve the data necessary to construct the pressure table in this subregion

      integer const maxNumEquilIterations = fs.getMaxNumEquilibrationIterations();
      real64 const equilTolerance = fs.getEquilibrationTolerance();
      real64 const datumElevation = fs.getDatumElevation();
      real64 const datumPressure = fs.getDatumPressure();
      string const initPhaseName = fs.getInitPhaseName(); // will go away when GOC/WOC are implemented

      localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
      real64 const minElevation = LvArray::math::min( globalMinElevation[equilIndex], datumElevation );
      real64 const maxElevation = LvArray::math::max( globalMaxElevation[equilIndex], datumElevation );
      real64 const elevationIncrement = LvArray::math::min( fs.getElevationIncrement(), maxElevation - minElevation );
      localIndex const numPointsInTable = ( elevationIncrement > 0 ) ? std::ceil( (maxElevation - minElevation) / elevationIncrement ) + 1 : 1;

      real64 const eps = 0.1 * (maxElevation - minElevation); // we add a small buffer to only log in the pathological cases
      GEOS_LOG_RANK_0_IF( ( (datumElevation > globalMaxElevation[equilIndex]+eps)  || (datumElevation < globalMinElevation[equilIndex]-eps) ),
                          getCatalogName() << " " << getDataContext() <<
                          ": By looking at the elevation of the cell centers in this model, GEOS found that " <<
                          "the min elevation is " << globalMinElevation[equilIndex] << " and the max elevation is " <<
                          globalMaxElevation[equilIndex] << "\nBut, a datum elevation of " << datumElevation <<
                          " was specified in the input file to equilibrate the model.\n " <<
                          "The simulation is going to proceed with this out-of-bound datum elevation," <<
                          " but the initial condition may be inaccurate." );

      array1d< array1d< real64 > > elevationValues;
      array1d< real64 > pressureValues;
      elevationValues.resize( 1 );
      elevationValues[0].resize( numPointsInTable );
      pressureValues.resize( numPointsInTable );

      // Step 3.2: retrieve the user-defined tables (temperature and comp fraction)

      FunctionManager & functionManager = FunctionManager::getInstance();

      array1d< TableFunction::KernelWrapper > compFracTableWrappers;
      arrayView1d< string const > compFracTableNames = fs.getComponentFractionVsElevationTableNames();
      for( integer ic = 0; ic < numComps; ++ic )
      {
        TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( compFracTableNames[ic] );
        compFracTableWrappers.emplace_back( compFracTable.createKernelWrapper() );
      }

      string const tempTableName = fs.getTemperatureVsElevationTableName();
      TableFunction const & tempTable = functionManager.getGroup< TableFunction >( tempTableName );
      TableFunction::KernelWrapper tempTableWrapper = tempTable.createKernelWrapper();

      // Step 3.3: retrieve the fluid model to compute densities
      // we end up with the same issue as in applyDirichletBC: there is not a clean way to retrieve the fluid info

      Group const & region = subRegion.getParent().getParent();
      auto itRegionFilter = regionFilter.find( region.getName() );
      if( itRegionFilter == regionFilter.end() )
      {
        return; // the region is not in target, there is nothing to do
      }
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      arrayView1d< string const > componentNames = fs.getComponentNames();
      GEOS_THROW_IF( fluid.componentNames().size() != componentNames.size(),
                     "Mismatch in number of components between constitutive model "
                     << fluid.getDataContext() << " and the Equilibrium initial condition " << fs.getDataContext(),
                     InputError );
      for( integer ic = 0; ic < fluid.numFluidComponents(); ++ic )
      {
        GEOS_THROW_IF( fluid.componentNames()[ic] != componentNames[ic],
                       "Mismatch in component names between constitutive model "
                       << fluid.getDataContext() << " and the Equilibrium initial condition " << fs.getDataContext(),
                       InputError );
      }

      // Note: for now, we assume that the reservoir is in a single-phase state at initialization
      arrayView1d< string const > phaseNames = fluid.phaseNames();
      auto const itPhaseNames = std::find( std::begin( phaseNames ), std::end( phaseNames ), initPhaseName );
      GEOS_THROW_IF( itPhaseNames == std::end( phaseNames ),
                     getCatalogName() << " " << getDataContext() << ": phase name " <<
                     initPhaseName << " not found in the phases of " << fluid.getDataContext(),
                     InputError );
      integer const ipInit = std::distance( std::begin( phaseNames ), itPhaseNames );

      // Step 3.4: compute the hydrostatic pressure values

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        using FluidType = TYPEOFREF( castedFluid );
        typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

        // note: inside this kernel, serialPolicy is used, and elevation/pressure values don't go to the GPU
        isothermalCompositionalMultiphaseBaseKernels::
          HydrostaticPressureKernel::ReturnType const returnValue =
          isothermalCompositionalMultiphaseBaseKernels::
            HydrostaticPressureKernel::launch( numPointsInTable,
                                               numComps,
                                               numPhases,
                                               ipInit,
                                               maxNumEquilIterations,
                                               equilTolerance,
                                               gravVector,
                                               minElevation,
                                               elevationIncrement,
                                               datumElevation,
                                               datumPressure,
                                               fluidWrapper,
                                               compFracTableWrappers.toViewConst(),
                                               tempTableWrapper,
                                               elevationValues.toNestedView(),
                                               pressureValues.toView() );

        GEOS_THROW_IF( returnValue ==  isothermalCompositionalMultiphaseBaseKernels::HydrostaticPressureKernel::ReturnType::FAILED_TO_CONVERGE,
                       getCatalogName() << " " << getDataContext() <<
                       ": hydrostatic pressure initialization failed to converge in region " << region.getName() << "! \n" <<
                       "Try to loosen the equilibration tolerance, or increase the number of equilibration iterations. \n" <<
                       "If nothing works, something may be wrong in the fluid model, see <Constitutive> ",
                       std::runtime_error );

        GEOS_LOG_RANK_0_IF( returnValue == isothermalCompositionalMultiphaseBaseKernels::HydrostaticPressureKernel::ReturnType::DETECTED_MULTIPHASE_FLOW,
                            getCatalogName() << " " << getDataContext() <<
                            ": currently, GEOS assumes that there is only one mobile phase when computing the hydrostatic pressure. \n" <<
                            "We detected multiple phases using the provided datum pressure, temperature, and component fractions. \n" <<
                            "Please make sure that only one phase is mobile at the beginning of the simulation. \n" <<
                            "If this is not the case, the problem will not be at equilibrium when the simulation starts" );

      } );

      // Step 3.5: create hydrostatic pressure table

      string const tableName = fs.getName() + "_" + subRegion.getName() + "_" + phaseNames[ipInit] + "_table";
      TableFunction * const presTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
      presTable->setTableCoordinates( elevationValues, { units::Distance } );
      presTable->setTableValues( pressureValues, units::Pressure );
      presTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
      TableFunction::KernelWrapper presTableWrapper = presTable->createKernelWrapper();

      // Step 4: assign pressure, temperature, and component fraction as a function of elevation
      // TODO: this last step should probably be delayed to wait for the creation of FaceElements
      // TODO: this last step should be modified to account for GOC and WOC
      arrayView2d< real64 const > const elemCenter =
        subRegion.getReference< array2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

      arrayView1d< real64 > const pres = subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView1d< real64 > const temp = subRegion.getReference< array1d< real64 > >( fields::flow::temperature::key() );
      arrayView2d< real64, compflow::USD_COMP > const compFrac =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::globalCompFraction::key() );
      arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappersViewConst =
        compFracTableWrappers.toViewConst();

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minPressure( LvArray::NumericLimits< real64 >::max );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [targetSet,
                                                           elemCenter,
                                                           presTableWrapper,
                                                           tempTableWrapper,
                                                           compFracTableWrappersViewConst,
                                                           numComps,
                                                           minPressure,
                                                           pres,
                                                           temp,
                                                           compFrac] GEOS_HOST_DEVICE ( localIndex const i )
      {
        localIndex const k = targetSet[i];
        real64 const elevation = elemCenter[k][2];

        pres[k] = presTableWrapper.compute( &elevation );
        minPressure.min( pres[k] );
        temp[k] = tempTableWrapper.compute( &elevation );
        for( integer ic = 0; ic < numComps; ++ic )
        {
          compFrac[k][ic] = compFracTableWrappersViewConst[ic].compute( &elevation );
        }
      } );

      GEOS_ERROR_IF( minPressure.get() < 0.0,
                     GEOS_FMT( "{}: A negative pressure of {} Pa was found during hydrostatic initialization in region/subRegion {}/{}",
                               getDataContext(), minPressure.get(), region.getName(), subRegion.getName() ) );
    } );
  } );
}

} // namespace geos
