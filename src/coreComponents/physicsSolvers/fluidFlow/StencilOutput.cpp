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
 * @file StencilOutput.cpp
 */

#include "StencilOutput.hpp"

#include "common/Units.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

StencilOutput::StencilOutput( const string & name,
                              Group * const parent ):
  Base( name, parent )
{}


class StencilOutput::StencilOutputKernel
{
public:
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using PermeabilityAccessors = StencilMaterialAccessors< PermeabilityBase,
                                                          fields::permeability::permeability >;


  /**
   * @brief launch a kernel to compute the element-element connection data using a given StencilWrapper.
   * @tparam POLICY the Raja host / device launching policy
   * @tparam STENCILWRAPPER_T
   * @param connData The buffer where the connection data will be written
   * @param stencilWrapper the StencilWrapper to use to compute the element-element connection data
   * @param permeability the permeability data that will be used to compute the transmissibility
   */
  template< typename POLICY, typename STENCILWRAPPER_T >
  static void launch( arrayView1d< StencilOutput::KernelConnectionData > const & connData,
                      STENCILWRAPPER_T const & stencilWrapper,
                      ElementViewConst< arrayView3d< real64 const > > const & permeability )
  {
    forAll< POLICY >( stencilWrapper.size(), [stencilWrapper,
                                              permeability,
                                              connData] GEOS_HOST_DEVICE ( localIndex const iConn )
    {
      real64 transmissibility[1][2];
      real64 dummy[1][2];

      stencilWrapper.computeWeights( iConn,
                                     permeability,
                                     permeability,      // we don't need derivatives
                                     transmissibility,
                                     dummy );           // we don't need derivatives

      for( localIndex i = 0; i < 2; ++i )
      {
        connData[iConn].transmissibility[i] = transmissibility[0][i];
        connData[iConn].regionId[i] =    stencilWrapper.getElementRegionIndices()( iConn, i );
        connData[iConn].subRegionId[i] = stencilWrapper.getElementSubRegionIndices()( iConn, i );
        connData[iConn].elementId[i] =   stencilWrapper.getElementIndices()( iConn, i );
      }
    } );
  }

private:
  StencilOutputKernel()
  {}
};

StencilOutput::ConnectionData StencilOutput::ToConnectionData( KernelConnectionData const & kernelData,
                                                               LocalToGlobalMap const & localToGlobalMap )
{
  ConnectionData conn;
  for( localIndex i = 0; i < 2; ++i )
  {
    conn.transmissibility[i] = kernelData.transmissibility[i];
    conn.globalId[i] = localToGlobalMap[kernelData.regionId[i]][kernelData.subRegionId[i]][kernelData.elementId[i]];
  }
  return conn;
}

template< typename STENCILWRAPPER_T >
array1d< StencilOutput::KernelConnectionData >
StencilOutput::gatherTimeStepData( MeshLevel const & mesh,
                                   STENCILWRAPPER_T const & stencilWrapper ) const
{
  ElementRegionManager const & elemManager = mesh.getElemManager();

  // allocate a big enough buffer to store all connection data
  array1d< KernelConnectionData > kernelData{ stencilWrapper.size() };
  typename StencilOutputKernel::PermeabilityAccessors accessor( elemManager, m_solver->getName() );

  StencilOutputKernel::launch< parallelDevicePolicy<> >( kernelData.toView(), stencilWrapper,
                                                         accessor.get<
                                                           fields::permeability::permeability >() );

  return kernelData;
}


string StencilOutput::getOutputFileName( string_view stencilName ) const
{
  return GEOS_FMT( "{}/{}.csv",
                   m_outputDir, stencilName );
}

void writeHeader( std::ofstream & outputFile )
{
  outputFile << GEOS_FMT( "{},A global id,B global id,{},{}",
                          units::getDescription( units::Unit::Time ),
                          GEOS_FMT( "A {}", units::getDescription( units::Unit::Transmissibility ) ),
                          GEOS_FMT( "B {}", units::getDescription( units::Unit::Transmissibility ) ) );
  outputFile << std::endl;
}

template< typename CONN_DATA_CONT >
void writeData( std::ofstream & outputFile, real64 const time, CONN_DATA_CONT const & data )
{
  for( StencilOutput::ConnectionData const & conn : data )
  {
    outputFile << time << ","
               << conn.globalId[0] << ","
               << conn.globalId[1] << ","
               << conn.transmissibility[0] << ","
               << conn.transmissibility[1] << std::endl;
  }
}

void StencilOutput::outputTimeStepData( MeshLevel const & mesh,
                                        string_view stencilName,
                                        real64 const outputTime,
                                        arrayView1d< KernelConnectionData > const & kernelData )
{
  string const fileName = getOutputFileName( stencilName );
  bool const newFile = outputFiles.count( fileName ) == 0;
  std::ofstream outputFile( fileName, newFile ? std::ios_base::out : std::ios_base::app );

  if( newFile )
  {
    writeHeader( outputFile );
    outputFiles.insert( fileName );
    GEOS_LOG_RANK_IF( getLogLevel() > 0, GEOS_FMT( "{}: initialized file \"{}\".",
                                                   getDataContext(), fileName ) );
  }

  { // output timestep data
    std::set< ConnectionData > fileData;

    LocalToGlobalMap localToGlobalMap = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >(
      ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );

    std::transform( kernelData.begin(), kernelData.end(), std::inserter( fileData, fileData.end() ),
                    [&]( auto const & kernelConn )
    {
      return ToConnectionData( kernelConn, localToGlobalMap );
    } );

    writeData( outputFile, outputTime, fileData );
    GEOS_LOG_RANK_IF( getLogLevel() > 0, GEOS_FMT( "{}: {} connections writen in file \"{}\".",
                                                   getDataContext(), kernelData.size(), fileName ) );
  }

  outputFile.close();
}


bool StencilOutput::execute( real64 const time_n,
                             real64 const dt,
                             integer const GEOS_UNUSED_PARAM( cycleNumber ),
                             integer const GEOS_UNUSED_PARAM( eventCounter ),
                             real64 const GEOS_UNUSED_PARAM( eventProgress ),
                             DomainPartition & domain )
{
  if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
  {
    m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                            MeshLevel & mesh,
                                                                            arrayView1d< string const > const & )
    {
      NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
      FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
      FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_solver->getDiscretizationName() );

      fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
      {
        // gather
        auto const stencilWrapper = stencil.createKernelWrapper();
        array1d< KernelConnectionData > const kernelData = gatherTimeStepData( mesh, stencilWrapper );
        // output
        outputTimeStepData( mesh, fluxApprox.getName(), time_n + dt, kernelData.toView() );
      } );
    } );
  }
  return false;
}


// bool StencilOutput::execute( real64 const time_n,
//                              real64 const dt,
//                              integer const GEOS_UNUSED_PARAM( cycleNumber ),
//                              integer const GEOS_UNUSED_PARAM( eventCounter ),
//                              real64 const GEOS_UNUSED_PARAM( eventProgress ),
//                              DomainPartition & domain )
// {

//   if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
//   {
//     m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
//                                                                             MeshLevel & mesh,
//                                                                             arrayView1d< string const > const & )
//     {
//       array1d< KernelConnectionData > kernelData;
//       ElementRegionManager const & elemManager = mesh.getElemManager();

//       { // data aquisition
//         NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
//         FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
//         FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_solver->getDiscretizationName() );

//         {
//           int connectionCount = 0;
//           fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
//           {
//             connectionCount += stencil.size();
//           } );
//           kernelData.resize( connectionCount );
//         }

//         fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
//         {
//           auto const stencilWrapper = stencil.createKernelWrapper();
//           using StencilOutputKernel;
//           typename PermeabilityAccessors accessor( elemManager, m_solver->getName() );

//           launch< parallelDevicePolicy<> >( stencilWrapper,
//                                             accessor.get< fields::permeability::permeability >(),
//                                             kernelData.toView() );
//         } );
//       }


//       { // write output file
//         string const fileName = getOutputFileName( mesh.getName() );
//         real64 const outputTime = time_n + dt;
//         bool const newFile = outputFiles.count( fileName ) == 0;
//         std::ofstream outputFile( fileName, newFile ? std::ios_base::out : std::ios_base::app );

//         if( newFile )
//         {
//           writeHeader( outputFile );
//           outputFiles.insert( fileName );
//           GEOS_LOG_RANK_IF( getLogLevel() > 0, GEOS_FMT( "{}: initialized file \"{}\".",
//                                                          getDataContext(), fileName ) );

//         }

//         { // output timestep data
//           LocalToGlobalMap localToGlobalMap =
//             elemManager.constructArrayViewAccessor< globalIndex, 1 >(
//               ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );
//           std::set< ConnectionData > fileData;

//           std::transform( kernelData.begin(), kernelData.end(), std::inserter( fileData, fileData.end() ),
//                           [&]( auto const & kernelConn )
//           {
//             return ToConnectionData( kernelConn, localToGlobalMap );
//           } );

//           writeData( outputFile, outputTime, fileData );
//           GEOS_LOG_RANK_IF( getLogLevel() > 0, GEOS_FMT( "{}: {} connections writen in file \"{}\".",
//                                                          getDataContext(), kernelData.size(), fileName ) );
//         }


//         outputFile.close();
//       }
//     } );
//   }
//   return false;
// }


REGISTER_CATALOG_ENTRY( TaskBase,
                        StencilOutput,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
