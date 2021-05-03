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
 * @file TriaxialDriver.hpp
 */

#ifndef SRC_COMPONENTS_TASKS_TRIAXIALDRIVER_HPP_
#define SRC_COMPONENTS_TASKS_TRIAXIALDRIVER_HPP_

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/GeosxState.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/Functions/TableFunction.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/Tasks/TaskBase.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"

namespace geosx
{

using namespace constitutive;

/**
 * @class TriaxialDriver
 *
 * Class to allow for triaxial tests of the solid constitutive models without the
 * complexity of setting up a single element test.
 *
 */
class TriaxialDriver : public TaskBase
{
public:
  TriaxialDriver( const string & name,
                  Group * const parent );
  ~TriaxialDriver() override;

  static string catalogName() { return "TriaxialDriver"; }

  virtual bool execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                        real64 const GEOSX_UNUSED_PARAM( dt ),
                        integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                        DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override
  {
    GEOSX_THROW_IF( MpiWrapper::commRank() > 0,
                    "Triaxial Driver should only be run in serial",
                    std::runtime_error );

    DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
    ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
    SolidBase & solid = constitutiveManager.getGroup< SolidBase >( m_solidMaterialName );

    GEOSX_LOG_RANK_0( "Launching Triaxial Driver" );
    GEOSX_LOG_RANK_0( "  Material .......... " << m_solidMaterialName );
    GEOSX_LOG_RANK_0( "  Type .............. " << solid.getCatalogName() );
    GEOSX_LOG_RANK_0( "  Mode .............. " << m_mode );
    GEOSX_LOG_RANK_0( "  Strain Function ... " << m_strainFunctionName );
    GEOSX_LOG_RANK_0( "  Stress Function ... " << m_stressFunctionName );
    GEOSX_LOG_RANK_0( "  Steps ............. " << m_numSteps );
    GEOSX_LOG_RANK_0( "  Output ............ " << m_outputFileName );

    conduit::Node node;
    dataRepository::Group rootGroup( "root", node );
    dataRepository::Group discretization( "discretization", &rootGroup );

    discretization.resize( 1 ); // one element
    solid.allocateConstitutiveData( discretization, 1 ); // one quadrature point

    ConstitutivePassThru< SolidBase >::execute( solid, [&]( auto & castedSolid )
    {
      using CONSTITUTIVE_TYPE = TYPEOFREF( castedSolid );
      typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = castedSolid.createKernelUpdates();

      constitutiveUpdate.m_oldStress( 0, 0, 0 ) = m_axialStress[0];
      constitutiveUpdate.m_oldStress( 0, 0, 1 ) = m_radialStress[0];
      constitutiveUpdate.m_oldStress( 0, 0, 2 ) = m_radialStress[0];

      real64 stress[6] = {};
      real64 strainIncrement[6] = {};
      real64 stiffness[6][6] = {{}};

      if( m_mode == "triaxial" )
      {
        printf( "Timestep | Newton | Residual\n" );
        for( localIndex n=1; n<m_time.size(); ++n )
        {
          strainIncrement[0] = m_axialStrain[n]-m_axialStrain[n-1];
          strainIncrement[1] = 0;
          strainIncrement[2] = 0;

          printf( "\n" );
          for( localIndex k=0; k<25; ++k )
          {
            constitutiveUpdate.smallStrainUpdate( 0, 0, strainIncrement, stress, stiffness );
            real64 norm = fabs( stress[1]-m_radialStress[n] )/(fabs( m_radialStress[n] )+1);

            printf( "%8ld   %6ld   %.2e\n", n, k, norm );
            if( norm < 1e-12 )
            {
              break;
            }
            else
            {
              real64 jacobian = stiffness[1][1]+stiffness[1][2];
              strainIncrement[1] -= (stress[1]-m_radialStress[n]) / jacobian;
              strainIncrement[2] = strainIncrement[1];
            }
          }
          castedSolid.saveConvergedState();

          m_axialStress[n]  = stress[0];
          m_radialStrain[n] = m_radialStrain[n-1]+strainIncrement[1];
        }
      }
      else if( m_mode == "volumetric" || m_mode == "oedometer" )
      {
        for( localIndex n=1; n<m_time.size(); ++n )
        {
          strainIncrement[0] = m_axialStrain[n]-m_axialStrain[n-1];
          strainIncrement[1] = m_radialStrain[n]-m_radialStrain[n-1];
          strainIncrement[2] = strainIncrement[1];

          constitutiveUpdate.smallStrainUpdate( 0, 0, strainIncrement, stress, stiffness );
          castedSolid.saveConvergedState();

        }
      }
      else
      {
        GEOSX_THROW( "Test mode \'" << m_mode << "\' not recognized.", InputError );
      }

    } );

    outputResults();
    return false;
  }

  void postProcessInput() override;
  void outputResults();

private:

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * solidMaterialNameString() { return "material"; }
    constexpr static char const * modeString() { return "mode"; }
    constexpr static char const * strainFunctionString() { return "strainFunction"; }
    constexpr static char const * stressFunctionString() { return "stressFunction"; }
    constexpr static char const * numStepsString() { return "steps"; }
    constexpr static char const * outputString() { return "output"; }
    constexpr static char const * tangentString() { return "useNumericalTangent"; }
  };

  string m_solidMaterialName; ///< Material identifier
  string m_mode; ///< Test mode: triaxial, volumetric, oedometer
  string m_strainFunctionName; ///< Time-dependent function controlling strain (role depends on test mode)
  string m_stressFunctionName; ///< Time-dependent function controlling stress (role depends on test mode)
  int m_numSteps;    ///< Number of load steps
  string m_outputFileName; ///< Output file name
  int m_useNumericalTangent; ///< Flag to avoid using analytical tangent in driver

  array1d< real64 > m_time;
  array1d< real64 > m_axialStrain;
  array1d< real64 > m_axialStress;
  array1d< real64 > m_radialStrain;
  array1d< real64 > m_radialStress;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_TASKS_TRIAXIALDRIVER_HPP_ */
