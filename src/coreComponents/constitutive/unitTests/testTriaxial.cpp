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

#include "gtest/gtest.h"

#include "managers/ProblemManager.hpp"
#include "managers/initialization.hpp"
#include "managers/GeosxState.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/Functions/TableFunction.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/CamClay.hpp"
#include "constitutive/solid/InvariantDecompositions.hpp"

#include "dataRepository/xmlWrapper.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

struct StrainData
{
  real64 strainIncrement[6] = { };
};


template< typename CMW >
void getStress( CMW const cmw,
                real64 (& stress)[6] )
{
  forAll< serialPolicy >( 1, [&stress, cmw] ( localIndex const k )
  {
    stress[0] = cmw.m_newStress( k, 0, 0 );
    stress[1] = cmw.m_newStress( k, 0, 1 );
    stress[2] = cmw.m_newStress( k, 0, 2 );
    stress[3] = cmw.m_newStress( k, 0, 3 );
    stress[4] = cmw.m_newStress( k, 0, 4 );
    stress[5] = cmw.m_newStress( k, 0, 5 );
  } );
}


template< typename POLICY >
void triaxialDriver()
{
  getGlobalState().getProblemManager().parseCommandLineInput();
  getGlobalState().getProblemManager().parseInputFile();

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  FunctionBase & function = functionManager.getGroup< FunctionBase >( "timeFunction" );

  TableFunction & table = dynamicCast< TableFunction & >( function );

  table.initializeFunction();

  real64 time = 0.01;
  real64 result = table.evaluate( &time );

  std::cout << result << std::endl;
/*
   // create a Cam-Clay model, and test xml input
   conduit::Node node;
   dataRepository::Group rootGroup( "root", node );
   ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

   //real64 const friction = 30.0; // will use later in checks

   string const inputStream =
    "<Constitutive>"
    "   <CamClay"
    "      name=\"granite\" "
    "      defaultDensity=\"2700\" "
    "      defaultRefPressure=\"-10.0\" "
    "      defaultRefStrainVol=\"-1e-4\" "
    "      defaultShearModulus=\"5400.0\" "
    "      defaultBulkModulus=\"5400.0\" "
    "      defaultPreConsolidationPressure =\"-90.0\" "
    "      defaultShapeParameter=\"1.0\" "
    "      defaultCslSlope=\"1.05\" "
    "      defaultVirginCompressionIndex=\"0.13\" "
    "      defaultRecompressionIndex=\"0.018\"/>"
    "</Constitutive>";

   xmlWrapper::xmlDocument xmlDocument;
   xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(),
                                                             inputStream.size() );
   if( !xmlResult )
   {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
   }

   xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
   constitutiveManager.processInputFileRecursive( xmlConstitutiveNode );
   constitutiveManager.postProcessInputRecursive();

   localIndex constexpr numElem = 2;
   localIndex constexpr numQuad = 4;

   dataRepository::Group disc( "discretization", &rootGroup );
   disc.resize( numElem );

   CamClay & cm = constitutiveManager.getConstitutiveRelation< CamClay >( "granite" );
   cm.allocateConstitutiveData( disc, numQuad );

   // confirm allocation sizes

   EXPECT_EQ( cm.size(), numElem );
   EXPECT_EQ( cm.numQuadraturePoints(), numQuad );

   // use updates class to run a uniaxial compression test.  the material
   // begins with some cohesion, but it quickly degrades to zero under
   // plastic loading.  at the end of the loading, we confirm the
   // stress point lies on the correctly-positioned yield surface.

   CamClay::KernelWrapper cmw = cm.createKernelUpdates();

   StrainData data;
   real64 inc = -1e-3;
   real64 total = 0;
   data.strainIncrement[0] = inc;
   data.strainIncrement[1] = 0;
   data.strainIncrement[2] = 0;
   data.strainIncrement[3] = 0;
   data.strainIncrement[4] = 0;
   data.strainIncrement[5] = 0;


   for( localIndex loadstep=0; loadstep < 50; ++loadstep )
   {
    forAll< parallelDevicePolicy<> >( 1, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 stress[6] = {0};
      real64 stiffness[6][6] = {{0}};
      cmw.smallStrainUpdate( k, 0, data.strainIncrement, stress, stiffness );
    } );
    cm.saveConvergedState();
    total = loadstep;
   }



   real64 stress[6] = {0};
   getStress( cmw, stress );

   // loading was set up to drive to total cohesion loss (c=0), at which
   // point the Q/P ratio should equal the slope of the DP yield surface:

   real64 invariantP, invariantQ;
   real64 deviator[6];

   twoInvariant::stressDecomposition( stress,
                                     invariantP,
                                     invariantQ,
                                     deviator );

   //real64 phi = friction * M_PI / 180;
   //real64 slope = -6 * sin( phi ) / ( 3 - sin( phi ) );
   //std::cout << invariantP << " " << invariantQ << " " << total << std::endl;



   //EXPECT_TRUE( fabs( invariantQ / invariantP / slope - 1 ) < 1e-8 );

   // we now use a finite-difference check of tangent stiffness to confirm
   // the analytical form is working properly.

   cmw.checkSmallStrainStiffness( 0, 0, data.strainIncrement, true );

   EXPECT_TRUE( cmw.checkSmallStrainStiffness( 0, 0, data.strainIncrement ) );
 */
}


#ifdef USE_CUDA
TEST( TriaxialTests, TriaxialDevice )
{
  triaxialDriver< geosx::parallelDevicePolicy< > >();
}
#endif
TEST( TriaxialTests, TriaxialHost )
{
  triaxialDriver< serialPolicy >();
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv, true ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
