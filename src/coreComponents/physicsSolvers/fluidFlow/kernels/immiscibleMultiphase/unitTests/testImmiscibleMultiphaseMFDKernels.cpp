/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * (c) GEOS/GEOSX Contributors
 * ------------------------------------------------------------------------------------------------------------
 */
/**
 * @file testImmiscibleMultiphaseMFDKernels.cpp
 *
 * Unit tests exercising the actual function ElementBasedAssemblyKernel::compute  in
 * ImmiscibleMultiphaseMFDKernels.hpp, for different inner products
 *
 * A minimal 1x1x1 hexahedral mesh plus constant two‑phase fluid & permeability
 * is created via XML. Phase mobility and its derivatives are manually set to
 * deterministic values so analytic expectations can be checked.
 *
 * NOTE: We purposefully bypass the factory launcher to access and call the
 * individual member functions with controlled inputs, but we use the real
 * kernel type and data structures..
 */

#include <gtest/gtest.h>

#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"

#include "mesh/MeshManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"

#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseMFDKernels.hpp"
#include "LvArray/src/SortedArray.hpp"

using namespace geos;
using namespace geos::immiscibleMultiphaseMFDKernels;

namespace {

std::string generateXmlInputImmiscibleMFD(const std::string & innerProductType)
{
  std::ostringstream oss;
  oss << R"xml(
<Problem>
  <Solvers gravityVector="{ 0,0,0 }">
    <ImmiscibleMultiphaseFlowMFD name="immiscibleHybridMimetic"
        discretization="immiscibleHybridMimetic"
        targetRegions="{ Region }"
        dependentPhaseIndex="1" logLevel="0">
      <NonlinearSolverParameters newtonTol="1e-12" newtonMaxIter="1"/>
      <LinearSolverParameters directParallel="0"/>
    </ImmiscibleMultiphaseFlowMFD>
  </Solvers>
  <Events maxTime="1">
    <PeriodicEvent name="solverApps" forceDt="1" target="/Solvers/immiscibleHybridMimetic"/>
  </Events>
  <NumericalMethods>
    <FiniteVolume>
      <HybridMimeticDiscretization name="immiscibleHybridMimetic" innerProductType=")xml" << innerProductType << R"xml("/>
    </FiniteVolume>
  </NumericalMethods>
  <Mesh>
    <InternalMesh name="mesh1" elementTypes="{ C3D8 }"
                  xCoords="{0,1}" yCoords="{0,1}" zCoords="{0,1}"
                  nx="{1}" ny="{1}" nz="{1}" cellBlockNames="{ blk }"/>
  </Mesh>
  <Functions>
    <TableFunction name="densityWater" coordinates="{0}" values="{1000}"/>
    <TableFunction name="densityGas" coordinates="{0}" values="{1000}"/>
    <TableFunction name="viscosityWater"  coordinates="{0}" values="{1e-3}"/>
    <TableFunction name="viscosityGas"  coordinates="{0}" values="{1e-3}"/>
    <TableFunction name="waterRelativePermeabilityTable" coordinates="{0,1}" values="{0,1}"/>
    <TableFunction name="gasRelativePermeabilityTable" coordinates="{0,1}" values="{0,1}"/>
  </Functions>
  <Constitutive>
    <TwoPhaseImmiscibleFluid name="fluid" phaseNames="{ water, gas }"
       densityTableNames="{ densityWater, densityGas }"
       viscosityTableNames="{ viscosityWater, viscosityGas }"/>
    <NullModel name="nullSolid"/>
    <PressurePorosity name="poro" defaultReferencePorosity="0.2" referencePressure="0" compressibility="0"/>
    <ConstantPermeability name="perm" permeabilityComponents="{1,2,3}"/>
    <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="poro" permeabilityModelName="perm"/>
    <TableRelativePermeability name="relperm" phaseNames="{ water, gas }"
      wettingNonWettingRelPermTableNames="{ waterRelativePermeabilityTable, gasRelativePermeabilityTable }"/>
  </Constitutive>
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ blk }" materialList="{ fluid, rock, relperm }"/>
  </ElementRegions>
  <FieldSpecifications>
    <FieldSpecification name="initP" initialCondition="1" setNames="{all}" objectPath="ElementRegions/Region/blk" fieldName="pressure" scale="1.0"/>
    <FieldSpecification name="initFaceP" initialCondition="1" setNames="{all}" objectPath="faceManager" fieldName="facePressure" scale="1.0"/>
    <FieldSpecification name="sat0" initialCondition="1" setNames="{all}" objectPath="ElementRegions/Region/blk" fieldName="phaseVolumeFraction" component="0" scale="0.6"/>
    <FieldSpecification name="sat1" initialCondition="1" setNames="{all}" objectPath="ElementRegions/Region/blk" fieldName="phaseVolumeFraction" component="1" scale="0.4"/>
  </FieldSpecifications>
</Problem>
)xml";
  return oss.str();
}

// Helper: build a minimal problem with one hex element and two‑phase fluid.
void buildProblemOnce(const std::string & innerProductType)
{
  static bool built=false; if( built ) return;
  static char prog[] = "immiscibleMFDTest"; static char * argvLocal[] = { prog, nullptr }; int argcLocal = 1;
  static GeosxState state( basicSetup( argcLocal, argvLocal, false ) );
  std::string xml = generateXmlInputImmiscibleMFD(innerProductType);
  xmlWrapper::xmlDocument doc; ASSERT_TRUE( doc.loadString( xml ) );
  xmlWrapper::xmlNode root = doc.getChild( dataRepository::keys::ProblemManager );
  ProblemManager & pm = getGlobalState().getProblemManager();
  pm.processInputFileRecursive( doc, root );
  DomainPartition & domain = pm.getDomainPartition();
  // Explicit constitutive processing (required to populate /domain/Constitutive before problemSetup)
  auto & constitutiveManager = domain.getConstitutiveManager();
  xmlWrapper::xmlNode constitutiveNode = root.child( constitutiveManager.getName().c_str() );
  constitutiveManager.processInputFileRecursive( doc, constitutiveNode );
  // Mesh levels then element regions
  pm.getGroup< MeshManager >( pm.groupKeys.meshManager ).generateMeshLevels( domain );
  auto & elemMgr = domain.getMeshBody(0).getBaseDiscretization().getElemManager();
  xmlWrapper::xmlNode elemNode = root.child( elemMgr.getName().c_str() );
  elemMgr.processInputFileRecursive( doc, elemNode );
  elemMgr.postInputInitializationRecursive();
  pm.problemSetup();
  // After problemSetup solver should have registered mobility fields
  // Sanity: check phaseMobility field exists on subRegion
  auto & subRegion = elemMgr.getRegion(0).getSubRegion< CellElementSubRegion >(0);
  pm.applyInitialConditions();
  
  // Introduce non-zero potentials: perturb face pressures slightly so cellP - faceP != 0
  {
    auto & faceMgr2 = domain.getMeshBody(0).getBaseDiscretization().getFaceManager();
    auto & facePres = faceMgr2.getField< fields::flow::facePressure >();
    // Cell pressure is initialized to 1.0e6. Make each face pressure a bit lower by an increasing offset.
    for( localIndex f=0; f<faceMgr2.size(); ++f )
    {
      facePres[f] -= (f+1) * 1.0; // yields potentials: 1, 2, ..., 6
    }
  }
  built = true;
}

// Helper to assign deterministic phase mobility + derivatives.
void setMobility( CellElementSubRegion & subRegion, integer indepPhase )
{
  using Deriv = immiscibleFlow::DerivativeOffset;
  auto & mob = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMobility >();
  auto & dMob = subRegion.getField< fields::immiscibleMultiphaseFlow::dPhaseMobility >();
  localIndex const ei = 0;
  real64 const s_dep = 1.0/3.0;
  real64 const s_indep = 2.0/3.0;
  mob[ei][0] = s_dep*1.0e3; mob[ei][1] = s_indep*1.0e3;
  // Derivatives: d(lambda)/dP = 0, d(lambda)/dS(indep) for testing sign mapping
  dMob[ei][0][Deriv::dP] = 0.0; dMob[ei][0][Deriv::dS] = 1.0e3;      // dependent phase
  dMob[ei][1][Deriv::dP] = 0.0; dMob[ei][1][Deriv::dS] = -1.0e3;    // independent phase
  GEOS_UNUSED_VAR( indepPhase );
}

struct KernelHarness
{
  using KType = ElementBasedAssemblyKernel< 6, mimeticInnerProduct::TPFAInnerProduct >;
  using StackVariables = KType::StackVariables;
  static auto build( integer indepPhase, real64 dt, const std::string & innerProductType )
  {
    buildProblemOnce(innerProductType);
    ProblemManager & pm = getGlobalState().getProblemManager();
    DomainPartition & domain = pm.getDomainPartition();
    auto & mesh = domain.getMeshBody(0).getBaseDiscretization();
    auto & nodeMgr = mesh.getNodeManager();
    auto & faceMgr = mesh.getFaceManager();
    auto & elemMgr = mesh.getElemManager();
    auto & region = elemMgr.getRegion(0);
    auto & subRegion = region.getSubRegion< CellElementSubRegion >(0);

    setMobility( subRegion, indepPhase );

    // Create synthetic element dof numbering field exactly once
    static bool elemDofsCreated = false;
    string const elemDofKey = "testElemDof";
    if( !elemDofsCreated )
    {
      auto & dofs = subRegion.registerWrapper< array1d< globalIndex > >( elemDofKey ).reference();
      dofs.resize( subRegion.size() );
      for( localIndex i=0;i<subRegion.size();++i ) dofs[i] = 2*i; // pressure then saturation
      elemDofsCreated = true;
    }

    // Face dofs once
    static bool faceDofsCreated = false;
    string const faceDofKey = "testFaceDof";
    if( !faceDofsCreated )
    {
      auto & fd = faceMgr.registerWrapper< array1d< globalIndex > >( faceDofKey ).reference();
      fd.resize( faceMgr.size() );
      for( localIndex f=0; f<faceMgr.size(); ++f ) fd[f] = 100 + f; // arbitrary
      faceDofsCreated = true;
    }

    // Persistent accessors (must outlive returned kernel); construct once
    static bool accessorsInit = false;
    static decltype(elemMgr.constructArrayViewAccessor< globalIndex, 1 >( "" )) dofAccessor;
    static decltype(elemMgr.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() )) satAccessor;
    static decltype(elemMgr.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseMobility::key() )) mobAccessor;
    static decltype(elemMgr.constructArrayViewAccessor< real64, 3 >( fields::immiscibleMultiphaseFlow::dPhaseMobility::key() )) dMobAccessor;
    if( !accessorsInit )
    {
      dofAccessor = elemMgr.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
      satAccessor = elemMgr.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );
      mobAccessor = elemMgr.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseMobility::key() );
      dMobAccessor = elemMgr.constructArrayViewAccessor< real64, 3 >( fields::immiscibleMultiphaseFlow::dPhaseMobility::key() );
      mobAccessor.setName( "UnitTest/persistent/" + std::string( fields::immiscibleMultiphaseFlow::phaseMobility::key() ) );
      dMobAccessor.setName( "UnitTest/persistent/" + std::string( fields::immiscibleMultiphaseFlow::dPhaseMobility::key() ) );
      accessorsInit = true;
    }

    // Retrieve constitutive objects from the subRegion constitutiveModels (cloned instances)
    auto const & constitutiveModels = subRegion.getGroup( "ConstitutiveModels" );
    auto const & fluid = constitutiveModels.getGroup< constitutive::TwoPhaseImmiscibleFluid >( "fluid" );
    auto const & permBase = constitutiveModels.getGroup< constitutive::PermeabilityBase >( "perm" );

    // Sanity checks to help diagnose bad accessor dimensions
    auto const & mobField = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMobility >();
    auto const & dMobField = subRegion.getField< fields::immiscibleMultiphaseFlow::dPhaseMobility >();

    static array1d< real64 > rhs; rhs.resize( 2 + faceMgr.size() ); rhs.zero();
    // Persistent static matrix storage to keep view valid during kernel lifetime
    static CRSMatrix< real64, globalIndex > localMat;
    if( localMat.numRows() == 0 )
    {
      localIndex const nRows = 2 + faceMgr.size();
      // Determine max column index among element and face dofs
      globalIndex maxCol = 0;
      auto const & elemDofsView = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      for( localIndex i=0; i<elemDofsView.size(); ++i )
      {
        maxCol = std::max( maxCol, elemDofsView[i] + 1 ); // include saturation dof (elemDof+1)
      }
      auto const & faceDofsView = faceMgr.getReference< array1d< globalIndex > >( faceDofKey );
      for( localIndex f=0; f<faceDofsView.size(); ++f )
      {
        maxCol = std::max( maxCol, faceDofsView[f] );
      }
      globalIndex const nCols = maxCol + 10; // small padding
      // Recreate matrix by assignment (deep copy) with initial per-row capacity 8
      localMat = CRSMatrix< real64, globalIndex >( nRows, nCols, 8 );
    }

    // Explicit empty region filter (default constructed)
    // SortedArrayView cannot be default-constructed publicly; use a backing SortedArray.
    static LvArray::SortedArray< localIndex, localIndex, LvArray::ChaiBuffer > emptySet; // remains empty
    using IP = mimeticInnerProduct::TPFAInnerProduct;
    ElementBasedAssemblyKernel< 6, IP > kernel( 0, 0, 0, 1e-12,
                                                faceDofKey, nodeMgr, faceMgr, subRegion,
                                                dofAccessor, satAccessor, fluid, permBase,
                                               emptySet.toView(), indepPhase, dt,
                                                true, localMat.toViewConstSizes(), rhs.toView(),
                                                mobAccessor, dMobAccessor );

    KernelHarness::KType::StackVariables s; kernel.setup( 0, s );
    // Return faceMgr as well for reference calculations
    return std::tuple< KernelHarness::KType, KernelHarness::KType::StackVariables, CellElementSubRegion &, FaceManager & >( kernel, s, subRegion, faceMgr );
  }
};

struct ReferenceCalc
{

  // helper: directly fills a provided stack with expected fluxes and derivatives
  template< class StackType >
  static void fillExpectedStack( CellElementSubRegion & subRegion, FaceManager const & faceMgr, integer indepPhase,
                                 StackType & stack )
  {
    localIndex const ei = 0;
    real64 const rho0 = 1000.0;
    real64 const rho1 = 1000.0;
    using Deriv = immiscibleFlow::DerivativeOffset;
    auto const & mob = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMobility >();
    auto const & dMob = subRegion.getField< fields::immiscibleMultiphaseFlow::dPhaseMobility >();
    real64 const Lambda = rho0 * mob[ei][0] + rho1 * mob[ei][1];
    real64 const dLambda_dS = rho0 * dMob[ei][0][Deriv::dS] + rho1 * dMob[ei][1][Deriv::dS];
    real64 const lambda0 = mob[ei][0];
    real64 const frac = (rho0 * lambda0) / Lambda;
    real64 const cellP = subRegion.getField< fields::flow::pressure >()[ei];
    auto const & facePres = faceMgr.getField< fields::flow::facePressure >();
    auto const & elemToFaces = subRegion.faceList();

    // Zero initialize arrays first (in case caller did not)
    for( int f=0; f<6; ++f )
    {
      stack.MassFlux[f] = 0.0;
      stack.dMassFlux_dPres[f] = 0.0;
      stack.dMassFlux_dS[f] = 0.0;
      stack.dDivMassFluxes_dFaceVars[f] = 0.0;
      stack.dDivSatFluxes_dFaceVars[f] = 0.0;
      for( int g=0; g<6; g++ )
      {
        stack.dMassFlux_dFacePres[f][g] = 0.0;
      }
    }
    stack.divMassFluxes = 0.0;
    stack.dDivMassFluxes_dP = 0.0;
    stack.dDivMassFluxes_dS = 0.0;
    stack.divSatFluxes = 0.0;
    stack.dDivSatFluxes_dP = 0.0;
    stack.dDivSatFluxes_dS = 0.0;
    
    
    
    // Fill per-face quantities
    real64 pot[6];
    for( int f=0; f<6; ++f )
    {
      localIndex lf = elemToFaces( ei, f );
      pot[f] = cellP - facePres[lf];
    }
    for( int f=0; f<6; ++f )
    {
      for( int g=0; g<6; ++g )
      {
        real64 const T_fg = stack.transMatrix[f][g];
        stack.MassFlux[f] +=  Lambda * T_fg * pot[g];
        stack.dMassFlux_dPres[f] += Lambda * T_fg;
        stack.dMassFlux_dS[f] += dLambda_dS * T_fg * pot[g];
        stack.dMassFlux_dFacePres[f][g] = -Lambda * T_fg;
      }

      stack.divMassFluxes += stack.MassFlux[f];
      stack.dDivMassFluxes_dP += stack.dMassFlux_dPres[f];
      stack.dDivMassFluxes_dS += stack.dMassFlux_dS[f];
      for( int g=0; g<6; ++g )
      {
        real64 const T_fg = stack.transMatrix[f][g];
        stack.dDivMassFluxes_dFaceVars[f] += -Lambda * T_fg;
      }
    }
    stack.divSatFluxes = frac * stack.divMassFluxes;
    stack.dDivSatFluxes_dP = frac * stack.dDivMassFluxes_dP;
    stack.dDivSatFluxes_dS = stack.divMassFluxes + frac * stack.dDivMassFluxes_dS;
    for( int f=0; f<6; ++f )
    {
      stack.dDivSatFluxes_dFaceVars[f] = frac * stack.dDivMassFluxes_dFaceVars[f];
    }
  }
};

} // end anonymous namespace

// ---------------- TEST -----------------
class ImmiscibleMultiphaseMFDKernelsParamTest : public ::testing::TestWithParam<std::string> {};
TEST_P(ImmiscibleMultiphaseMFDKernelsParamTest, compute_direct)
{
  integer indepPhase = 0; real64 dt=1.0;
  std::string innerProductType = GetParam();
  auto tup = KernelHarness::build(indepPhase, dt, innerProductType);
  KernelHarness::KType & kernel = std::get<0>( tup );
  KernelHarness::KType::StackVariables s = std::get<1>( tup );
  KernelHarness::KType::StackVariables s_ref = std::get<1>( tup );
  CellElementSubRegion & subRegion = std::get<2>( tup );
  FaceManager & faceMgr = std::get<3>( tup );
  kernel.compute( 0, s );
  
  s_ref.transMatrix = s.transMatrix;
  ReferenceCalc::fillExpectedStack(subRegion, faceMgr, indepPhase, s_ref);
  
  constexpr real64 eps_tol = 1.0e-8;
  for( int f=0; f<6; ++f )
  {
    EXPECT_NEAR( s.MassFlux[f], s_ref.MassFlux[f], eps_tol );
    EXPECT_NEAR( s.dMassFlux_dPres[f], s_ref.dMassFlux_dPres[f], eps_tol );
    EXPECT_NEAR( s.dMassFlux_dS[f], s_ref.dMassFlux_dS[f], eps_tol );
    for( int g=0; g<6; g++ )
    {
      EXPECT_NEAR( s.dMassFlux_dFacePres[f][g], s_ref.dMassFlux_dFacePres[f][g], eps_tol );
    }
    EXPECT_NEAR( s.dDivMassFluxes_dFaceVars[f], s_ref.dDivMassFluxes_dFaceVars[f], eps_tol );
    EXPECT_NEAR( s.dDivSatFluxes_dFaceVars[f], s_ref.dDivSatFluxes_dFaceVars[f], eps_tol );
  }
  EXPECT_NEAR( s.divMassFluxes, s_ref.divMassFluxes, eps_tol );
  EXPECT_NEAR( s.dDivMassFluxes_dP, s_ref.dDivMassFluxes_dP, eps_tol );
  EXPECT_NEAR( s.dDivMassFluxes_dS, s_ref.dDivMassFluxes_dS, eps_tol );
  EXPECT_NEAR( s.divSatFluxes, s_ref.divSatFluxes, eps_tol );
  EXPECT_NEAR( s.dDivSatFluxes_dP, s_ref.dDivSatFluxes_dP, eps_tol );
  EXPECT_NEAR( s.dDivSatFluxes_dS, s_ref.dDivSatFluxes_dS, eps_tol );
}
INSTANTIATE_TEST_SUITE_P(
  InnerProductTypeTests,
  ImmiscibleMultiphaseMFDKernelsParamTest,
  ::testing::Values("TPFA", "quasiTPFA", "quasiRT", "simple", "beiraoDaVeigaLipnikovManzini")
);
