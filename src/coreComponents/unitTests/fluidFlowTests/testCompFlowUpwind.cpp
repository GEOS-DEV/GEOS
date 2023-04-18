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
 * @file testCompFlowUpwind.cpp
 */
// Source includes
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVMUpwindUtilities.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "unitTests/fluidFlowTests/testFlowKernelHelpers.hpp"
// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::CompositionalMultiphaseFVMUpwindUtilities;
using UHelpers = geosx::CompositionalMultiphaseFVMUpwindUtilities::UpwindHelpers;

using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const *xmlInput =
        "<Problem>\n"
        "  <Solvers gravityVector=\"{-9.81, 0.0, 0.0}\">\n"
        "    <CompositionalMultiphaseFVM name=\"compflow\"\n"
        "                                 logLevel=\"0\"\n"
        "                                 discretization=\"fluidTPFA\"\n"
        "                                 targetRegions=\"{Region2}\"\n"
        "                                 temperature=\"297.15\"\n"
        "                                 useMass=\"1\">\n"
        "                                 \n"
        "      <NonlinearSolverParameters newtonTol=\"1.0e-6\"\n"
        "                                 newtonMaxIter=\"2\"/>\n"
        "      <LinearSolverParameters solverType=\"gmres\"\n"
        "                              krylovTol=\"1.0e-10\"/>\n"
        "    </CompositionalMultiphaseFVM>\n"
        "  </Solvers>\n"
        "  <Mesh>\n"
        "    <InternalMesh name=\"mesh1\"\n"
        "                  elementTypes=\"{C3D8}\" \n"
        "                  xCoords=\"{0, 2}\"\n"
        "                  yCoords=\"{0, 1}\"\n"
        "                  zCoords=\"{0, 1}\"\n"
        "                  nx=\"{2}\"\n"
        "                  ny=\"{1}\"\n"
        "                  nz=\"{1}\"\n"
        "                  cellBlockNames=\"{cb1}\"/>\n"
        "  </Mesh>\n"
        "  <NumericalMethods>\n"
        "    <FiniteVolume>\n"
        "      <TwoPointFluxApproximation name=\"fluidTPFA\""
        "         upwindSchemeName=\"hybrid\"/>\n"
        "    </FiniteVolume>\n"
        "  </NumericalMethods>\n"
        "  <ElementRegions>\n"
        "    <CellElementRegion name=\"Region2\" cellBlocks=\"{cb1}\" materialList=\"{fluid1, rock, relperm, cappressure}\" />\n"
        "  </ElementRegions>\n"
        "  <Constitutive>\n"
        "    <CompositionalMultiphaseFluid name=\"fluid1\"\n"
        "                                  phaseNames=\"{gas, oil}\"\n"
        "                                  equationsOfState=\"{PR, PR}\"\n"
        "                                  componentNames=\"{N2, C10, C20, H2O}\"\n"
        "                                  componentCriticalPressure=\"{34e5, 25.3e5, 14.6e5, 220.5e5}\"\n"
        "                                  componentCriticalTemperature=\"{126.2, 622.0, 782.0, 647.0}\"\n"
        "                                  componentAcentricFactor=\"{0.04, 0.443, 0.816, 0.344}\"\n"
        "                                  componentMolarWeight=\"{28e-3, 134e-3, 275e-3, 18e-3}\"\n"
        "                                  componentVolumeShift=\"{0, 0, 0, 0}\"\n"
        "                                  componentBinaryCoeff=\"{ {0, 0, 0, 0},\n"
        "                                                          {0, 0, 0, 0},\n"
        "                                                          {0, 0, 0, 0},\n"
        "                                                          {0, 0, 0, 0} }\"/>\n"
        "    <CompressibleSolidConstantPermeability name=\"rock\"\n"
        "                                           solidModelName=\"nullSolid\"\n"
        "                                           porosityModelName=\"rockPorosity\"\n"
        "                                           permeabilityModelName=\"rockPerm\"/>\n"
        "    <NullModel name=\"nullSolid\"/> \n"
        "   <PressurePorosity name=\"rockPorosity\"\n"
        "                     defaultReferencePorosity=\"0.05\"\n"
        "                     referencePressure = \"0.0\"\n"
        "                     compressibility=\"1.0e-9\"/>\n"
        "    <BrooksCoreyRelativePermeability name=\"relperm\"\n"
        "                                     phaseNames=\"{gas, oil}\"\n"
        "                                     phaseMinVolumeFraction=\"{0.15, 0.1}\"\n"
        "                                     phaseRelPermExponent=\"{2.0, 2.0}\"\n"
        "                                     phaseRelPermMaxValue=\"{0.9, 0.9}\"/>\n"
        "    <BrooksCoreyCapillaryPressure name=\"cappressure\"\n"
        "                                  phaseNames=\"{gas, oil}\"\n"
        "                                  phaseMinVolumeFraction=\"{0.05, 0.2}\"\n"
        "                                  phaseCapPressureExponentInv=\"{3.5,4.25}\"\n"
        "                                  phaseEntryPressure=\"{1e8, 0.}\"\n"
        "                                  capPressureEpsilon=\"0.0\"/> \n"
        "    <ConstantPermeability name=\"rockPerm\""
        "                          permeabilityComponents=\"{ 1.0e-17, 1.0e-17, 1.0e-17 }\" />\n"
        "  </Constitutive>\n"
        "  <FieldSpecifications>\n"
        "    <FieldSpecification name=\"initialPressure\"\n"
        "               initialCondition=\"1\"\n"
        "               setNames=\"{all}\"\n"
        "               objectPath=\"ElementRegions/Region2/cb1\"\n"
        "               fieldName=\"pressure\"\n"
        "               functionName=\"initialPressureFunc\"\n"
        "               scale=\"5e6\"/>\n"
        "    <FieldSpecification name=\"initialComposition_C10\"\n"
        "               initialCondition=\"1\"\n"
        "               setNames=\"{all}\"\n"
        "               objectPath=\"ElementRegions/Region2/cb1\"\n"
        "               fieldName=\"globalCompFraction\"\n"
        "               component=\"0\"\n"
        "               scale=\"0.099\"/>\n"
        "    <FieldSpecification name=\"initialComposition_N2\"\n"
        "               initialCondition=\"1\"\n"
        "               setNames=\"{all}\"\n"
        "               objectPath=\"ElementRegions/Region2/cb1\"\n"
        "               fieldName=\"globalCompFraction\"\n"
        "               component=\"1\"\n"
        "               functionName=\"initialN2\"\n"
        "               scale=\"1\"/>\n"
        "    <FieldSpecification name=\"initialComposition_C20\"\n"
        "               initialCondition=\"1\"\n"
        "               setNames=\"{all}\"\n"
        "               objectPath=\"ElementRegions/Region2/cb1\"\n"
        "               fieldName=\"globalCompFraction\"\n"
        "               component=\"2\"\n"
        "               functionName=\"initialC20\"\n"
        "               scale=\"1\"/>\n"
        "    <FieldSpecification name=\"initialComposition_H20\"\n"
        "               initialCondition=\"1\"\n"
        "               setNames=\"{all}\"\n"
        "               objectPath=\"ElementRegions/Region2/cb1\"\n"
        "               fieldName=\"globalCompFraction\"\n"
        "               component=\"3\"\n"
        "               scale=\"0.001\"/>\n"
        "  </FieldSpecifications>\n"
        "  <Functions>\n"
        "    <TableFunction name=\"initialPressureFunc\"\n"
        "                   inputVarNames=\"{elementCenter}\"\n"
        "                   coordinates=\"{0.0, 2.0}\"\n"
        "                   values=\"{.501, .5}\"/>\n"
        "    <TableFunction name=\"initialN2\"\n"
        "                   inputVarNames=\"{elementCenter}\"\n"
        "                   coordinates=\"{0.0, 2.0}\"\n"
        "                   values=\"{0.6, 0.3}\"/>\n"
        "    <TableFunction name=\"initialC20\"\n"
        "                   inputVarNames=\"{elementCenter}\"\n"
        "                   coordinates=\"{0.0, 2.0}\"\n"
        "                   values=\"{0.3, 0.6}\"/>\n"
        "  </Functions>"
        "</Problem>";
// Useful macro
#ifndef DEFINE_CAPB_FIELDS
#define DEFINE_CAP_FIELDS(active) \
  int const capPressureFlag = active; \
  auto phaseCapPressureView = getElementAccessor< true, 3, real64 >( &cap, \
                                                                     fields::cappres::phaseCapPressure::key(), \
                                                                     numFluxSupportPoints, \
                                                                     seri[iconn], \
                                                                     sesri[iconn], \
                                                                     sei[iconn], \
                                                                     1, \
                                                                     NP ); \
  auto dPhaseCapPressure_dPhaseVolFracView = getElementAccessor< true, 4, real64 >( &cap,                     \
                                                                                    fields::cappres::dPhaseCapPressure_dPhaseVolFraction::key(), \
                                                                                    numFluxSupportPoints, \
                                                                                    seri[iconn], \
                                                                                    sesri[iconn], \
                                                                                    sei[iconn], \
                                                                                    1, \
                                                                                    NP, \
                                                                                    NP );
#endif


#ifndef DEFINE_FLUID_FIELDS
#define DEFINE_FLUID_FIELDS() \
  auto phaseMassDensView = getElementAccessor< true, 3, real64 >( &fluid, \
                                                                  fields::multifluid::phaseMassDensity::key(), \
                                                                  numFluxSupportPoints, \
                                                                  seri[iconn], \
                                                                  sesri[iconn], \
                                                                  sei[iconn], 1, \
                                                                  NP ); \
  auto dPhaseMassDensView = getElementAccessor< true, 4, real64 >( &fluid, \
                                                                   fields::multifluid::dPhaseMassDensity::key(), \
                                                                   numFluxSupportPoints, \
                                                                   seri[iconn], \
                                                                   sesri[iconn], \
                                                                   sei[iconn], \
                                                                   1, \
                                                                   NP, NC+2 );
#endif

#ifndef DEFINE_SUBR_FIELDS
#define DEFINE_SUBR_FIELDS(A) \
  auto presView = getElementAccessor< true, 1, real64 >( &subRegion, \
                                                         fields::flow::pressure::key(), \
                                                         numFluxSupportPoints, \
                                                         seri[iconn], \
                                                         sesri[iconn], \
                                                         sei[iconn] ); \
  Array< Array< Array< double, 1 >, 1 >, 1 > gravCoefView; \
  if( A ) \
  { \
    gravCoefView = getElementAccessor< true, 1, real64 >( &subRegion, \
                                                          fields::flow::gravityCoefficient::key(), \
                                                          numFluxSupportPoints, \
                                                          seri[iconn], \
                                                          sesri[iconn], \
                                                          sei[iconn] ); \
  } \
  else \
  { \
    real64 const temp[1] = {0.0}; \
    gravCoefView = AccessorHelper< true >::template makeElementAccessor< 1 >( &(temp[0]), \
                                                                              numFluxSupportPoints, \
                                                                              seri[iconn], \
                                                                              sesri[iconn], \
                                                                              sei[iconn] ); \
  } \
  auto phaseMobView = getElementAccessor< true, 2, real64 >( &subRegion, \
                                                             fields::flow::phaseMobility::key(), \
                                                             numFluxSupportPoints, \
                                                             seri[iconn], \
                                                             sesri[iconn], \
                                                             sei[iconn], NP ); \
  auto dPhaseMobView = getElementAccessor< true, 3, real64 >( &subRegion, \
                                                              fields::flow::dPhaseMobility::key(), \
                                                              numFluxSupportPoints, \
                                                              seri[iconn], \
                                                              sesri[iconn], \
                                                              sei[iconn], \
                                                              NP, NC+2 ); \
  auto \
    dCompFrac_dCompDensView = getElementAccessor< true, 3, real64 >( &subRegion, \
                                                                     fields::flow::dGlobalCompFraction_dGlobalCompDensity::key(), \
                                                                     numFluxSupportPoints, \
                                                                     seri[iconn], \
                                                                     sesri[iconn], \
                                                                     sei[iconn], \
                                                                     NC, \
                                                                     NC ); \
  auto dPhaseVolFracView = getElementAccessor< true, 3, real64 >( &subRegion, \
                                                                  fields::flow::dPhaseVolumeFraction::key(), \
                                                                  numFluxSupportPoints, \
                                                                  seri[iconn], \
                                                                  sesri[iconn], \
                                                                  sei[iconn],        \
                                                                  NP, NC+2 );
#endif

//// Sphinx end before input XML
template<bool FULL, localIndex N, typename T, typename ... DIMS>
auto getElementAccessor(Group const *const group,
                        std::string const &key,
                        localIndex const stencilSize,
                        arraySlice1d<localIndex const> const &stencilRegIndices,
                        arraySlice1d<localIndex const> const &stencilSubRegIndices,
                        arraySlice1d<localIndex const> const &stencilElemIndices,
                        DIMS... otherDims) {

    ArrayView<T const, N> data;
    if (ElementSubRegionBase const *const subRegion = dynamicCast<ElementSubRegionBase const *const>(group)) {
        data = subRegion->getReference<Array<T, N> >(key);
    } else if (MultiFluidBase const *const fluid = dynamicCast<MultiFluidBase const *const>(group)) {
        data = fluid->getReference<Array<T, N> >(key);
    } else if (CapillaryPressureBase const *const cap = dynamicCast<CapillaryPressureBase const *const>(group)) {
        data = cap->getReference<Array<T, N> >(key);
    }

    data.move(LvArray::MemorySpace::host, false);
    auto view = AccessorHelper<FULL>::template makeElementAccessor<N, T>(data.data() + stencilElemIndices[0],
                                                                         stencilSize,
                                                                         stencilRegIndices,
                                                                         stencilSubRegIndices,
                                                                         stencilElemIndices,
                                                                         otherDims ...);
    return view;

}

using PermeabilityAccessors =
        StencilMaterialAccessors<PermeabilityBase,
                fields::permeability::permeability,
                fields::permeability::dPerm_dPressure>;

template<localIndex NC, localIndex NP>
void testCompositionalStandardUpwind(CompositionalMultiphaseFVM &solver,
                                     DomainPartition &domain) {
    solver.forDiscretizationOnMeshTargets(domain.getMeshBodies(), [&](string const,
                                                                      MeshLevel &mesh,
                                                                      arrayView1d<string const> const &regionNames) {
        ElementRegionManager &elementRegionManager = mesh.getElemManager();

        elementRegionManager.forElementSubRegions(regionNames,
                                                  [&](localIndex const,
                                                      ElementSubRegionBase &subRegion) {
                                                      SCOPED_TRACE(subRegion.getParent().getParent().getName() + "/" + subRegion.getName());

                                                      string const &fluidName = subRegion.getReference<string>(CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString());
                                                      string const &capName = subRegion.getReference<string>(
                                                              CompositionalMultiphaseBase::viewKeyStruct::capPressureNamesString());

                                                      PermeabilityAccessors permeabilityAccessors(elementRegionManager,
                                                                                                  solver.getName());
                                                      auto const &permeability = permeabilityAccessors.get(fields::permeability::permeability{});
                                                      auto const &dPerm_dPres = permeabilityAccessors.get(fields::permeability::dPerm_dPressure{});

                                                      MultiFluidBase &fluid = subRegion.getConstitutiveModel<MultiFluidBase>(fluidName);
                                                      CapillaryPressureBase &cap = subRegion.getConstitutiveModel<CapillaryPressureBase>(capName);

                                                      // reset the solver state to zero out variable updates
                                                      solver.resetStateToBeginningOfStep(domain);

                                                      FluxApproximationBase const &fluxApprox =
                                                              domain.getNumericalMethodManager().getFiniteVolumeManager().getFluxApproximation(solver.getDiscretizationName());
                                                      fluxApprox.forStencils<CellElementStencilTPFA>(mesh,
                                                                                                     [&](auto const &stencil) {

                                                                                                         auto const &seri = stencil.getElementRegionIndices();
                                                                                                         auto const &sesri = stencil.getElementSubRegionIndices();
                                                                                                         auto const &sei = stencil.getElementIndices();
                                                                                                         localIndex constexpr numFluxSupportPoints = 2;
                                                                                                         localIndex constexpr maxNumConn = 1;
                                                                                                         real64 trans[maxNumConn][numFluxSupportPoints]{};
                                                                                                         real64 dTrans_dP[maxNumConn][numFluxSupportPoints]{};

                                                                                                         for (localIndex iconn = 0; iconn < stencil.size(); ++iconn) {

                                                                                                             auto wrapper = stencil.createKernelWrapper();
                                                                                                             wrapper.computeWeights(
                                                                                                                     iconn,
                                                                                                                     permeability,
                                                                                                                     dPerm_dPres,
                                                                                                                     trans,
                                                                                                                     dTrans_dP);

                                                                                                             DEFINE_SUBR_FIELDS(1)
                                                                                                             //include cap pressure
                                                                                                             DEFINE_CAP_FIELDS(0);
                                                                                                             //fluid fetched fields
                                                                                                             DEFINE_FLUID_FIELDS()

                                                                                                             real64 potGrad = 0.;

                                                                                                             localIndex k_up[NP]{};
                                                                                                             real64 phaseFlux[NP]{0.0, 0.0};
                                                                                                             real64 dPhaseFlux_dP[NP][numFluxSupportPoints]{};
                                                                                                             real64 dPhaseFlux_dC[NP][numFluxSupportPoints][NC]{};

                                                                                                             localIndex k_up_expected[NP]{0,1};
                                                                                                             //classical mass balance equation way of computing potential and finding upwinding direction
                                                                                                             for (localIndex ip = 0; ip < NP; ++ip) {

                                                                                                                 UHelpers::computePPUPhaseFlux<NC, numFluxSupportPoints>(
                                                                                                                         NP,
                                                                                                                         ip,
                                                                                                                         {seri[iconn][0], seri[iconn][1]},
                                                                                                                         {sesri[iconn][0], sesri[iconn][1]},
                                                                                                                         {sei[iconn][0], sei[iconn][1]},
                                                                                                                         {trans[iconn][0], trans[iconn][1]},
                                                                                                                         {dTrans_dP[iconn][0], dTrans_dP[iconn][1]},
                                                                                                                         presView.toNestedViewConst(),
                                                                                                                         gravCoefView.toNestedViewConst(),
                                                                                                                         phaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         dCompFrac_dCompDensView.toNestedViewConst(),
                                                                                                                         phaseMassDensView.toNestedViewConst(),
                                                                                                                         dPhaseMassDensView.toNestedViewConst(),
                                                                                                                         phaseCapPressureView.toNestedViewConst(),
                                                                                                                         dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         capPressureFlag,
                                                                                                                         k_up[ip],
                                                                                                                         potGrad,
                                                                                                                         phaseFlux[ip],
                                                                                                                         dPhaseFlux_dP[ip],
                                                                                                                         dPhaseFlux_dC[ip]);

                                                                                                                 GEOSX_LOG_RANK(GEOSX_FMT( "Phase index :{}, Phase Flux :{}", ip, phaseFlux[ip]));
                                                                                                                 EXPECT_EQ(k_up_expected[ip],k_up[ip]);
                                                                                                             }
                                                                                                         }
                                                                                                     });
                                                  });
    });

}//EOfunc

template<localIndex NC, localIndex NP>
void testCompositionalUpwindHU(CompositionalMultiphaseFVM &solver,
                               DomainPartition &domain) {
    solver.forDiscretizationOnMeshTargets(domain.getMeshBodies(), [&](string const,
                                                                      MeshLevel &mesh,
                                                                      arrayView1d<string const> const &regionNames) {
        ElementRegionManager &elementRegionManager = mesh.getElemManager();

        elementRegionManager.forElementSubRegions(regionNames,
                                                  [&](localIndex const,
                                                      ElementSubRegionBase &subRegion) {
                                                      SCOPED_TRACE(subRegion.getParent().getParent().getName() + "/" + subRegion.getName());

                                                      string const &fluidName = subRegion.getReference<string>(CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString());
                                                      string const &capName = subRegion.getReference<string>(CompositionalMultiphaseBase::viewKeyStruct::capPressureNamesString());

                                                      PermeabilityAccessors permeabilityAccessors(elementRegionManager, solver.getName());
                                                      auto const &permeability = permeabilityAccessors.get(fields::permeability::permeability{});
                                                      auto const &dPerm_dPres = permeabilityAccessors.get(fields::permeability::dPerm_dPressure{});

                                                      MultiFluidBase &fluid = subRegion.getConstitutiveModel<MultiFluidBase>(fluidName);
                                                      CapillaryPressureBase &cap = subRegion.getConstitutiveModel<CapillaryPressureBase>(capName);

                                                      // reset the solver state to zero out variable updates
                                                      solver.resetStateToBeginningOfStep(domain);

                                                      FluxApproximationBase const &fluxApprox = domain.getNumericalMethodManager().getFiniteVolumeManager().getFluxApproximation(
                                                              solver.getDiscretizationName());
                                                      fluxApprox.forStencils<CellElementStencilTPFA>(mesh,
                                                                                                     [&](auto const &stencil) {

                                                                                                         auto const &seri = stencil.getElementRegionIndices();
                                                                                                         auto const &sesri = stencil.getElementSubRegionIndices();
                                                                                                         auto const &sei = stencil.getElementIndices();
                                                                                                         localIndex constexpr numFluxSupportPoints = 2;
                                                                                                         localIndex constexpr maxNumConn = 1;

                                                                                                         real64 trans[maxNumConn][numFluxSupportPoints]{};
                                                                                                         real64 dTrans_dP[maxNumConn][numFluxSupportPoints]{};

                                                                                                         real64 totFlux = 0.;
                                                                                                         real64 potGrad = 0.;

                                                                                                         localIndex k_up_hu_mu[NP]{}, k_up_hu_g[NP]{}, k_up_ppu[NP]{};
                                                                                                         real64 phaseFlux[NP]{0.0, 0.0};
                                                                                                         real64 dPhaseFlux_dP[NP][numFluxSupportPoints]{};
                                                                                                         real64 dPhaseFlux_dC[NP][numFluxSupportPoints][NC]{};

                                                                                                         //unused for now
                                                                                                         real64 fflow{};
                                                                                                         real64 dFflow_dP[numFluxSupportPoints]{};
                                                                                                         real64 dFflow_dC[numFluxSupportPoints][NC]{};
                                                                                                         real64 phaseFluxHU{};
                                                                                                         real64 dPhaseFluxHU_dP[numFluxSupportPoints]{};
                                                                                                         real64 dPhaseFluxHU_dC[numFluxSupportPoints][NC]{};

                                                                                                         for (localIndex iconn = 0; iconn < stencil.size(); ++iconn) {

                                                                                                             auto wrapper = stencil.createKernelWrapper();
                                                                                                             wrapper.computeWeights(
                                                                                                                     iconn,
                                                                                                                     permeability,
                                                                                                                     dPerm_dPres,
                                                                                                                     trans,
                                                                                                                     dTrans_dP);

                                                                                                             DEFINE_SUBR_FIELDS(1)
                                                                                                             //include cap pressure
                                                                                                             DEFINE_CAP_FIELDS(0);
                                                                                                             //fluid fetched fields
                                                                                                             DEFINE_FLUID_FIELDS()

                                                                                                             //classical mass balance equation way of computing potential and finding upwinding direction
                                                                                                             for (localIndex ip = 0; ip < NP; ++ip) {

                                                                                                                 UHelpers::computePPUPhaseFlux<NC, numFluxSupportPoints>(
                                                                                                                         NP,
                                                                                                                         ip,
                                                                                                                         {seri[iconn][0], seri[iconn][1]},
                                                                                                                         {sesri[iconn][0], sesri[iconn][1]},
                                                                                                                         {sei[iconn][0], sei[iconn][1]},
                                                                                                                         {trans[iconn][0], trans[iconn][1]},
                                                                                                                         {dTrans_dP[iconn][0], dTrans_dP[iconn][1]},
                                                                                                                         presView.toNestedViewConst(),
                                                                                                                         gravCoefView.toNestedViewConst(),
                                                                                                                         phaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         dCompFrac_dCompDensView.toNestedViewConst(),
                                                                                                                         phaseMassDensView.toNestedViewConst(),
                                                                                                                         dPhaseMassDensView.toNestedViewConst(),
                                                                                                                         phaseCapPressureView.toNestedViewConst(),
                                                                                                                         dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         capPressureFlag,
                                                                                                                         k_up_ppu[ip],
                                                                                                                         potGrad,
                                                                                                                         phaseFlux[ip],
                                                                                                                         dPhaseFlux_dP[ip],
                                                                                                                         dPhaseFlux_dC[ip]);

                                                                                                                 GEOSX_LOG_RANK( GEOSX_FMT( "Phase index :{}, Phase Flux :{}", ip, phaseFlux[ip]));

                                                                                                                 totFlux += phaseFlux[ip];
                                                                                                             }

                                                                                                             localIndex k_exp_hu_mu[NP] = {1, 1};
                                                                                                             localIndex k_exp_hu_g[NP] = {0, 1};

                                                                                                             for (localIndex ip = 0; ip < NP; ++ip) {

                                                                                                                 UpwindHelpers::computeFractionalFlow<NC, numFluxSupportPoints,
                                                                                                                         DrivingForces::Viscous,
                                                                                                                         HybridUpwind>(
                                                                                                                         NP,
                                                                                                                         ip,
                                                                                                                         {seri[iconn][0], seri[iconn][1]},
                                                                                                                         {sesri[iconn][0], sesri[iconn][1]},
                                                                                                                         {sei[iconn][0], sei[iconn][1]},
                                                                                                                         {trans[iconn][0], trans[iconn][1]},
                                                                                                                         {dTrans_dP[iconn][0], dTrans_dP[iconn][1]},
                                                                                                                         totFlux,
                                                                                                                         presView.toNestedViewConst(),
                                                                                                                         gravCoefView.toNestedViewConst(),
                                                                                                                         dCompFrac_dCompDensView.toNestedViewConst(),
                                                                                                                         phaseMassDensView.toNestedViewConst(),
                                                                                                                         dPhaseMassDensView.toNestedViewConst(),
                                                                                                                         phaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         phaseCapPressureView.toNestedViewConst(),
                                                                                                                         dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         capPressureFlag,
                                                                                                                         k_up_hu_mu[ip],
                                                                                                                         fflow,
                                                                                                                         dFflow_dP,
                                                                                                                         dFflow_dC);

                                                                                                                 phaseFluxHU += fflow * totFlux;

                                                                                                                 EXPECT_EQ(k_up_hu_mu[ip], k_exp_hu_mu[ip]);

                                                                                                                 GEOSX_LOG_RANK(GEOSX_FMT("1/ ip {} , phaseFlux(PPU) {}, phaseFlux(HU) {}",
                                                                                                                                          ip, phaseFlux[ip], phaseFluxHU ));

                                                                                                                 localIndex k_up_hu_og = -1;
                                                                                                                 UpwindHelpers::computePotentialFluxes<NC,
                                                                                                                         DrivingForces::Gravity,
                                                                                                                         numFluxSupportPoints, HybridUpwind>(
                                                                                                                         NP,
                                                                                                                         ip,
                                                                                                                         {seri[iconn][0], seri[iconn][1]},
                                                                                                                         {sesri[iconn][0], sesri[iconn][1]},
                                                                                                                         {sei[iconn][0], sei[iconn][1]},
                                                                                                                         {trans[iconn][0], trans[iconn][1]},
                                                                                                                         {dTrans_dP[iconn][0], dTrans_dP[iconn][1]},
                                                                                                                         totFlux,
                                                                                                                         presView.toNestedViewConst(),
                                                                                                                         gravCoefView.toNestedViewConst(),
                                                                                                                         phaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseMobView.toNestedViewConst(),
                                                                                                                         dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         dCompFrac_dCompDensView.toNestedViewConst(),
                                                                                                                         phaseMassDensView.toNestedViewConst(),
                                                                                                                         dPhaseMassDensView.toNestedViewConst(),
                                                                                                                         phaseCapPressureView.toNestedViewConst(),
                                                                                                                         dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                                                                         capPressureFlag,
                                                                                                                         k_up_hu_g[ip],
                                                                                                                         k_up_hu_og,
                                                                                                                         phaseFluxHU,
                                                                                                                         dPhaseFluxHU_dP,
                                                                                                                         dPhaseFluxHU_dC
                                                                                                                 );

                                                                                                                 GEOSX_LOG_RANK(GEOSX_FMT("2/ ip {} , phaseFlux(PPU) {}, phaseFlux(HU) {}",
                                                                                                                                          ip, phaseFlux[ip], phaseFluxHU ));
                                                                                                                 EXPECT_EQ( k_up_hu_g[ip], k_exp_hu_g[ip]);
                                                                                                             }
                                                                                                         }
                                                                                                     });
                                                  });
    });

}


class CompositionalMultiphaseFlowUpwindHelperKernelsTest : public ::testing::Test {
public:

    CompositionalMultiphaseFlowUpwindHelperKernelsTest()
            :
            state(std::make_unique<CommandLineOptions>(g_commandLineOptions)) {}

protected:

    void SetUp() override {
        setupProblemFromXML(state.getProblemManager(), xmlInput);
        solver = &state.getProblemManager().getPhysicsSolverManager().getGroup<CompositionalMultiphaseFVM>(
                "compflow");

        DomainPartition &domain = state.getProblemManager().getDomainPartition();

        solver->setupSystem(domain,
                            solver->getDofManager(),
                            solver->getLocalMatrix(),
                            solver->getSystemRhs(),
                            solver->getSystemSolution());

        solver->implicitStepSetup(time, dt, domain);
    }

    static real64 constexpr time = 0.0;
    static real64 constexpr dt = 1e-2;
    static real64 constexpr eps = std::numeric_limits<real64>::epsilon();

    GeosxState state;
    CompositionalMultiphaseFVM *solver;
};

real64 constexpr CompositionalMultiphaseFlowUpwindHelperKernelsTest::time;
real64 constexpr CompositionalMultiphaseFlowUpwindHelperKernelsTest::dt;
real64 constexpr CompositionalMultiphaseFlowUpwindHelperKernelsTest::eps;

TEST_F(CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_standardPPU) {
    localIndex constexpr NP = 2;
    localIndex constexpr NC = 4;
    DomainPartition &domain = state.getProblemManager().getDomainPartition();
    testCompositionalStandardUpwind<NC, NP>(*solver, domain);
}

TEST_F(CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_HU) {
    localIndex constexpr NP = 2;
    localIndex constexpr NC = 4;
    DomainPartition &domain = state.getProblemManager().getDomainPartition();
    testCompositionalUpwindHU<NC, NP>(*solver, domain);
}


int main(int argc,
         char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    g_commandLineOptions = *geosx::basicSetup(argc, argv);
    int const result = RUN_ALL_TESTS();
    geosx::basicCleanup();
    return result;
}
