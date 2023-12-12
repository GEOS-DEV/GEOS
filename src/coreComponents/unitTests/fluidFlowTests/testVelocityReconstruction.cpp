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

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/Logger.hpp"
#include "mainInterface/initialization.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "testFlowKernelHelpers.hpp"


// TPL includes
#include <gtest/gtest.h>

using namespace geos;

enum class CASE {
    ALIGNED,
    TILTED
};

void compareVelocities(real64 const (&computedVel)[3],
                       real64 const (&expectedVel)[3]) {
    for (int dir = 0; dir < 3; ++dir) {
        geos::testing::checkRelativeError(computedVel[dir], expectedVel[dir], 1e-15);
    }
}


void setReferences(CASE const &tested, real64 (&expectedVel)[3]) {
    if (tested == CASE::ALIGNED)//faces are aligned with x-y-z
    {
        expectedVel[0] = 1.0;
        expectedVel[1] = 1.0;
        expectedVel[2] = 1.0;
    }
//    else if(tested == CASE::TILTED)//case is tilted through normals
//    {
//        expectedVel[0] = xxx;
//        expectedVel[1] = xxx;
//        expectedVel[2] = xxx;
//    }
}


TEST(testAligned, Velocity_aligned) {
    constexpr int nfaces = 1;  //2-hexa problem
    constexpr int numElemtInStencil = 2;

    real64 const faceNormal[nfaces][3] = {{1., 0., 0.}};
    real64 const cellToFaceVec[nfaces][numElemtInStencil][3] = {{{-1., 0., 0.}, {1., 0., 0}}};
    real64 const transMult[nfaces] = {1.};
    real64 const geomStabSum = 1.;  //irrelevant for now

    localIndex elementRegionIndices[] = {0, 0};
    localIndex elementSubRegionIndices[] = {0, 0};
    localIndex elementIndices[] = {0, 1};
    real64 weight[] = {1, -1};

    localIndex kf = 0;   //aka face index
    CellElementStencilTPFA tpfa;
    tpfa.add(2, elementRegionIndices, elementSubRegionIndices, elementIndices, weight, kf);
    tpfa.addVectors(transMult[kf], geomStabSum, faceNormal[kf], cellToFaceVec[kf]);

    CellElementStencilTPFA::KernelWrapper wrapper = tpfa.createKernelWrapper();

    const real64 flux = 0.001;
    constexpr localIndex nPhases = 2;
    constexpr localIndex nDirs = 3;
    const real64 phaseVelocity[numElemtInStencil][nPhases][nDirs] = { {{0.,0.,0.}, {0.,0.,0.}}, {{0.,0.,0.}, {0.,0.,0.}}, };

    array2d<real64> arr(2,3);
//    std::ptrdiff_t const sizes[2] = {1, 3};
//    arr.resize(2, sizes);
    arr[0][0] = 2.;
    arr[1][0] = -1.;
    arr[0][1] = 0.;
    arr[1][1] = 0.;
    arr[0][2] = 0.;
    arr[1][2] = 0.;
    arrayView2d<const real64> const view = arr.toViewConst();

    CellElementStencilTPFA::IndexContainerViewConstType const &seri = tpfa.getElementRegionIndices();
    CellElementStencilTPFA::IndexContainerViewConstType const &sesri = tpfa.getElementSubRegionIndices();
    CellElementStencilTPFA::IndexContainerViewConstType const &sei = tpfa.getElementIndices();

    auto const phaseVelocityView = AccessorHelper<true>::makeElementAccessor<3>( phaseVelocity[0][0],
                                                                                 tpfa.stencilSize(0),
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0],
                                                                                 nPhases,
                                                                                 nDirs);

    wrapper.computeVelocity(0/*iconn*/, 0/*ip*/, flux, {view[0], view[1]}, phaseVelocityView.toNestedView());
    wrapper.computeVelocity(0/*iconn*/, 1/*ip*/, 100*flux, {view[0], view[1]}, phaseVelocityView.toNestedView());

    GEOS_LOG_RANK(GEOS_FMT("phase velocity : {}",phaseVelocityView));

}
//TEST( testTilter, Velocity_tilted )
//{
//}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    geos::basicSetup(argc, argv);

    int const result = RUN_ALL_TESTS();

    geos::basicCleanup();

    return result;
}
