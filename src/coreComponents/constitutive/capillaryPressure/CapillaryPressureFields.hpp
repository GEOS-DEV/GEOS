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
 * @file CapillaryPressureFields.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREFIELDS_HPP_
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREFIELDS_HPP_

#include "constitutive/capillaryPressure/Layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geosx {

    namespace fields {

        namespace cappres {

            using array2dLayoutPhase = array2d<real64, compflow::LAYOUT_PHASE>;
            using array3dLayoutCapPressure = array3d<real64, constitutive::cappres::LAYOUT_CAPPRES>;
            using array4dLayoutCapPressure_dS = array4d<real64, constitutive::cappres::LAYOUT_CAPPRES_DS>;



            enum ModeIndexType : integer {
                DRAINAGE = 0,//to be used in array of Kernels
                IMBIBITION = 1,
                DRAINAGE_TO_IMBIBITION = 2,
                IMBIBITION_TO_DRAINAGE = 3
            };

            DECLARE_FIELD(phaseCapPressure,
                          "phaseCapPressure",
                          array3dLayoutCapPressure,
                          0,
                          LEVEL_0,
                          WRITE_AND_READ,
                          "Phase capillary pressure");

            DECLARE_FIELD(dPhaseCapPressure_dPhaseVolFraction,
                          "dPhaseCapPressure_dPhaseVolFraction",
                          array4dLayoutCapPressure_dS,
                          0,
                          NOPLOT,
                          WRITE_AND_READ,
                          "Derivative of phase capillary pressure with respect to phase volume fraction");

            DECLARE_FIELD(jFuncMultiplier,
                          "jFuncMultiplier",
                          array2d<real64>,
                          0,
                          NOPLOT,
                          WRITE_AND_READ,
                          "Multiplier for the Leverett J-function");

            DECLARE_FIELD(phaseTrappedVolFraction,
                          "phaseTrappedVolumeFraction",
                          array3dLayoutCapPressure,
                          0,
                          LEVEL_0,
                          WRITE_AND_READ,
                          "Phase Trapped Volume Fraction");

            DECLARE_FIELD(mode,
                          "Hysteresis Mode",
                          array1d<integer>,
                          0,
                          LEVEL_0,
                          WRITE_AND_READ,
                          "Hysteresis mode");


            DECLARE_FIELD(phaseMaxHistoricalVolFraction,
                          "phaseMaxHistoricalVolFraction",
                          array2dLayoutPhase,
                          0,
                          LEVEL_0,
                          WRITE_AND_READ,
                          "Phase max historical phase volume fraction");

            DECLARE_FIELD(phaseMinHistoricalVolFraction,
                          "phaseMinHistoricalVolFraction",
                          array2dLayoutPhase,
                          0,
                          LEVEL_0,
                          WRITE_AND_READ,
                          "Phase min historical phase volume fraction");


        }

    }

}

#endif // GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREFIELDS_HPP_
