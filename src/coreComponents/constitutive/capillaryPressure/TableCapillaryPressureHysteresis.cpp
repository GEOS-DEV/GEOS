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


#include "TableCapillaryPressureHysteresis.hpp"

#include "constitutive/capillaryPressure/TableCapillaryPressureHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos {

    using namespace dataRepository;

    namespace constitutive {

        TableCapillaryPressureHysteresis::TableCapillaryPressureHysteresis(const std::string &name,
                                                                           dataRepository::Group *const parent)
                : CapillaryPressureBase(name, parent) {

            registerWrapper(viewKeyStruct::phaseHasHysteresisString(), &m_phaseHasHysteresis).
                            setInputFlag(InputFlags::FALSE)
                    .                         // will be deduced from tables
                            setSizedFromParent(0);

            registerWrapper(viewKeyStruct::landParameterString(), &m_landParam).
                    setInputFlag(InputFlags::FALSE).                       // will be deduced from tables
                    setSizedFromParent(0);

            //2phase
            registerWrapper(viewKeyStruct::drainageWettingNonWettingCapPresTableNameString(),
                            &m_drainageWettingNonWettingCapPresTableName).
                    setInputFlag(InputFlags::OPTIONAL).
                    setDescription("Name of the drainage two-phase table for capillary pressure curve. \n"
                                   "If you want to use 3-phase flow please use instead " +
                                   string(viewKeyStruct::drainageWettingIntermediateCapPresTableNameString()) +
                                   " and " +
                                   string(viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString()) +
                                   "to specify the tables names");
            registerWrapper(viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(),
                            &m_imbibitionWettingNonWettingCapPresTableName).
                    setInputFlag(InputFlags::OPTIONAL).
                    setDescription("Name of the drainage two-phase table for capillary pressure curve. \n"
                                   "If you want to use 3-phase flow please use instead " +
                                   string(viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString()) +
                                   " and " +
                                   string(viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString()) +
                                   "to specify the tables names");
            //3phase
            registerWrapper(viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),
                            &m_drainageWettingIntermediateCapPresTableName).
                    setInputFlag(InputFlags::OPTIONAL).
                    setDescription(
                    "Drainage wetting/intermediate (e.g. w/o) capillary pressure table name for the wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves");
            registerWrapper(viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString(),
                            &m_drainageNonWettingIntermediateCapPresTableName).
                    setInputFlag(InputFlags::OPTIONAL).
                    setDescription(
                    "Drainage non-wetting/intermediate (e.g. o/g) capillary pressure table name for the non-wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves");
            registerWrapper(viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(),
                            &m_imbibitionWettingIntermediateCapPresTableName).
                    setInputFlag(InputFlags::OPTIONAL).
                    setDescription("Imbibition wetting/intermediate (e.g. w/o) table name for the wetting phase.\n"
                                   "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves");
            registerWrapper(viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString(),
                            &m_imbibitionNonWettingIntermediateCapPresTableName).
                    setInputFlag(InputFlags::OPTIONAL).
                    setDescription("Imbibition non-wetting/intermediate (e.g. o/g) table name for the wetting phase.\n"
                                   "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves");

            // kernels
            //2p
            registerWrapper(viewKeyStruct::wettingNonWettingCapillaryPressureKernelWrappersString(),
                            &m_wettingNonWettingCapillaryPressureKernelWrappers)
                    .setSizedFromParent(0).setRestartFlags(RestartFlags::NO_WRITE);
            //3p
            registerWrapper(viewKeyStruct::wettingIntermediateCapillaryPressureKernelWrappersString(),
                            &m_wettingIntermediateCapillaryPressureKernelWrappers)
                    .setSizedFromParent(0).setRestartFlags(RestartFlags::NO_WRITE);
            registerWrapper(viewKeyStruct::nonWettingIntermediateCapillaryPressureKernelWrappersString(),
                            &m_nonWettingIntermediateCapillaryPressureKernelWrappers)
                    .setSizedFromParent(0).setRestartFlags(RestartFlags::NO_WRITE);


            registerWrapper(viewKeyStruct::wettingCurveString(), &m_wettingCurve).
                            setInputFlag(
                            InputFlags::FALSE).         // will be deduced from tables
                            setSizedFromParent(
                            0)
                    .setRestartFlags(RestartFlags::NO_WRITE);

            registerWrapper(viewKeyStruct::nonWettingCurveString(), &m_nonWettingCurve).
                            setInputFlag(
                            InputFlags::FALSE).         // will be deduced from tables
                            setSizedFromParent(
                            0)
                    .setRestartFlags(RestartFlags::NO_WRITE);

            //Forwarded to KilloughHysteresis
            registerWrapper(KilloughHysteresis::viewKeyStruct::jerauldParameterAString(), &m_jerauldParam_a).
                    setInputFlag(InputFlags::OPTIONAL).
                    setApplyDefaultValue(0.1).
                    setDescription(
                    "First parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation).");

            registerWrapper(KilloughHysteresis::viewKeyStruct::jerauldParameterBString(), &m_jerauldParam_b).
                    setInputFlag(InputFlags::OPTIONAL).
                    setApplyDefaultValue(0.0).
                    setDescription(
                    "Second parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation).");


            registerWrapper(KilloughHysteresis::viewKeyStruct::killoughCurvatureParameterPcString(),
                            &m_killoughCurvatureParamCapPres).
                    setInputFlag(
                    InputFlags::OPTIONAL).
                    setApplyDefaultValue(
                    .1).
                    setDescription(
                    "Curvature parameter introduced by Killough for wetting-phase hysteresis (see RTD documentation).");

            //misc
            registerWrapper(viewKeyStruct::phaseIntermediateMinVolFractionString(), &m_phaseIntermediateMinVolFraction).
                    setInputFlag(InputFlags::FALSE).setDescription("min vol fraction of intermediate if exist").
                    // will be deduced from tables
                    setSizedFromParent(0);

            registerField< fields::cappres::mode >( &m_mode );


            registerField< fields::cappres::phaseMaxHistoricalVolFraction >(
                           &m_phaseMaxHistoricalVolFraction );
            registerField< fields::cappres::phaseMinHistoricalVolFraction >(
                           &m_phaseMinHistoricalVolFraction );

        }

/// usual utils

        void TableCapillaryPressureHysteresis::postProcessInput() {

            using TPP = ThreePhasePairPhaseType;

            integer const numPhases = m_phaseNames.size();
            GEOS_THROW_IF(numPhases != 2 && numPhases != 3,
                           GEOS_FMT("{}: the expected number of fluid phases is either two, or three",
                                     getFullName()),
                           InputError);

            m_phaseHasHysteresis.resize(2);

            if (numPhases == 2) {
                GEOS_THROW_IF(m_drainageWettingNonWettingCapPresTableName.empty(),
                               GEOS_FMT(
                                       "{}: for a two-phase flow simulation, we must use {} to specify the capillary pressure table for the drainage pair (wetting phase, non-wetting phase)",
                                       getFullName(),
                                       viewKeyStruct::drainageWettingNonWettingCapPresTableNameString()),
                               InputError);


                m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = (m_imbibitionWettingNonWettingCapPresTableName.empty() ||
                                                                   m_imbibitionWettingNonWettingCapPresTableName ==
                                                                   m_drainageWettingNonWettingCapPresTableName)
                                                                  ? 0 : 1;
                m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING];


            } else if (numPhases == 3) {


                GEOS_THROW_IF(m_drainageWettingIntermediateCapPresTableName.empty() ||
                               m_drainageNonWettingIntermediateCapPresTableName.empty(),
                               GEOS_FMT(
                                       "{}: for a three-phase flow simulation, we must use {} to specify the capillary pressure table "
                                       "for the pair (wetting phase, intermediate phase), and {} to specify the capillary pressure table "
                                       "for the pair (non-wetting phase, intermediate phase)",
                                       getFullName(),
                                       viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),
                                       viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString()),
                               InputError);

                m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = (m_imbibitionWettingIntermediateCapPresTableName.empty() ||
                                                                   m_imbibitionWettingIntermediateCapPresTableName ==
                                                                   m_drainageWettingIntermediateCapPresTableName)
                                                                  ? 0 : 1;

                m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = (m_imbibitionNonWettingIntermediateCapPresTableName.empty() ||
                                                                      m_imbibitionNonWettingIntermediateCapPresTableName ==
                                                                      m_drainageNonWettingIntermediateCapPresTableName)
                                                                     ? 0 : 1;
            }
            //Killough section
            //TODO improve hard coded default
            KilloughHysteresis::postProcessInput(m_jerauldParam_a, m_jerauldParam_b, 0,
                                                 m_killoughCurvatureParamCapPres);

            GEOS_THROW_IF(m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] == 0 &&
                           m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] == 0,
                           GEOS_FMT(
                                   "{}: we must use {} (2-phase) / {} or {} (3-phase) to specify at least one imbibition relative permeability table",
                                   getFullName(),
                                   viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(),
                                   viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(),
                                   viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString()),
                           InputError);

        }

        void TableCapillaryPressureHysteresis::initializePreSubGroups() {
            CapillaryPressureBase::initializePreSubGroups();

            integer const numPhases = m_phaseNames.size();
            FunctionManager const &functionManager = FunctionManager::getInstance();

            //equivalent to oil/gas - a.k.a two phase flow ordered by non wetting
            bool const capPresMustBeIncreasing = (m_phaseOrder[PhaseType::WATER] < 0)
                                                 ? true   // pc on the gas phase, function must be increasing
                                                 : false; // pc on the water phase, function must be decreasing


            // Step 1: check sanity of drainage tables
            if (numPhases == 2) {

                real64 drainageWettingPhaseMaxVolumeFraction, drainageWettingMinCapPres,
                        drainageNonWettingPhaseMinVolumeFraction, drainageNonWettingMinCapPres,
                        imbibitionWettingPhaseMaxVolumeFraction, imbibitionWettingMinCapPres,
                        imbibitionNonWettingPhaseMinVolumeFraction, imbibitionNonWettingMinCapPres,
                        wettingPhaseMinVolumeFraction, wettingMaxCapPres,
                        nonWettingPhaseMaxVolumeFraction, nonWettingMaxCapPres;

                {
		  
		  imbibitionNonWettingMinCapPres = 0.0; 

                    GEOS_THROW_IF(!functionManager.hasGroup(m_drainageWettingNonWettingCapPresTableName),
                                   GEOS_FMT("{}: the table function named {} could not be found",
                                             getFullName(),
                                             m_drainageWettingNonWettingCapPresTableName),
                                   InputError);
                    TableFunction const
                            &capPresTable = functionManager.getGroup<TableFunction>(
                            m_drainageWettingNonWettingCapPresTableName);

                    //w/o  or  w/g pair
                    if (!capPresMustBeIncreasing) {
                        TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTable, getFullName(),
                                                                                      capPresMustBeIncreasing,
                                                                                      drainageWettingPhaseMaxVolumeFraction,
                                                                                      wettingPhaseMinVolumeFraction,
                                                                                      drainageWettingMinCapPres,
                                                                                      wettingMaxCapPres);

                        drainageNonWettingPhaseMinVolumeFraction = 1. - drainageWettingPhaseMaxVolumeFraction;
                        nonWettingPhaseMaxVolumeFraction = 1. - wettingPhaseMinVolumeFraction;

                    } else { // o/g pair
                        TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTable, getFullName(),
                                                                                      capPresMustBeIncreasing,
                                                                                      nonWettingPhaseMaxVolumeFraction,
                                                                                      drainageNonWettingPhaseMinVolumeFraction,
                                                                                      nonWettingMaxCapPres,
                                                                                      drainageNonWettingMinCapPres );

                        drainageWettingPhaseMaxVolumeFraction = 1. - drainageNonWettingPhaseMinVolumeFraction;
                        wettingPhaseMinVolumeFraction = 1. - nonWettingPhaseMaxVolumeFraction;
                    }

                }

                {
                    GEOS_THROW_IF(!functionManager.hasGroup(m_imbibitionWettingNonWettingCapPresTableName),
                                   GEOS_FMT("{}: the table function named {} could not be found",
                                             getFullName(),
                                             m_imbibitionWettingNonWettingCapPresTableName),
                                   InputError);
                    TableFunction const
                            &capPresTable = functionManager.getGroup<TableFunction>(
                            m_imbibitionWettingNonWettingCapPresTableName);

                    //w/o  or  w/g pair
                    if (!capPresMustBeIncreasing) {
                        TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTable, getFullName(),
                                                                                      capPresMustBeIncreasing,
                                                                                      imbibitionWettingPhaseMaxVolumeFraction,
                                                                                      wettingPhaseMinVolumeFraction,
                                                                                      imbibitionWettingMinCapPres,
                                                                                      wettingMaxCapPres);

                        imbibitionNonWettingPhaseMinVolumeFraction = 1. - imbibitionWettingPhaseMaxVolumeFraction;
                        nonWettingPhaseMaxVolumeFraction = 1. - wettingPhaseMinVolumeFraction;

                    } else { // o/g pair
                        TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTable, getFullName(),
                                                                                      capPresMustBeIncreasing,
                                                                                      nonWettingPhaseMaxVolumeFraction,
                                                                                      imbibitionNonWettingPhaseMinVolumeFraction,
                                                                                      nonWettingMaxCapPres,
                                                                                      imbibitionWettingMinCapPres );

                        imbibitionWettingPhaseMaxVolumeFraction = 1. - imbibitionNonWettingPhaseMinVolumeFraction;
                        wettingPhaseMinVolumeFraction = 1. - nonWettingPhaseMaxVolumeFraction;
                    }
                }

                //constructing wetting/nonwetting curves

                if(!capPresMustBeIncreasing) {
                    m_wettingCurve.setPoints(
                            {wettingPhaseMinVolumeFraction, wettingMaxCapPres},   // same as imbibition min
                            {imbibitionWettingPhaseMaxVolumeFraction, imbibitionWettingMinCapPres},
                            {drainageWettingPhaseMaxVolumeFraction, drainageWettingMinCapPres});
                }
                else {
                    m_nonWettingCurve.setPoints(
                            {nonWettingPhaseMaxVolumeFraction,nonWettingMaxCapPres},
                            {imbibitionNonWettingPhaseMinVolumeFraction,imbibitionNonWettingMinCapPres},
                            {drainageNonWettingPhaseMinVolumeFraction,drainageNonWettingMinCapPres}
                            );
                }

            } else if (numPhases == 3) {

	      real64 drainageWettingPhaseMaxVolumeFraction, drainageWettingMinCapPres,
		drainageNonWettingPhaseMinVolumeFraction, drainageNonWettingMinCapPres,
		imbibitionWettingPhaseMaxVolumeFraction, imbibitionWettingMinCapPres,
		imbibitionNonWettingPhaseMinVolumeFraction, imbibitionNonWettingMinCapPres,
		wettingPhaseMinVolumeFraction, wettingMaxCapPres,
		nonWettingPhaseMaxVolumeFraction, nonWettingMaxCapPres;

	      GEOS_UNUSED_VAR( drainageWettingMinCapPres );	      
//define scope to avoid differentiate temp var (lazy)
                {
                    GEOS_THROW_IF(!functionManager.hasGroup(m_drainageWettingIntermediateCapPresTableName),
                                   GEOS_FMT("{}: the table function named {} could not be found",
                                             getFullName(),
                                             m_drainageWettingIntermediateCapPresTableName),
                                   InputError);
                    TableFunction const
                            &capPresTableWI = functionManager.getGroup<TableFunction>(
                            m_drainageWettingIntermediateCapPresTableName);
                    TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTableWI, getFullName(), false,
                                                                                  drainageWettingPhaseMaxVolumeFraction,
                                                                                  wettingPhaseMinVolumeFraction,
                                                                                  drainageNonWettingMinCapPres,
                                                                                  wettingMaxCapPres);

                    GEOS_THROW_IF(!functionManager.hasGroup(m_drainageNonWettingIntermediateCapPresTableName),
                                   GEOS_FMT("{}: the table function named {} could not be found",
                                             getFullName(),
                                             m_drainageNonWettingIntermediateCapPresTableName),
                                   InputError);
                    TableFunction const &capPresTableNWI =
                            functionManager.getGroup<TableFunction>(m_drainageNonWettingIntermediateCapPresTableName);
                    TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTableNWI, getFullName(), true,
                                                                                  nonWettingPhaseMaxVolumeFraction,
                                                                                  drainageNonWettingPhaseMinVolumeFraction,
                                                                                  nonWettingMaxCapPres,
                                                                                  drainageWettingPhaseMaxVolumeFraction
                    );

                    m_phaseIntermediateMinVolFraction =
                            1.0 - drainageWettingPhaseMaxVolumeFraction - drainageWettingPhaseMaxVolumeFraction;
                }

                if (!m_imbibitionWettingIntermediateCapPresTableName.empty()) {

                    GEOS_THROW_IF(!functionManager.hasGroup(m_imbibitionWettingIntermediateCapPresTableName),
                                   GEOS_FMT("{}: the table function named {} could not be found",
                                             getFullName(),
                                             m_imbibitionWettingIntermediateCapPresTableName),
                                   InputError);
                    TableFunction const
                            &capPresTableWI = functionManager.getGroup<TableFunction>(
                            m_imbibitionWettingIntermediateCapPresTableName);
                    TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTableWI, getFullName(), false,
                                                                                  imbibitionWettingPhaseMaxVolumeFraction,
                                                                                  wettingPhaseMinVolumeFraction,
                                                                                  imbibitionWettingMinCapPres,
                                                                                  wettingMaxCapPres
                                                                                  );


                }

                if (!m_imbibitionNonWettingIntermediateCapPresTableName.empty()) {

                    GEOS_THROW_IF(!functionManager.hasGroup(m_imbibitionNonWettingIntermediateCapPresTableName),
                                   GEOS_FMT("{}: the table function named {} could not be found",
                                             getFullName(),
                                             m_imbibitionNonWettingIntermediateCapPresTableName),
                                   InputError);
                    TableFunction const &capPresTableNWI =
                            functionManager.getGroup<TableFunction>(m_imbibitionNonWettingIntermediateCapPresTableName);
                    TableCapillaryPressureHelpers::validateCapillaryPressureTable(capPresTableNWI, getFullName(), true,
                                                                                  nonWettingPhaseMaxVolumeFraction,
                                                                                  imbibitionNonWettingPhaseMinVolumeFraction,
                                                                                  nonWettingMaxCapPres,
                                                                                  imbibitionNonWettingMinCapPres);


                }
            }

            // Step 2: check the sanity btw drainage and imbibition
            auto const eps = 1e-15;
            if (numPhases == 2) {
                //TODO weak make stronger
                GEOS_THROW_IF(
                        m_wettingCurve.isZero() && m_nonWettingCurve.isZero(),
                        GEOS_FMT(
                                "{}: Inconsistent data for capillary pressure hysteresis. No hysteresis curve is defined.",
                                getFullName()),
                        InputError);

                GEOS_THROW_IF(
                        !m_wettingCurve.isZero() && !m_nonWettingCurve.isZero(),
                        GEOS_FMT(
                                "{}: Inconsistent data for capillary pressure hysteresis. Both non wetting and wetting hysteresis curve are defined in two phase flow setting.",
                                getFullName()),
                        InputError);


            } else if (numPhases == 3) {

                GEOS_THROW_IF(std::fabs(m_wettingCurve.oppositeBoundPhaseVolFraction - (1. - m_nonWettingCurve.oppositeBoundPhaseVolFraction - m_phaseIntermediateMinVolFraction)) > eps,
                               GEOS_FMT(
                                       "{}: Inconsistent data for capillary pressure hysteresis. {}, {} and {} should sum up to 1.",
                                       getFullName(), "Sw_min", "Snw_max", "Sinter_min"),
                               InputError);
                GEOS_THROW_IF(std::fabs(m_wettingCurve.drainageExtremaPhaseVolFraction - (1. - m_nonWettingCurve.drainageExtremaPhaseVolFraction - m_phaseIntermediateMinVolFraction)) > eps,
                               GEOS_FMT(
                                       "{}: Inconsistent data for capillary pressure hysteresis. {}, {} and {} should sum up to 1.",
                                       getFullName(), "Sw_min", "Snw_max", "Sinter_min"),
                               InputError);
                GEOS_THROW_IF(std::fabs(m_wettingCurve.imbibitionExtremaPhaseVolFraction - (1. - m_nonWettingCurve.imbibitionExtremaPhaseVolFraction - m_phaseIntermediateMinVolFraction)) > eps,
                               GEOS_FMT(
                                       "{}: Inconsistent data for capillary pressure hysteresis. {}, {} and {} should sum up to 1.",
                                       getFullName(), "Sw_min", "Snw_max", "Sinter_min"),
                               InputError);

            }


            // Step 3: compute the Land coefficient
            computeLandCoefficient();
            
            // Step 4: Ensure arrays are properly resized if they weren't resized earlier
            // This can happen if resizeFields() was called before phase names were set
            if (m_phaseMaxHistoricalVolFraction.size(1) == 0 && numPhases > 0) {
                localIndex const currentSize = m_phaseMaxHistoricalVolFraction.size(0);
                if (currentSize > 0) {
                    m_phaseMaxHistoricalVolFraction.resize(currentSize, numPhases);
                    m_phaseMinHistoricalVolFraction.resize(currentSize, numPhases);
                    m_phaseMaxHistoricalVolFraction.setValues<parallelDevicePolicy<> >(0.0);
                    m_phaseMinHistoricalVolFraction.setValues<parallelDevicePolicy<> >(1.0);
                }
            }
        }

/// Land coeff (tb refactored out in KilloughHysteresis) and saved cvgd

        void TableCapillaryPressureHysteresis::computeLandCoefficient() {
            // For now, we keep two separate Land parameters for the wetting and non-wetting phases
            // For two-phase flow, we make sure that they are equal
            m_landParam.resize(2);

            // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

            // Step 1: Land parameter for the wetting phase

            integer ipWetting, ipNonWetting;
            std::tie(ipWetting, ipNonWetting) = phaseIndex(m_phaseOrder);

            KilloughHysteresis::computeLandCoefficient( m_wettingCurve, m_landParam[ipWetting] );
            KilloughHysteresis::computeLandCoefficient( m_nonWettingCurve, m_landParam[ipNonWetting] );

        }

/// common utils
        void TableCapillaryPressureHysteresis::resizeFields(localIndex const size, localIndex const numPts) {
            CapillaryPressureBase::resizeFields(size, numPts);

            integer const numPhases = numFluidPhases();
            
            // If phase names haven't been set yet, we still need to resize m_mode
            // The phase arrays will be resized properly once phase names are set
            m_mode.resize(size);
            
            if (numPhases > 0) {
                m_phaseMaxHistoricalVolFraction.resize(size, numPhases);
                m_phaseMinHistoricalVolFraction.resize(size, numPhases);
                m_phaseMaxHistoricalVolFraction.setValues<parallelDevicePolicy<> >(0.0);
                m_phaseMinHistoricalVolFraction.setValues<parallelDevicePolicy<> >(1.0);
            }
            // If numPhases == 0, the arrays will remain uninitialized until resizeFields is called again
            // after phase names are set. This should happen automatically during initialization.
        }

        void TableCapillaryPressureHysteresis::saveConvergedPhaseVolFractionState(
                arrayView2d<real64 const, compflow::USD_PHASE> const &phaseVolFraction) const {
            CapillaryPressureBase::saveConvergedState();

            arrayView2d<real64, compflow::USD_PHASE> phaseMaxHistoricalVolFraction = m_phaseMaxHistoricalVolFraction.toView();
            arrayView2d<real64, compflow::USD_PHASE> phaseMinHistoricalVolFraction = m_phaseMinHistoricalVolFraction.toView();

            localIndex const numElems = phaseVolFraction.size(0);
            integer const numPhases = numFluidPhases();

            forAll<parallelDevicePolicy<> >(numElems, [=] GEOS_HOST_DEVICE(localIndex const ei) {
                for (integer ip = 0; ip < numPhases; ++ip) {
                    phaseMaxHistoricalVolFraction[ei][ip] = LvArray::math::max(phaseVolFraction[ei][ip],
                                                                               phaseMaxHistoricalVolFraction[ei][ip]);
                    phaseMinHistoricalVolFraction[ei][ip] = LvArray::math::min(phaseVolFraction[ei][ip],
                                                                               phaseMinHistoricalVolFraction[ei][ip]);
                }
            });

        }

        void
        TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionWettingCapillaryPressure(
                const arrayView1d<const TableFunction::KernelWrapper> &wettingKernelWapper,
                const KilloughHysteresis::HysteresisCurve &wettingCurve,
                const KilloughHysteresis::HysteresisCurve &nonWettingCurve, //discard if not needed
                const geos::real64 &landParam,
                const geos::real64 &phaseVolFraction,
                const geos::real64 &phaseMinHistoricalVolFraction,
                geos::real64 &phaseTrappedVolFrac,
                geos::real64 &phaseCapPressure,
                geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                const ModeIndexType &mode) const {
            GEOS_ASSERT(wettingCurve.isWetting());
            real64 const S = phaseVolFraction;
            real64 const Smxi = wettingCurve.imbibitionExtremaPhaseVolFraction;
            real64 const Smxd = wettingCurve.drainageExtremaPhaseVolFraction;
            real64 const Smin = wettingCurve.oppositeBoundPhaseVolFraction;

	    GEOS_UNUSED_VAR( Smxi, Smxd, phaseTrappedVolFrac );	      
	    
//  if( S <= Smin )
//  {
//    //below accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else if( S >= Smxd )
//  {
//    //above accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else
            {
                //drainage to imbibition
                real64 dpci_dS, dpcd_dS;
                real64 const pci = wettingKernelWapper[ModeIndexType::IMBIBITION].compute(&S, &dpci_dS);
                real64 const pcd = wettingKernelWapper[ModeIndexType::DRAINAGE].compute(&S, &dpcd_dS);
                real64 const Somin = m_phaseIntermediateMinVolFraction;

                // Step 1: get the trapped from wetting data
                real64 const Shy = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;

                real64 const E = m_killoughCurvatureParamCapPres;

                //Step 2. compute F as in (EQ 34.15) F = (1/(Sw-Shy+E)-1/E) / (1/(Swma-Shy+E)-1/E)
                //drainage to imbibition branch
                if (mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {

                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(nonWettingCurve,
                                                                               Shy,
                                                                               landParam,
                                                                               m_jerauldParam_a,
                                                                               m_jerauldParam_b,
                                                                               Scrt);
                    real64 const Swma = 1 - Scrt - Somin;
                    real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Swma - Shy + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. Eventually assemble everything following (EQ. 34.14)
                    phaseCapPressure = pcd + F * (pci - pcd);
                    dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpci_dS - dpcd_dS);
                }
                    //imbibition to drainage
                else if (mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(wettingCurve,
                                                                               Shy,
                                                                               landParam,
                                                                               m_jerauldParam_a,
                                                                               m_jerauldParam_b,
                                                                               Scrt);

                    real64 F = (1. / (Shy - S + E) - 1. / E) / (1. / (Shy - Scrt + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. Eventually assemble everything following (EQ. 34.14)
                    phaseCapPressure = pci + F * (pcd - pci);
                    dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpcd_dS - dpci_dS);


                } else {
                    GEOS_THROW(GEOS_FMT("{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                                          "TableCapillaryPressureHysteresis",
                                          (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                                             ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                                InputError);
                }


            }

        }

        void
        TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionWettingCapillaryPressure(
                const arrayView1d<const TableFunction::KernelWrapper> &wettingKernelWapper,
                const KilloughHysteresis::HysteresisCurve &wettingCurve,
                const geos::real64 &landParam,
                const geos::real64 &phaseVolFraction,
                const geos::real64 &phaseMinHistoricalVolFraction,
                geos::real64 &phaseTrappedVolFrac,
                geos::real64 &phaseCapPressure,
                geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                const ModeIndexType &mode) const {
            GEOS_ASSERT(wettingCurve.isWetting());
            real64 const S = phaseVolFraction;
            real64 const Smxi = wettingCurve.imbibitionExtremaPhaseVolFraction;
            real64 const Smxd = wettingCurve.drainageExtremaPhaseVolFraction;
            real64 const Smin = wettingCurve.oppositeBoundPhaseVolFraction;

	    GEOS_UNUSED_VAR( Smxi, Smxd, phaseTrappedVolFrac );	      
	    
//  if( S <= Smin )
//  {
//    //below accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else if( S >= Smxd )
//  {
//    //above accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else
            {
                //drainage to imbibition
                real64 dpci_dS, dpcd_dS;
                real64 const pci = wettingKernelWapper[ModeIndexType::IMBIBITION].compute(&S, &dpci_dS);
                real64 const pcd = wettingKernelWapper[ModeIndexType::DRAINAGE].compute(&S, &dpcd_dS);

                // Step 1: get the trapped from wetting data
                real64 const Shy = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;

                real64 const E = m_killoughCurvatureParamCapPres;

                //Step 2. compute F as in (EQ 34.15) F = (1/(Sw-Shy+E)-1/E) / (1/(Swma-Shy+E)-1/E)
                //drainage to imbibition branch
                if (mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {

                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(wettingCurve,
                                                                               Shy,
                                                                               landParam,
                                                                               m_jerauldParam_a,
                                                                               m_jerauldParam_b,
                                                                               Scrt);



                    //should be the pore space accessible to the two wetting phase
                    real64 const Swma = 1 - (1 -  Scrt);
                    real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Swma - Shy + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. Eventually assemble everything following (EQ. 34.14)
                    phaseCapPressure = pcd + F * (pci - pcd);
                    dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpci_dS - dpcd_dS);
                }
                    //imbibition to drainage
                else if (mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(wettingCurve,
                                                                               Shy,
                                                                               landParam,
                                                                               m_jerauldParam_a,
                                                                               m_jerauldParam_b,
                                                                               Scrt);

                    real64 F = (1. / (Shy - S + E) - 1. / E) / (1. / (Shy - Scrt + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. Eventually assemble everything following (EQ. 34.14)
                    phaseCapPressure = pci + F * (pcd - pci);
                    dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpcd_dS - dpci_dS);


                } else {
                    GEOS_THROW(GEOS_FMT("{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                                          "TableCapillaryPressureHysteresis",
                                          (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                                             ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                                InputError);
                }


            }

        }
        void TableCapillaryPressureHysteresis::KernelWrapper::computeTwoPhaseWetting(const geos::integer ipWetting,
                                                                                     const geos::integer GEOS_UNUSED_PARAM( ipNonWetting ),
                                                                                     const arraySlice1d<const geos::real64,
                                                                                             compflow::USD_PHASE -
                                                                                             1> &phaseVolFraction,
                                                                                     const arraySlice1d<const geos::real64,
                                                                                             compflow::USD_PHASE
                                                                                             -
                                                                                             1> &phaseMaxHistoricalVolFraction,
                                                                                     const arraySlice1d<const geos::real64,
                                                                                             compflow::USD_PHASE
                                                                                             -
                                                                                             1> &phaseMinHistoricalVolFraction,
                                                                                     const arraySlice1d<geos::real64,
                                                                                             relperm::USD_RELPERM
                                                                                             - 2> &phaseTrappedVolFrac,
                                                                                     arraySlice1d<geos::real64,
                                                                                             relperm::USD_RELPERM
                                                                                             -
                                                                                             2> const &phaseCapPressure,
                                                                                     arraySlice2d<geos::real64,
                                                                                             relperm::USD_RELPERM_DS
                                                                                             -
                                                                                             2> const &dPhaseCapPressure_dPhaseVolFrac,
                                                                                     ModeIndexType &mode) const {
            using TTP = ThreePhasePairPhaseType;

            // Validate array sizes and indices before accessing
            GEOS_ASSERT_MSG(ipWetting >= 0, "ipWetting must be non-negative");
            GEOS_ASSERT_MSG(static_cast<integer>(phaseVolFraction.size()) > ipWetting,
                        GEOS_FMT("phaseVolFraction array too small: size={}, ipWetting={}. "
                                 "This usually means the arrays haven't been properly resized. "
                                 "Ensure resizeFields() has been called before using the KernelWrapper.",
                                 phaseVolFraction.size(), ipWetting));
            GEOS_ASSERT_MSG(static_cast<integer>(phaseMaxHistoricalVolFraction.size()) > ipWetting,
                        GEOS_FMT("phaseMaxHistoricalVolFraction array too small: size={}, ipWetting={}. "
                                 "This usually means the arrays haven't been properly resized. "
                                 "Ensure resizeFields() has been called before using the KernelWrapper.",
                                 phaseMaxHistoricalVolFraction.size(), ipWetting));
            GEOS_ASSERT_MSG(static_cast<integer>(phaseMinHistoricalVolFraction.size()) > ipWetting,
                        GEOS_FMT("phaseMinHistoricalVolFraction array too small: size={}, ipWetting={}. "
                                 "This usually means the arrays haven't been properly resized. "
                                 "Ensure resizeFields() has been called before using the KernelWrapper.",
                                 phaseMinHistoricalVolFraction.size(), ipWetting));
            GEOS_ASSERT_MSG(static_cast<integer>(m_wettingNonWettingCapillaryPressureKernelWrappers.size()) >= 2,
                        GEOS_FMT("m_wettingNonWettingCapillaryPressureKernelWrappers must have at least 2 elements, but got {}. "
                                 "This usually means createAllTableKernelWrappers() failed to populate the arrays.",
                                 m_wettingNonWettingCapillaryPressureKernelWrappers.size()));

            // Determine mode based on saturation condition (matching relative permeability logic)
            // Use DRAINAGE when saturation is at or below minimum, use scanning curves when above minimum
            // This matches the relative permeability logic: drainage when S <= S_min, imbibition (scanning) when S > S_min
            bool const useDrainage = !m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] ||
                                     phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer;
            
            //--- wetting  cap pressure -- W/O or W/G two phase flow
            // Use drainage curve when S <= S_min (matching relative permeability logic)
            // Use scanning curves when S > S_min (matching relative permeability logic)
            // DEBUG: Print mode for capillary pressure
            // printf("CapPressure: mode=%d, S_w=%.6e, S_min=%.6e, S_max=%.6e, hasHyst=%d\n",
            //        static_cast<int>(mode),
            //        phaseVolFraction[ipWetting],
            //        phaseMinHistoricalVolFraction[ipWetting],
            //        phaseMaxHistoricalVolFraction[ipWetting],
            //        static_cast<int>(m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING]));
            
            if (useDrainage) {
                // Use simple drainage curve (matching relative permeability)
                mode = ModeIndexType::DRAINAGE;
                phaseTrappedVolFrac[ipWetting] = LvArray::math::min(phaseVolFraction[ipWetting],
                                                                    m_wettingCurve.oppositeBoundPhaseVolFraction);
                // printf("CapPressure: Using DRAINAGE curve (arrayIndex=0)\n");
                GEOS_ASSERT_MSG(static_cast<integer>(m_wettingNonWettingCapillaryPressureKernelWrappers.size()) > ModeIndexType::DRAINAGE,
                                "Invalid array index for kernel wrapper access");
                computeBoundCapillaryPressure(
                        m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE],
                        phaseVolFraction[ipWetting],
                        phaseCapPressure[ipWetting],
                        dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting]);
            }
            else {
                // Use scanning curves (matching relative permeability)
                // When S > S_min, saturation is increasing (moving away from minimum), 
                // so we're transitioning from drainage to imbibition -> use DRAINAGE_TO_IMBIBITION mode
                mode = ModeIndexType::DRAINAGE_TO_IMBIBITION;
                // printf("CapPressure: Using IMBIBITION (scanning) curves (mode=%d)\n", static_cast<int>(mode));
                computeImbibitionWettingCapillaryPressure(m_wettingNonWettingCapillaryPressureKernelWrappers,
                                                          m_wettingCurve,
                                                          m_landParam[ipWetting],
                                                          phaseVolFraction[ipWetting],
                                                          phaseMinHistoricalVolFraction[ipWetting],
                                                          phaseTrappedVolFrac[ipWetting],
                                                          phaseCapPressure[ipWetting],
                                                          dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting],
                                                          mode);
            }

// trapped vol fraction
            if (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {



                    real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
                    real64 const Shy = (phaseMinHistoricalVolFraction[ipWetting] < Smin)
                                       ? phaseMinHistoricalVolFraction[ipWetting]
                                       : Smin;
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_wettingCurve, Shy,
                                                                               m_landParam[ipWetting],
                                                                               m_jerauldParam_a,
                                                                               m_jerauldParam_b,
                                                                               Scrt);
                phaseTrappedVolFrac[ipWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipWetting]);


                //keep the same Land coeff as two phase only
                KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_wettingCurve.toNonWetting(), Shy,
                                                                           m_landParam[ipWetting],
                                                                           m_jerauldParam_a,
                                                                           m_jerauldParam_b,
                                                                           Scrt);
                phaseTrappedVolFrac[ipWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipWetting]);




            }
            else if (mode == ModeIndexType::IMBIBITION || mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {

                    real64 const Smax = m_wettingCurve.imbibitionExtremaPhaseVolFraction;
                    real64 const Shy = (phaseMaxHistoricalVolFraction[ipWetting] < Smax)
                                       ? phaseMaxHistoricalVolFraction[ipWetting]
                                       : Smax;
                    real64 Scrt = 0.0;
                    //TODO (jacques) check if still accurate
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_wettingCurve,
                                                                               Shy,
                                                                               m_landParam[ipWetting],
                                                                               m_jerauldParam_a,
                                                                               m_jerauldParam_b,
                                                                               Scrt);

                    phaseTrappedVolFrac[ipWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipWetting]);

                //keep the same Land coeff as two phase only
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_wettingCurve.toNonWetting(), Shy,
                                                                           m_landParam[ipWetting],
                                                                           m_jerauldParam_a,
                                                                           m_jerauldParam_b,
                                                                           Scrt);
                phaseTrappedVolFrac[ipWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipWetting]);

            }

        }

        void TableCapillaryPressureHysteresis::KernelWrapper::computeTwoPhaseNonWetting(const geos::integer ipWetting,
                                                                                        const geos::integer ipNonWetting,
                                                                                        const arraySlice1d<const geos::real64,
                                                                                                compflow::USD_PHASE -
                                                                                                1> &phaseVolFraction,
                                                                                        const arraySlice1d<const geos::real64,
                                                                                                compflow::USD_PHASE
                                                                                                -
                                                                                                1> &phaseMaxHistoricalVolFraction,
                                                                                        const arraySlice1d<const geos::real64,
                                                                                                compflow::USD_PHASE
                                                                                                -
                                                                                                1> &phaseMinHistoricalVolFraction,
                                                                                        const arraySlice1d<geos::real64,
                                                                                                relperm::USD_RELPERM
                                                                                                -
                                                                                                2> &phaseTrappedVolFrac,
                                                                                        arraySlice1d<geos::real64,
                                                                                                relperm::USD_RELPERM
                                                                                                -
                                                                                                2> const &phaseCapPressure,
                                                                                        arraySlice2d<geos::real64,
                                                                                                relperm::USD_RELPERM_DS
                                                                                                -
                                                                                                2> const &dPhaseCapPressure_dPhaseVolFrac,
                                                                                        ModeIndexType &mode) const {
            using TTP = ThreePhasePairPhaseType;
            //update state
            // TODO check if we can get rid of  DRAINAGE_TO_IMBIBITION && IMBIBITION_TO_DRAINAGE
            if (mode == ModeIndexType::DRAINAGE_TO_IMBIBITION &&
                phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] + flowReversalBuffer)
                mode = ModeIndexType::DRAINAGE;
            if (mode == ModeIndexType::IMBIBITION_TO_DRAINAGE &&
                phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer)
                mode = ModeIndexType::IMBIBITION;

            // Update mode based on flow direction if we're in a pure state
            // For non-wetting phase: if saturation is increasing (imbibition), mode should be IMBIBITION
            // If saturation is decreasing (drainage), mode should be DRAINAGE
            if (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::IMBIBITION) {
                // If current saturation is significantly above min historical, we're in imbibition
                if (phaseVolFraction[ipNonWetting] > phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer) {
                    if (mode == ModeIndexType::DRAINAGE) {
                        mode = ModeIndexType::DRAINAGE_TO_IMBIBITION;
                    }
                }
                // If current saturation is significantly below max historical, we're in drainage
                else if (phaseVolFraction[ipNonWetting] < phaseMaxHistoricalVolFraction[ipNonWetting] - flowReversalBuffer) {
                    if (mode == ModeIndexType::IMBIBITION) {
                        mode = ModeIndexType::IMBIBITION_TO_DRAINAGE;
                    }
                }
            }

            // Use simple drainage/imbibition curves when:
            // 1. No hysteresis enabled, OR
            // 2. We're in pure DRAINAGE mode (use drainage curve), OR
            // 3. We're in pure IMBIBITION mode (use imbibition curve)
            // Use scanning curves only during transitions (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
            if (!m_phaseHasHysteresis[TTP::INTERMEDIATE_NONWETTING] ||
                mode == ModeIndexType::DRAINAGE ||
                mode == ModeIndexType::IMBIBITION) {
                phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(phaseVolFraction[ipNonWetting],
                                                                       (mode == ModeIndexType::DRAINAGE)
                                                                       ? m_nonWettingCurve.drainageExtremaPhaseVolFraction
                                                                       : m_nonWettingCurve.imbibitionExtremaPhaseVolFraction);
                // Ensure mode is a valid array index (0 or 1)
                integer const arrayIndex = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE : ModeIndexType::IMBIBITION;
                GEOS_ASSERT_MSG(arrayIndex >= 0 && arrayIndex < static_cast<integer>(m_wettingNonWettingCapillaryPressureKernelWrappers.size()),
                                "Invalid array index for kernel wrapper access");
                computeBoundCapillaryPressure(
                        m_wettingNonWettingCapillaryPressureKernelWrappers[arrayIndex],
                        phaseVolFraction[ipNonWetting],
                        phaseCapPressure[ipNonWetting],
                        dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting]);
                // when pc is on the gas phase, we need to multiply user input by -1
                // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
                phaseCapPressure[ipNonWetting] *= -1;
                dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;

            } else {
                // We're in a transition state (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
                // Use scanning curves

                computeImbibitionNonWettingCapillaryPressure(m_wettingNonWettingCapillaryPressureKernelWrappers,
                                                             m_nonWettingCurve,
                                                             m_landParam[ipNonWetting],
                                                             phaseVolFraction[ipNonWetting],
                                                             phaseMaxHistoricalVolFraction[ipNonWetting],
                                                             phaseTrappedVolFrac[ipNonWetting],
                                                             phaseCapPressure[ipNonWetting],
                                                             dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting],
                                                             mode);

                // when pc is on the gas phase, we need to multiply user input by -1
                // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
                phaseCapPressure[ipNonWetting] *= -1;
                dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
            }

// trapped vol fraction
            if (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {

                {
                    real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
                    real64 const Shy = (phaseMaxHistoricalVolFraction[ipNonWetting] > Smax)
                                       ? phaseMaxHistoricalVolFraction[ipNonWetting]
                                       : Smax;
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_nonWettingCurve, Shy,
                                                                               m_landParam[ipNonWetting],
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);
                    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipNonWetting]);

                    //keep the same Land coeff as two phase only
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_nonWettingCurve.toWetting(), Shy,
                                                                               m_landParam[ipNonWetting],
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);
                    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipNonWetting]);
                }
            } else if (mode == ModeIndexType::IMBIBITION || mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {

                {
                    real64 const Smin = m_nonWettingCurve.imbibitionExtremaPhaseVolFraction;;
                    real64 const Shy = (phaseMinHistoricalVolFraction[ipNonWetting] > Smin)
                                       ? phaseMinHistoricalVolFraction[ipNonWetting]
                                       : Smin;
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_nonWettingCurve, Shy,
                                                                               m_landParam[ipNonWetting],
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);
                    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipNonWetting]);

                    //keep the same Land coeff as two phase only
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_nonWettingCurve.toWetting(), Shy,
                                                                               m_landParam[ipNonWetting],
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);
                    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipNonWetting]);
                }
            }


        }

        void TableCapillaryPressureHysteresis::KernelWrapper::computeThreePhase(const geos::integer ipWetting,
                                                                                const geos::integer GEOS_UNUSED_PARAM( ipInter ),
                                                                                const geos::integer ipNonWetting,
                                                                                const arraySlice1d<const geos::real64,
                                                                                        compflow::USD_PHASE -
                                                                                        1> &phaseVolFraction,
                                                                                const arraySlice1d<const geos::real64,
                                                                                        compflow::USD_PHASE
                                                                                        -
                                                                                        1> &phaseMaxHistoricalVolFraction,
                                                                                const arraySlice1d<const geos::real64,
                                                                                        compflow::USD_PHASE
                                                                                        -
                                                                                        1> &phaseMinHistoricalVolFraction,
                                                                                const arraySlice1d<geos::real64,
                                                                                        relperm::USD_RELPERM
                                                                                        - 2> &phaseTrappedVolFrac,
                                                                                const arraySlice1d<geos::real64,
                                                                                        relperm::USD_RELPERM -
                                                                                        2> &phaseCapPressure,
                                                                                const arraySlice2d<geos::real64,
                                                                                        relperm::USD_RELPERM_DS
                                                                                        -
                                                                                        2> &dPhaseCapPressure_dPhaseVolFrac,
                                                                                ModeIndexType &mode) const {


            LvArray::forValuesInSlice(dPhaseCapPressure_dPhaseVolFrac, [](real64 &val) { val = 0.0; });
            using TTP = ThreePhasePairPhaseType;

            // -- wetting curve if drainage only
            if (!m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] ||
                (mode == ModeIndexType::DRAINAGE &&
                 phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer) ||
                (mode == ModeIndexType::IMBIBITION &&
                 phaseVolFraction[ipWetting] >= phaseMaxHistoricalVolFraction[ipWetting] + flowReversalBuffer)) {
                // water-oil capillary pressure
                phaseTrappedVolFrac[ipWetting] = LvArray::math::min(phaseVolFraction[ipWetting],
                                                                    m_wettingCurve.oppositeBoundPhaseVolFraction);
                phaseCapPressure[ipWetting] =
                        m_wettingIntermediateCapillaryPressureKernelWrappers[mode].compute(
                                &(phaseVolFraction)[ipWetting],
                                &(dPhaseCapPressure_dPhaseVolFrac)[ipWetting][ipWetting]);
            } else {
                mode = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE_TO_IMBIBITION
                                                         : ModeIndexType::IMBIBITION_TO_DRAINAGE;
                computeImbibitionWettingCapillaryPressure(m_wettingIntermediateCapillaryPressureKernelWrappers,
                                                          m_wettingCurve,
                                                          m_nonWettingCurve,
                                                          m_landParam[ipWetting],
                                                          phaseVolFraction[ipWetting],
                                                          phaseMinHistoricalVolFraction[ipWetting],
                                                          phaseTrappedVolFrac[ipWetting],
                                                          phaseCapPressure[ipWetting],
                                                          dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting],
                                                          mode);


            }


            // -- non-wetting cure if drainage only
            // gas-oil capillary pressure
            if (!m_phaseHasHysteresis[TTP::INTERMEDIATE_NONWETTING] ||
                (mode == ModeIndexType::DRAINAGE &&
                 phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] + flowReversalBuffer) ||
                (mode == ModeIndexType::IMBIBITION &&
                 phaseVolFraction[ipNonWetting] <= phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer)) {
                phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(phaseVolFraction[ipNonWetting],
                                                                       (mode == ModeIndexType::DRAINAGE)
                                                                       ? m_nonWettingCurve.drainageExtremaPhaseVolFraction
                                                                       : m_nonWettingCurve.imbibitionExtremaPhaseVolFraction);
                phaseCapPressure[ipNonWetting] =
                        m_nonWettingIntermediateCapillaryPressureKernelWrappers[mode].compute(
                                &(phaseVolFraction)[ipNonWetting],
                                &(dPhaseCapPressure_dPhaseVolFrac)[ipNonWetting][ipNonWetting]);


                // when pc is on the gas phase, we need to multiply user input by -1
                // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
                phaseCapPressure[ipNonWetting] *= -1;
                dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
            } else {

                mode = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE_TO_IMBIBITION
                                                         : ModeIndexType::IMBIBITION_TO_DRAINAGE;

                computeImbibitionNonWettingCapillaryPressure(m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                                                             m_nonWettingCurve,
                                                             m_wettingCurve,
                                                             m_landParam[ipNonWetting],
                                                             phaseVolFraction[ipNonWetting],
                                                             phaseMinHistoricalVolFraction[ipNonWetting],
                                                             phaseTrappedVolFrac[ipNonWetting],
                                                             phaseCapPressure[ipNonWetting],
                                                             dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting],
                                                             mode);
                // when pc is on the gas phase, we need to multiply user input by -1
                // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
                phaseCapPressure[ipNonWetting] *= -1;
                dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;

                //update trapped fraction
                if (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {


                    {
                        real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
                        real64 const Shy = (phaseMinHistoricalVolFraction[ipWetting] < Smin)
                                           ? phaseMinHistoricalVolFraction[ipWetting]
                                           : Smin;
                        real64 Scrt = 0.0;
                        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_wettingCurve, Shy,
                                                                                   m_landParam[ipWetting],
                                                                                   m_jerauldParam_a, m_jerauldParam_b,
                                                                                   Scrt);
                        phaseTrappedVolFrac[ipWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipWetting]);
                    }

                    {
                        real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
                        real64 const Shy = (phaseMaxHistoricalVolFraction[ipNonWetting] > Smax)
                                           ? phaseMaxHistoricalVolFraction[ipNonWetting]
                                           : Smax;
                        real64 Scrt = 0.0;
                        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_nonWettingCurve, Shy,
                                                                                   m_landParam[ipNonWetting],
                                                                                   m_jerauldParam_a, m_jerauldParam_b,
                                                                                   Scrt);
                        phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipNonWetting]);
                    }
                } else if (mode == ModeIndexType::IMBIBITION || mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                    {
                        real64 const Smax = m_wettingCurve.imbibitionExtremaPhaseVolFraction;
                        real64 const Shy = (phaseMaxHistoricalVolFraction[ipWetting] < Smax)
                                           ? phaseMaxHistoricalVolFraction[ipWetting]
                                           : Smax;
                        real64 Scrt = 0.0;
                        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_wettingCurve, Shy,
                                                                                   m_landParam[ipWetting],
                                                                                   m_jerauldParam_a, m_jerauldParam_b,
                                                                                   Scrt);
                        phaseTrappedVolFrac[ipWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipWetting]);
                    }

                    {
                        real64 const Smin = m_nonWettingCurve.imbibitionExtremaPhaseVolFraction;;
                        real64 const Shy = (phaseMinHistoricalVolFraction[ipNonWetting] > Smin)
                                           ? phaseMinHistoricalVolFraction[ipNonWetting]
                                           : Smin;
                        real64 Scrt = 0.0;
                        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(m_nonWettingCurve, Shy,
                                                                                   m_landParam[ipNonWetting],
                                                                                   m_jerauldParam_a, m_jerauldParam_b,
                                                                                   Scrt);
                        phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min(Scrt, phaseVolFraction[ipNonWetting]);
                    }
                }


            }


        }


        void
        TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionNonWettingCapillaryPressure(
                const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
                const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                const KilloughHysteresis::HysteresisCurve &wettingCurve,
                const geos::real64 &landParam,
                const geos::real64 &phaseVolFraction,
                const geos::real64 &phaseMaxHistoricalVolFraction,
                geos::real64 &phaseTrappedVolFrac,
                geos::real64 &phaseCapPressure,
                geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                const ModeIndexType &mode) const {
            GEOS_ASSERT(!nonWettingCurve.isWetting());
            real64 const S = phaseVolFraction;
            real64 const Smii = nonWettingCurve.imbibitionExtremaPhaseVolFraction;
            real64 const Smid = nonWettingCurve.drainageExtremaPhaseVolFraction;
            real64 const Smax = nonWettingCurve.oppositeBoundPhaseVolFraction;

	    GEOS_UNUSED_VAR( Smii, Smid );
	    
//  if( S >= Smax )
//  {
//    //above accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else if( S <= Smid )
//  {
//    //below accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else
            {
                //drainage to imbibition
                real64 dpci_dS, dpcd_dS;
                real64 const pci = nonWettingKernelWrapper[ModeIndexType::IMBIBITION].compute(&S, &dpci_dS);
                real64 const pcd = nonWettingKernelWrapper[ModeIndexType::DRAINAGE].compute(&S, &dpcd_dS);

                // Step 1: get the trapped from wetting data
                real64 const Shy = (phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax;

                //drainage to imbibition
                if (mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(nonWettingCurve, Shy, landParam,
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);

                    real64 const E = m_killoughCurvatureParamCapPres;

                    //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
                    real64 F = (1. / (Shy - S + E) - 1. / E) / (1. / (Shy - Scrt + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. compute dF_dS
                    real64 dF_dS = (1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

                    //Step 4. Eventually assemble everything following (EQ. 34.20)
                    phaseCapPressure = pcd + F * (pci - pcd);
                    dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpci_dS - dpcd_dS);
                    dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pci - pcd);

                    //update trapped fraction
                    phaseTrappedVolFrac = LvArray::math::min(Scrt, S);

                }
                    //imbibition to drainage
                else if (mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(wettingCurve, Shy, landParam,
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);
                    real64 Sgma = 1. - Scrt - m_phaseIntermediateMinVolFraction;

                    real64 const E = m_killoughCurvatureParamCapPres;

                    //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
                    real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Sgma - Shy + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. compute dF_dS
                    real64 dF_dS = (-1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

                    //Step 4. Eventually assemble everything following (EQ. 34.20)
                    phaseCapPressure = pci + F * (pcd - pci);
                    dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpcd_dS - dpci_dS);
                    dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pcd - pci);
                } else {
                    GEOS_THROW(GEOS_FMT("{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                                          "TableCapillaryPressureHysteresis",
                                          (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                                             ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                                InputError);
                }


            }
        }


        void
        TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionNonWettingCapillaryPressure(
                const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
                const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                const geos::real64 &landParam,
                const geos::real64 &phaseVolFraction,
                const geos::real64 &phaseMaxHistoricalVolFraction,
                geos::real64 &phaseTrappedVolFrac,
                geos::real64 &phaseCapPressure,
                geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                const ModeIndexType &mode) const {

            GEOS_ASSERT(!nonWettingCurve.isWetting());
            real64 const S = phaseVolFraction;
            real64 const Smii = nonWettingCurve.imbibitionExtremaPhaseVolFraction;
            real64 const Smid = nonWettingCurve.drainageExtremaPhaseVolFraction;
            real64 const Smax = nonWettingCurve.oppositeBoundPhaseVolFraction;

	    GEOS_UNUSED_VAR( Smii, Smid );
	    
//  if( S >= Smax )
//  {
//    //above accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else if( S <= Smid )
//  {
//    //below accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else
            {
                //drainage to imbibition
                real64 dpci_dS, dpcd_dS;
                real64 const pci = nonWettingKernelWrapper[ModeIndexType::IMBIBITION].compute(&S, &dpci_dS);
                real64 const pcd = nonWettingKernelWrapper[ModeIndexType::DRAINAGE].compute(&S, &dpcd_dS);

                // Step 1: get the trapped from wetting data
                real64 const Shy = (phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax;

                //drainage to imbibition
                if (mode == ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(nonWettingCurve, Shy, landParam,
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);

                    real64 const E = m_killoughCurvatureParamCapPres;

                    //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
                    real64 F = (1. / (Shy - S + E) - 1. / E) / (1. / (Shy - Scrt + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. compute dF_dS
                    real64 dF_dS = (1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

                    //Step 4. Eventually assemble everything following (EQ. 34.20)
                    phaseCapPressure = pcd + F * (pci - pcd);
                    dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpci_dS - dpcd_dS);
                    dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pci - pcd);

                    //update trapped fraction
                    phaseTrappedVolFrac = LvArray::math::min(Scrt, S);

                }
                    //imbibition to drainage
                else if (mode == ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                    real64 Scrt = 0.0;
                    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(nonWettingCurve, Shy, landParam,
                                                                               m_jerauldParam_a, m_jerauldParam_b,
                                                                               Scrt);
                    real64 Sgma = 1. - (1. - Scrt);

                    real64 const E = m_killoughCurvatureParamCapPres;

                    //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
                    real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Sgma - Shy + E) - 1. / E);
                    //force bound
                    F = LvArray::math::max(F, 0.0);
                    F = LvArray::math::min(F, 1.0);

                    //Step 3. compute dF_dS
                    real64 dF_dS = (-1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

                    //Step 4. Eventually assemble everything following (EQ. 34.20)
                    phaseCapPressure = pci + F * (pcd - pci);
                    dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpcd_dS - dpci_dS);
                    dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pcd - pci);
                } else {
                    GEOS_THROW(GEOS_FMT("{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                                          "TableCapillaryPressureHysteresis",
                                          (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                                             ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                                InputError);
                }


            }
        }


        void TableCapillaryPressureHysteresis::KernelWrapper::computeBoundCapillaryPressure(
                const TableFunction::KernelWrapper &drainageRelpermWrapper,
                const geos::real64 &phaseVolFraction,
                geos::real64 &phaseCapPressure,
                geos::real64 &dPhaseCapPressure_dPhaseVolFrac) const {
            phaseCapPressure = drainageRelpermWrapper.compute(&phaseVolFraction,
                                                              &dPhaseCapPressure_dPhaseVolFrac);

        }

        /// Helper function to compute Pc(S) for given S, mode, and historical values
        /// Used for Newton-Raphson inversion in computeInv
        GEOS_HOST_DEVICE
        inline real64 TableCapillaryPressureHysteresis::KernelWrapper::computeCapillaryPressureForSaturation(
                real64 const S,
                fields::cappres::ModeIndexType const &mode,
                integer const ipPhase,
                real64 const &phaseMinHistoricalVolFraction,
                real64 const &phaseMaxHistoricalVolFraction,
                arrayView1d<TableFunction::KernelWrapper const> const &capPresKernelWrappers,
                KilloughHysteresis::HysteresisCurve const &wettingCurve,
                KilloughHysteresis::HysteresisCurve const &nonWettingCurve,
                real64 const &landParam,
                real64 const &phaseIntermediateMinVolFraction,
                real64 const &killoughCurvatureParam,
                real64 const &jerauldParam_a,
                real64 const &jerauldParam_b,
                bool const isWettingPhase,
                real64 const precomputedScrt,
                real64 const precomputedDenomF,
                real64 const precomputedShy) const
        {
            real64 pc = 0.0;
            real64 dpc_dS = 0.0;
            
            // For pure drainage or imbibition modes, use the table directly
            if (mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION) {
                integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? fields::cappres::ModeIndexType::DRAINAGE : fields::cappres::ModeIndexType::IMBIBITION;
                if (arrayIndex < static_cast<integer>(capPresKernelWrappers.size())) {
                    pc = capPresKernelWrappers[arrayIndex].compute(&S, &dpc_dS);
                }
            }
            // For scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
            else {
                // Precomputed values should be provided by caller (from computeFlux before local_solver)
                // If not provided, fall back to drainage curve to avoid calling computeTrappedCriticalPhaseVolFraction
                if (precomputedScrt < 0.0) {
                    // Precomputed values not provided - return drainage curve value
                    real64 const arrayIndex = fields::cappres::ModeIndexType::DRAINAGE;
                    if (arrayIndex < static_cast<integer>(capPresKernelWrappers.size())) {
                        pc = capPresKernelWrappers[arrayIndex].compute(&S, &dpc_dS);
                    }
                    return pc;
                }
                
                real64 dpci_dS, dpcd_dS;
                real64 const pci = capPresKernelWrappers[fields::cappres::ModeIndexType::IMBIBITION].compute(&S, &dpci_dS);
                real64 const pcd = capPresKernelWrappers[fields::cappres::ModeIndexType::DRAINAGE].compute(&S, &dpcd_dS);
                
                real64 const E = killoughCurvatureParam;
                
                if (isWettingPhase) {
                    real64 const Smin = wettingCurve.oppositeBoundPhaseVolFraction;
                    real64 const Shy = (precomputedShy >= 0.0) ? precomputedShy : 
                                      ((phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin);
                    
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                        // Use precomputed value - should have been computed before local_solver
                        real64 const Scrt = precomputedScrt;
                        real64 const Swma = 1 - Scrt - phaseIntermediateMinVolFraction;
                        real64 denomF = precomputedDenomF;
                        if (LvArray::math::abs(precomputedDenomF) < 1e-15) {
                            denomF = (1. / (Swma - Shy + E) - 1. / E);
                        }
                        // Guard against division by zero
                        real64 F = 0.0;
                        if (LvArray::math::abs(denomF) >= 1e-15) {
                            real64 const F_num = (1. / (S - Shy + E) - 1. / E);
                            F = F_num / denomF;
                        }
                        F = LvArray::math::max(F, 0.0);
                        F = LvArray::math::min(F, 1.0);
                        
                        pc = pcd + F * (pci - pcd);
                        dpc_dS = dpcd_dS + F * (dpci_dS - dpcd_dS);
                    }
                    else if (mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                        // Use precomputed value - should have been computed before local_solver
                        real64 const Scrt = precomputedScrt;
                        
                        real64 denomF = precomputedDenomF;
                        if (LvArray::math::abs(precomputedDenomF) < 1e-15) {
                            denomF = (1. / (Shy - Scrt + E) - 1. / E);
                        }
                        // Guard against division by zero
                        real64 F = 0.0;
                        if (LvArray::math::abs(denomF) >= 1e-15) {
                            real64 const F_num = (1. / (Shy - S + E) - 1. / E);
                            F = F_num / denomF;
                        }
                        F = LvArray::math::max(F, 0.0);
                        F = LvArray::math::min(F, 1.0);
                        
                        pc = pci + F * (pcd - pci);
                        dpc_dS = dpci_dS + F * (dpcd_dS - dpci_dS);
                    }
                }
                else {
                    // Non-wetting phase
                    real64 const Smax = nonWettingCurve.oppositeBoundPhaseVolFraction;
                    real64 const Shy = (precomputedShy >= 0.0) ? precomputedShy :
                                      ((phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax);
                    
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                        // Use precomputed value - should have been computed before local_solver
                        real64 const Scrt = precomputedScrt;
                        
                        real64 denomF = precomputedDenomF;
                        if (LvArray::math::abs(precomputedDenomF) < 1e-15) {
                            denomF = (1. / (Shy - Scrt + E) - 1. / E);
                        }
                        // Guard against division by zero
                        real64 F = 0.0;
                        if (LvArray::math::abs(denomF) >= 1e-15) {
                            real64 const F_num = (1. / (Shy - S + E) - 1. / E);
                            F = F_num / denomF;
                        }
                        F = LvArray::math::max(F, 0.0);
                        F = LvArray::math::min(F, 1.0);
                        
                        pc = pcd + F * (pci - pcd);
                        dpc_dS = dpci_dS + F * (dpci_dS - dpcd_dS);
                        // Note: there's also a dF/dS term in the non-wetting case, but for inversion we approximate
                    }
                    else if (mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                        // Use precomputed value - should have been computed before local_solver
                        real64 const Scrt = precomputedScrt;
                        real64 const Sgma = 1. - Scrt - phaseIntermediateMinVolFraction;
                        
                        real64 denomF = precomputedDenomF;
                        if (LvArray::math::abs(precomputedDenomF) < 1e-15) {
                            denomF = (1. / (Sgma - Shy + E) - 1. / E);
                        }
                        // Guard against division by zero
                        real64 F = 0.0;
                        if (LvArray::math::abs(denomF) >= 1e-15) {
                            real64 const F_num = (1. / (S - Shy + E) - 1. / E);
                            F = F_num / denomF;
                        }
                        F = LvArray::math::max(F, 0.0);
                        F = LvArray::math::min(F, 1.0);
                        
                        pc = pci + F * (pcd - pci);
                        dpc_dS = dpcd_dS + F * (dpcd_dS - dpci_dS);
                        // Note: there's also a dF/dS term in the non-wetting case, but for inversion we approximate
                    }
                }
            }
            
            return pc;
        }

        void
        TableCapillaryPressureHysteresis::KernelWrapper::computeInv(
                arraySlice1d<real64, compflow::USD_PHASE - 1> const &phaseVolFraction,
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                arraySlice1d<real64 const, cappres::USD_CAPPRES - 2> const &phaseTrappedVolFrac,
                arraySlice1d<real64 const, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                arraySlice2d<real64, cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                fields::cappres::ModeIndexType const &mode) const
        {
            LvArray::forValuesInSlice(dPhaseCapPressure_dPhaseVolFrac, [](real64 &val) { val = 0.0; });

            using PT = CapillaryPressureBase::PhaseType;
            integer const ipWater = (PT::WATER < m_phaseOrder.size()) ? m_phaseOrder[PT::WATER] : -1;
            integer const ipOil = (PT::OIL < m_phaseOrder.size()) ? m_phaseOrder[PT::OIL] : -1;
            integer const ipGas = (PT::GAS < m_phaseOrder.size()) ? m_phaseOrder[PT::GAS] : -1;

            // Newton-Raphson parameters
            constexpr real64 tol = 1e-9;
            constexpr integer maxIter = 20;
            constexpr real64 minS = 0.0;
            constexpr real64 maxS = 1.0;

            // Determine which phase pairs need inversion
            if (ipWater >= 0 && ipOil >= 0 && ipGas >= 0) {
                // Three-phase: invert wetting and non-wetting capillary pressures
                
                // 1. Invert wetting phase (water-oil capillary pressure)
                constexpr real64 pcEpsilon = 1e-10;
                if (ipWater >= 0 && LvArray::math::abs(phaseCapPressure[ipWater]) > pcEpsilon) {
                    real64 S_guess = phaseMaxHistoricalVolFraction[ipWater];
                    if (S_guess <= phaseMinHistoricalVolFraction[ipWater] || S_guess < minS || S_guess > maxS) {
                        // Fall back to drainage/imbibition curve evaluation if available
                        real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
                        real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
                        S_guess = 0.5 * (Smin + Smax);
                        S_guess = LvArray::math::max(S_guess, minS);
                        S_guess = LvArray::math::min(S_guess, maxS);
                    }
                    
                    // Precomputed values should be provided by caller (from computeFlux before local_solver)
                    // If not provided (sentinel values), fall back to simple drainage/imbibition curves
                    // This avoids calling computeTrappedCriticalPhaseVolFraction inside local_solver
                    real64 precomputedScrt_water = -1.0;
                    real64 precomputedDenomF_water = 0.0;
                    real64 precomputedShy_water = -1.0;
                    
                    // For scanning curves without precomputed values, use drainage curve as fallback
                    // This avoids the floating-point exception that would occur in computeTrappedCriticalPhaseVolFraction
                    
                    real64 S = S_guess;
                    for (integer iter = 0; iter < maxIter; ++iter) {
                        // Compute Pc(S) using the helper function
                        real64 pc_computed = this->computeCapillaryPressureForSaturation(
                            S, mode, ipWater,
                            phaseMinHistoricalVolFraction[ipWater],
                            phaseMaxHistoricalVolFraction[ipWater],
                            m_wettingIntermediateCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipWater],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            true, // isWettingPhase
                            precomputedScrt_water,
                            precomputedDenomF_water,
                            precomputedShy_water);
                        
                        real64 residual = pc_computed - phaseCapPressure[ipWater];
                        
                        if (LvArray::math::abs(residual) < tol) {
                            break;
                        }
                        
                        // Compute derivative dPc/dS (approximate using finite difference if needed)
                        // For simplicity, use the derivative from the table evaluation
                        // In a more sophisticated implementation, we could compute the analytical derivative
                        real64 const dS = 1e-6;
                        real64 const S_pert = LvArray::math::min(S + dS, maxS);
                        real64 const pc_pert = this->computeCapillaryPressureForSaturation(
                            S_pert, mode, ipWater,
                            phaseMinHistoricalVolFraction[ipWater],
                            phaseMaxHistoricalVolFraction[ipWater],
                            m_wettingIntermediateCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipWater],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            true, // isWettingPhase
                            precomputedScrt_water,
                            precomputedDenomF_water,
                            precomputedShy_water);
                        
                        real64 const dpc_dS = (pc_pert - pc_computed) / dS;
                        
                        // Avoid division by zero
                        if (LvArray::math::abs(dpc_dS) > 1e-12) {
                            S = S - residual / dpc_dS;
                            // Clamp to valid range
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        } else {
                            // If derivative is zero, try bisection or use a small step
                            S = S - residual * 1e-6; // Small fixed step
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        }
                    }
                    
                    phaseVolFraction[ipWater] = S;
                    phaseVolFraction[ipOil] = 1.0 - S - phaseVolFraction[ipGas];
                    
                    // Compute derivative at final S
                    real64 dpc_dS_final = 0.0;
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION) {
                        integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? fields::cappres::ModeIndexType::DRAINAGE : fields::cappres::ModeIndexType::IMBIBITION;
                        m_wettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute(&S, &dpc_dS_final);
                    }
                    dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] = dpc_dS_final;
                }
                
                // 2. Invert non-wetting phase (gas-oil capillary pressure)
                if (ipGas >= 0 && LvArray::math::abs(phaseCapPressure[ipGas]) > pcEpsilon) {
                    // Note: For non-wetting phase, the capillary pressure is typically negative
                    // and the relationship may be inverted (increasing Pc with decreasing S)
                    real64 S_guess = phaseMaxHistoricalVolFraction[ipGas];
                    if (S_guess <= phaseMinHistoricalVolFraction[ipGas] || S_guess < minS || S_guess > maxS) {
                        real64 const Smin = m_nonWettingCurve.imbibitionExtremaPhaseVolFraction;
                        real64 const Smax = m_nonWettingCurve.drainageExtremaPhaseVolFraction;
                        S_guess = 0.5 * (Smin + Smax);
                        S_guess = LvArray::math::max(S_guess, minS);
                        S_guess = LvArray::math::min(S_guess, maxS);
                    }
                    
                    // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
                    real64 precomputedScrt_gas = -1.0;
                    real64 precomputedDenomF_gas = 0.0;
                    real64 precomputedShy_gas = -1.0;
                    
                    // Only precompute for scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION || 
                        mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                        real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
                        precomputedShy_gas = (phaseMaxHistoricalVolFraction[ipGas] < Smax) ? 
                                            phaseMaxHistoricalVolFraction[ipGas] : Smax;
                        real64 const E = m_killoughCurvatureParamCapPres;
                        
                        if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                            // Precomputed values should be provided by caller - skip computation here
                            // to avoid calling computeTrappedCriticalPhaseVolFraction inside local_solver
                            precomputedScrt_gas = -1.0; // Indicates not precomputed
                            precomputedDenomF_gas = 0.0;
                        }
                        else if (mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                            // Precomputed values should be provided by caller - skip computation here
                            // to avoid calling computeTrappedCriticalPhaseVolFraction inside local_solver
                            precomputedScrt_gas = -1.0; // Indicates not precomputed
                            precomputedDenomF_gas = 0.0;
                        }
                    }
                    
                    real64 S = S_guess;
                    for (integer iter = 0; iter < maxIter; ++iter) {
                        real64 pc_computed = this->computeCapillaryPressureForSaturation(
                            S, mode, ipGas,
                            phaseMinHistoricalVolFraction[ipGas],
                            phaseMaxHistoricalVolFraction[ipGas],
                            m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipGas],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            false, // isWettingPhase = false
                            precomputedScrt_gas,
                            precomputedDenomF_gas,
                            precomputedShy_gas);
                        
                        // For non-wetting phase, capillary pressure is typically multiplied by -1
                        // Check which sign to use based on the input
                        real64 residual = pc_computed - phaseCapPressure[ipGas];
                        
                        if (LvArray::math::abs(residual) < tol) {
                            break;
                        }
                        
                        real64 const dS = 1e-6;
                        real64 const S_pert = LvArray::math::min(S + dS, maxS);
                        real64 const pc_pert = this->computeCapillaryPressureForSaturation(
                            S_pert, mode, ipGas,
                            phaseMinHistoricalVolFraction[ipGas],
                            phaseMaxHistoricalVolFraction[ipGas],
                            m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipGas],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            false, // isWettingPhase = false
                            precomputedScrt_gas,
                            precomputedDenomF_gas,
                            precomputedShy_gas);
                        
                        real64 const dpc_dS = (pc_pert - pc_computed) / dS;
                        
                        if (LvArray::math::abs(dpc_dS) > 1e-12) {
                            S = S - residual / dpc_dS;
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        } else {
                            S = S - residual * 1e-6;
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        }
                    }
                    
                    phaseVolFraction[ipGas] = S;
                    
                    // Update oil phase (ensure sum = 1)
                    if (ipWater >= 0 && ipOil >= 0) {
                        phaseVolFraction[ipOil] = 1.0 - phaseVolFraction[ipWater] - phaseVolFraction[ipGas];
                        phaseVolFraction[ipOil] = LvArray::math::max(phaseVolFraction[ipOil], 0.0);
                        phaseVolFraction[ipOil] = LvArray::math::min(phaseVolFraction[ipOil], 1.0);
                    }
                    
                    // Compute derivative at final S
                    real64 dpc_dS_final = 0.0;
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION) {
                        integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? fields::cappres::ModeIndexType::DRAINAGE : fields::cappres::ModeIndexType::IMBIBITION;
                        m_nonWettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute(&S, &dpc_dS_final);
                    }
                    dPhaseCapPressure_dPhaseVolFrac[ipGas][ipGas] = dpc_dS_final;
                }
            }
            else if (ipWater < 0) {
                // Two-phase: oil-gas (non-wetting phase)
                // Similar to above but simpler
                constexpr real64 pcEpsilon = 1e-10;
                if (ipGas >= 0 && LvArray::math::abs(phaseCapPressure[ipGas]) > pcEpsilon) {
                    real64 S_guess = phaseMaxHistoricalVolFraction[ipGas];
                    if (S_guess < minS || S_guess > maxS) {
                        S_guess = 0.5;
                    }
                    
                    // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
                    real64 precomputedScrt_gas = -1.0;
                    real64 precomputedDenomF_gas = 0.0;
                    real64 precomputedShy_gas = -1.0;
                    
                    // Only precompute for scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION || 
                        mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                        real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
                        precomputedShy_gas = (phaseMaxHistoricalVolFraction[ipGas] < Smax) ? 
                                            phaseMaxHistoricalVolFraction[ipGas] : Smax;
                        real64 const E = m_killoughCurvatureParamCapPres;
                        
                        if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                            // Precomputed values should be provided by caller - skip computation here
                            // to avoid calling computeTrappedCriticalPhaseVolFraction inside local_solver
                            precomputedScrt_gas = -1.0; // Indicates not precomputed
                            precomputedDenomF_gas = 0.0;
                        }
                        else if (mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                            // Precomputed values should be provided by caller - skip computation here
                            // to avoid calling computeTrappedCriticalPhaseVolFraction inside local_solver
                            precomputedScrt_gas = -1.0; // Indicates not precomputed
                            precomputedDenomF_gas = 0.0;
                        }
                    }
                    
                    real64 S = S_guess;
                    for (integer iter = 0; iter < maxIter; ++iter) {
                        real64 pc_computed = this->computeCapillaryPressureForSaturation(
                            S, mode, ipGas,
                            phaseMinHistoricalVolFraction[ipGas],
                            phaseMaxHistoricalVolFraction[ipGas],
                            m_wettingNonWettingCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipGas],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            false, // isWettingPhase = false
                            precomputedScrt_gas,
                            precomputedDenomF_gas,
                            precomputedShy_gas);
                        
                        real64 residual = pc_computed - phaseCapPressure[ipGas];
                        
                        if (LvArray::math::abs(residual) < tol) {
                            break;
                        }
                        
                        real64 const dS = 1e-6;
                        real64 const S_pert = LvArray::math::min(S + dS, maxS);
                        real64 const pc_pert = this->computeCapillaryPressureForSaturation(
                            S_pert, mode, ipGas,
                            phaseMinHistoricalVolFraction[ipGas],
                            phaseMaxHistoricalVolFraction[ipGas],
                            m_wettingNonWettingCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipGas],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            false, // isWettingPhase = false
                            precomputedScrt_gas,
                            precomputedDenomF_gas,
                            precomputedShy_gas);
                        
                        real64 const dpc_dS = (pc_pert - pc_computed) / dS;
                        
                        if (LvArray::math::abs(dpc_dS) > 1e-12) {
                            S = S - residual / dpc_dS;
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        } else {
                            S = S - residual * 1e-6;
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        }
                    }
                    
                    phaseVolFraction[ipGas] = S;
                    if (ipOil >= 0) {
                        phaseVolFraction[ipOil] = 1.0 - S;
                    }
                }
            }
            else {
                // Two-phase: water-oil or water-gas (wetting phase)
                constexpr real64 pcEpsilon = 1e-10;
                if (ipWater >= 0 && LvArray::math::abs(phaseCapPressure[ipWater]) > pcEpsilon) {
                    real64 S_guess = phaseMaxHistoricalVolFraction[ipWater];
                    if (S_guess < minS || S_guess > maxS) {
                        real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
                        real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
                        S_guess = 0.5 * (Smin + Smax);
                        S_guess = LvArray::math::max(S_guess, minS);
                        S_guess = LvArray::math::min(S_guess, maxS);
                    }
                    
                    // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
                    real64 precomputedScrt_water = -1.0;
                    real64 precomputedDenomF_water = 0.0;
                    real64 precomputedShy_water = -1.0;
                    
                    // Only precompute for scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
                    if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION || 
                        mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                        real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
                        precomputedShy_water = (phaseMinHistoricalVolFraction[ipWater] > Smin) ? 
                                              phaseMinHistoricalVolFraction[ipWater] : Smin;
                        real64 const E = m_killoughCurvatureParamCapPres;
                        
                        if (mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION) {
                            // Precomputed values should be provided by caller - skip computation here
                            // to avoid calling computeTrappedCriticalPhaseVolFraction inside local_solver
                            precomputedScrt_water = -1.0; // Indicates not precomputed
                            precomputedDenomF_water = 0.0;
                        }
                        else if (mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE) {
                            // Precomputed values should be provided by caller - skip computation here
                            // to avoid calling computeTrappedCriticalPhaseVolFraction inside local_solver
                            precomputedScrt_water = -1.0; // Indicates not precomputed
                            precomputedDenomF_water = 0.0;
                        }
                    }
                    
                    real64 S = S_guess;
                    for (integer iter = 0; iter < maxIter; ++iter) {
                        real64 pc_computed = this->computeCapillaryPressureForSaturation(
                            S, mode, ipWater,
                            phaseMinHistoricalVolFraction[ipWater],
                            phaseMaxHistoricalVolFraction[ipWater],
                            m_wettingNonWettingCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipWater],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            true, // isWettingPhase
                            precomputedScrt_water,
                            precomputedDenomF_water,
                            precomputedShy_water);
                        
                        real64 residual = pc_computed - phaseCapPressure[ipWater];
                        
                        if (LvArray::math::abs(residual) < tol) {
                            break;
                        }
                        
                        real64 const dS = 1e-6;
                        real64 const S_pert = LvArray::math::min(S + dS, maxS);
                        real64 const pc_pert = this->computeCapillaryPressureForSaturation(
                            S_pert, mode, ipWater,
                            phaseMinHistoricalVolFraction[ipWater],
                            phaseMaxHistoricalVolFraction[ipWater],
                            m_wettingNonWettingCapillaryPressureKernelWrappers,
                            m_wettingCurve,
                            m_nonWettingCurve,
                            m_landParam[ipWater],
                            m_phaseIntermediateMinVolFraction,
                            m_killoughCurvatureParamCapPres,
                            m_jerauldParam_a,
                            m_jerauldParam_b,
                            true, // isWettingPhase
                            precomputedScrt_water,
                            precomputedDenomF_water,
                            precomputedShy_water);
                        
                        real64 const dpc_dS = (pc_pert - pc_computed) / dS;
                        
                        if (LvArray::math::abs(dpc_dS) > 1e-12) {
                            S = S - residual / dpc_dS;
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        } else {
                            S = S - residual * 1e-6;
                            S = LvArray::math::max(S, minS);
                            S = LvArray::math::min(S, maxS);
                        }
                    }
                    
                    phaseVolFraction[ipWater] = S;
                    // Set non-wetting phase
                    if (ipGas >= 0) {
                        phaseVolFraction[ipGas] = 1.0 - S;
                    } else if (ipOil >= 0) {
                        phaseVolFraction[ipOil] = 1.0 - S;
                    }
                }
            }
        }

/// kernel creation

        TableCapillaryPressureHysteresis::KernelWrapper
        TableCapillaryPressureHysteresis::createKernelWrapper() {

            // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
            createAllTableKernelWrappers();

            // Validate that the arrays are properly populated
            integer const numPhases = m_phaseNames.size();
            if (numPhases == 2) {
                GEOS_THROW_IF(m_wettingNonWettingCapillaryPressureKernelWrappers.size() != 2,
                               GEOS_FMT("{}: Expected 2 kernel wrappers for two-phase flow, but got {}. "
                                        "This usually means createAllTableKernelWrappers() failed to populate the arrays. "
                                        "Check that table functions '{}' and '{}' exist and are properly defined.",
                                         getFullName(),
                                         m_wettingNonWettingCapillaryPressureKernelWrappers.size(),
                                         m_drainageWettingNonWettingCapPresTableName,
                                         m_imbibitionWettingNonWettingCapPresTableName.empty() ? m_drainageWettingNonWettingCapPresTableName : m_imbibitionWettingNonWettingCapPresTableName),
                               InputError);
            } else if (numPhases == 3) {
                GEOS_THROW_IF(m_wettingIntermediateCapillaryPressureKernelWrappers.size() != 2,
                               GEOS_FMT("{}: Expected 2 wetting-intermediate kernel wrappers for three-phase flow, but got {}",
                                         getFullName(),
                                         m_wettingIntermediateCapillaryPressureKernelWrappers.size()),
                               InputError);
                GEOS_THROW_IF(m_nonWettingIntermediateCapillaryPressureKernelWrappers.size() != 2,
                               GEOS_FMT("{}: Expected 2 non-wetting-intermediate kernel wrappers for three-phase flow, but got {}",
                                         getFullName(),
                                         m_nonWettingIntermediateCapillaryPressureKernelWrappers.size()),
                               InputError);
            }

            // Validate that m_phaseHasHysteresis is properly initialized
            GEOS_THROW_IF(m_phaseHasHysteresis.size() != 2,
                           GEOS_FMT("{}: m_phaseHasHysteresis must have size 2, but got {}",
                                     getFullName(),
                                     m_phaseHasHysteresis.size()),
                           InputError);

            // Validate that the historical volume fraction arrays have been resized
            // These arrays should be resized by resizeFields() before the KernelWrapper is used
            // If they're empty, it means resizeFields() hasn't been called yet
            GEOS_THROW_IF(m_phaseMaxHistoricalVolFraction.size(0) == 0 || m_phaseMaxHistoricalVolFraction.size(1) == 0,
                           GEOS_FMT("{}: phaseMaxHistoricalVolFraction array has not been resized (size=[{}, {}]). "
                                    "This usually means resizeFields() has not been called yet. "
                                    "The arrays must be resized before the KernelWrapper can be used.",
                                     getFullName(),
                                     m_phaseMaxHistoricalVolFraction.size(0),
                                     m_phaseMaxHistoricalVolFraction.size(1)),
                           InputError);
            GEOS_THROW_IF(m_phaseMinHistoricalVolFraction.size(0) == 0 || m_phaseMinHistoricalVolFraction.size(1) == 0,
                           GEOS_FMT("{}: phaseMinHistoricalVolFraction array has not been resized (size=[{}, {}]). "
                                    "This usually means resizeFields() has not been called yet. "
                                    "The arrays must be resized before the KernelWrapper can be used.",
                                     getFullName(),
                                     m_phaseMinHistoricalVolFraction.size(0),
                                     m_phaseMinHistoricalVolFraction.size(1)),
                           InputError);

            // then we create the actual TableRelativePermeabilityHysteresis::KernelWrapper
            return KernelWrapper(m_wettingNonWettingCapillaryPressureKernelWrappers,
                                 m_wettingIntermediateCapillaryPressureKernelWrappers,
                                 m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                                 m_phaseHasHysteresis,
                                 m_landParam,
                                 m_jerauldParam_a,
                                 m_jerauldParam_b,
                                 m_killoughCurvatureParamCapPres,
                                 m_phaseIntermediateMinVolFraction,
                                 m_wettingCurve,
                                 m_nonWettingCurve,
                                 m_phaseMinHistoricalVolFraction,
                                 m_phaseMaxHistoricalVolFraction,
                                 m_phaseTypes,
                                 m_phaseOrder,
                                 m_mode,
                                 m_phaseTrappedVolFrac,
                                 m_phaseCapPressure,
                                 m_dPhaseCapPressure_dPhaseVolFrac);
        }

        void TableCapillaryPressureHysteresis::createAllTableKernelWrappers() {
            using TPP = ThreePhasePairPhaseType;

            FunctionManager const &functionManager = FunctionManager::getInstance();

            integer const numPhases = m_phaseNames.size();

            // Ensure m_phaseHasHysteresis is initialized before accessing it
            // This can happen if createKernelWrapper is called before postProcessInput
            if (m_phaseHasHysteresis.size() == 0) {
                m_phaseHasHysteresis.resize(2);
                if (numPhases == 2) {
                    m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = 
                        (m_imbibitionWettingNonWettingCapPresTableName.empty() ||
                         m_imbibitionWettingNonWettingCapPresTableName == m_drainageWettingNonWettingCapPresTableName)
                        ? 0 : 1;
                    m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING];
                } else if (numPhases == 3) {
                    m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = 
                        (m_imbibitionWettingIntermediateCapPresTableName.empty() ||
                         m_imbibitionWettingIntermediateCapPresTableName == m_drainageWettingIntermediateCapPresTableName)
                        ? 0 : 1;
                    m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = 
                        (m_imbibitionNonWettingIntermediateCapPresTableName.empty() ||
                         m_imbibitionNonWettingIntermediateCapPresTableName == m_drainageNonWettingIntermediateCapPresTableName)
                        ? 0 : 1;
                }
            }

            // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

            m_wettingNonWettingCapillaryPressureKernelWrappers.clear();
            m_wettingIntermediateCapillaryPressureKernelWrappers.clear();
            m_nonWettingIntermediateCapillaryPressureKernelWrappers.clear();

            if (numPhases == 2) {
                GEOS_THROW_IF(m_drainageWettingNonWettingCapPresTableName.empty(),
                               GEOS_FMT("{}: drainageWettingNonWettingCapPressureTableName is empty for two-phase flow",
                                         getFullName()),
                               InputError);

                GEOS_THROW_IF(!functionManager.hasGroup(m_drainageWettingNonWettingCapPresTableName),
                               GEOS_FMT("{}: the table function named '{}' could not be found",
                                         getFullName(),
                                         m_drainageWettingNonWettingCapPresTableName),
                               InputError);
                TableFunction const &drainageCapPresTable = functionManager.getGroup<TableFunction>(
                        m_drainageWettingNonWettingCapPresTableName);
                m_wettingNonWettingCapillaryPressureKernelWrappers.emplace_back(
                        drainageCapPresTable.createKernelWrapper());

                string const &imbibitionTableName = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING]
                                                                     ? m_imbibitionWettingNonWettingCapPresTableName
                                                                     : m_drainageWettingNonWettingCapPresTableName;
                GEOS_THROW_IF(imbibitionTableName.empty(),
                               GEOS_FMT("{}: imbibition table name is empty for two-phase flow",
                                         getFullName()),
                               InputError);
                GEOS_THROW_IF(!functionManager.hasGroup(imbibitionTableName),
                               GEOS_FMT("{}: the table function named '{}' could not be found",
                                         getFullName(),
                                         imbibitionTableName),
                               InputError);
                TableFunction const &imbibitionWettingCapPresTable = functionManager.getGroup<TableFunction>(
                                imbibitionTableName);
                m_wettingNonWettingCapillaryPressureKernelWrappers.emplace_back(
                        imbibitionWettingCapPresTable.createKernelWrapper());
                
                GEOS_THROW_IF(m_wettingNonWettingCapillaryPressureKernelWrappers.size() != 2,
                               GEOS_FMT("{}: Expected 2 kernel wrappers after creation, but got {}",
                                         getFullName(),
                                         m_wettingNonWettingCapillaryPressureKernelWrappers.size()),
                               InputError);

            } else if (numPhases == 3) {
                TableFunction const &drainageWICapPres = functionManager.getGroup<TableFunction>(
                        m_drainageWettingIntermediateCapPresTableName);
                m_wettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
                        drainageWICapPres.createKernelWrapper());

                TableFunction const &drainageNWICapPres = functionManager.getGroup<TableFunction>(
                        m_drainageNonWettingIntermediateCapPresTableName);
                m_nonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
                        drainageNWICapPres.createKernelWrapper());

                TableFunction const &imbibitionWICapPres = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING]
                                                           ? functionManager.getGroup<TableFunction>(
                                m_imbibitionWettingIntermediateCapPresTableName)
                                                           : functionManager.getGroup<TableFunction>(
                                m_drainageWettingIntermediateCapPresTableName);
                m_wettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
                        imbibitionWICapPres.createKernelWrapper());

                TableFunction const &imbibitionNWICapPres = m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING]
                                                            ? functionManager.getGroup<TableFunction>(
                                m_imbibitionNonWettingIntermediateCapPresTableName)
                                                            : functionManager.getGroup<TableFunction>(
                                m_drainageNonWettingIntermediateCapPresTableName);
                m_nonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
                        imbibitionNWICapPres.createKernelWrapper());
            }

        }

///kernel ctor
        TableCapillaryPressureHysteresis::KernelWrapper::KernelWrapper(
                arrayView1d<const TableFunction::KernelWrapper> const &wettingNonWettingCapillaryPressureKernelWrappers,
                arrayView1d<const TableFunction::KernelWrapper> const &wettingIntermediateCapillaryPressureKernelWrappers,
                arrayView1d<const TableFunction::KernelWrapper> const &nonWettingIntermediateCapillaryPressureKernelWrappers,
                const arrayView1d<const geos::integer> &phaseHasHysteresis,
                const arrayView1d<const geos::real64> &landParam,
                const real64 & jerauldParam_a,
                const real64 & jerauldParam_b,
                const real64 & killoughCurvaturePcParam,
                const geos::real64 &phaseIntermediateMinVolFraction,
                const KilloughHysteresis::HysteresisCurve &wettingCurve,
                const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                const arrayView2d<const geos::real64, compflow::USD_PHASE> &phaseMinHistoricalVolFraction,
                const arrayView2d<const geos::real64, compflow::USD_PHASE> &phaseMaxHistoricalVolFraction,
                arrayView1d<integer const> const &phaseTypes,
                arrayView1d<integer const> const &phaseOrder,
                arrayView1d<integer> const &mode,
                arrayView3d<real64, cappres::USD_CAPPRES> const &phaseTrapped,
                const arrayView3d<geos::real64, relperm::USD_RELPERM> &phaseCapPressure,
                const arrayView4d<geos::real64, relperm::USD_RELPERM_DS> &dPhaseCapPressure_dPhaseVolFrac)
                :
                CapillaryPressureBaseUpdate(phaseTypes,
                                            phaseOrder,
                                            phaseTrapped,
                                            phaseCapPressure,
                                            dPhaseCapPressure_dPhaseVolFrac),
                m_wettingNonWettingCapillaryPressureKernelWrappers(wettingNonWettingCapillaryPressureKernelWrappers),
                m_wettingIntermediateCapillaryPressureKernelWrappers(
                        wettingIntermediateCapillaryPressureKernelWrappers),
                m_nonWettingIntermediateCapillaryPressureKernelWrappers(
                        nonWettingIntermediateCapillaryPressureKernelWrappers),
                m_phaseHasHysteresis(phaseHasHysteresis),
                m_landParam(landParam),
                m_jerauldParam_a(jerauldParam_a),
                m_jerauldParam_b(jerauldParam_b),
                m_killoughCurvatureParamCapPres(killoughCurvaturePcParam),
                m_phaseIntermediateMinVolFraction(phaseIntermediateMinVolFraction),
                m_wettingCurve(wettingCurve),
                m_nonWettingCurve(nonWettingCurve),
                m_phaseMinHistoricalVolFraction(phaseMinHistoricalVolFraction),
                m_phaseMaxHistoricalVolFraction(phaseMaxHistoricalVolFraction),
                m_mode(mode) {}


        REGISTER_CATALOG_ENTRY(ConstitutiveBase, TableCapillaryPressureHysteresis, std::string const &, Group * const)

    }
} // geos
