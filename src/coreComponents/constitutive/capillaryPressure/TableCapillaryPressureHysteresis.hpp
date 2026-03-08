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


#ifndef GEOS_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
#define GEOS_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "functions/TableFunction.hpp"

#include "constitutive/KilloughHysteresis.hpp"
#include "CapillaryPressureFields.hpp"
#include "common/DataLayouts.hpp"

namespace geos {


    namespace constitutive {

        class TableCapillaryPressureHysteresis : public CapillaryPressureBase {
//  /// useful constant
//  static constexpr real64 CAP_INF = 1e9;
////          std::numeric_limits< real64 >::max();
//  static constexpr real64 CAP_INF_DERIV = 1e9;
////          std::numeric_limits< real64 >::max();

            typedef fields::cappres::ModeIndexType ModeIndexType;

        public:

            /// order of the phase properties for three-phase flow
            struct ThreePhasePairPhaseType {
                enum : integer {
                    INTERMEDIATE_WETTING = 0,   ///< index for intermediate-wetting
                    INTERMEDIATE_NONWETTING = 1 ///< index for intermediate-non-wetting
                };
            };


            TableCapillaryPressureHysteresis(std::string const &name,
                                             dataRepository::Group *const parent);

            static std::string catalogName() { return "TableCapillaryPressureHysteresis"; }

            virtual string getCatalogName() const override { return catalogName(); }

            ///Kernel
            class KernelWrapper final : public CapillaryPressureBaseUpdate {
            public:

                KernelWrapper(
                        arrayView1d<TableFunction::KernelWrapper const> const &wettingNonWettingCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &inverseWettingNonWettingCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &wettingIntermediateCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &inverseWettingIntermediateCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &nonWettingIntermediateCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &inverseNonWettingIntermediateCapillaryPressureKernelWrappers,
                        arrayView1d<integer const> const &phaseHasHysteresis,
                        arrayView1d<real64 const> const &landParam,
                        real64 const &jerauldParam_a,
                        real64 const &jerauldParam_b,
                        real64 const &killoughCurvaturePcParameter,
                        real64 const &phaseIntermediateMinVolFraction,
                        KilloughHysteresis::HysteresisCurve const &wettingCurve,
                        KilloughHysteresis::HysteresisCurve const &nonWettingCurve,
                        arrayView2d<real64 const, compflow::USD_PHASE> const &phaseMinHistoricalVolFraction,
                        arrayView2d<real64 const, compflow::USD_PHASE> const &phaseMaxHistoricalVolFraction,
                        arrayView2d<real64, compflow::USD_PHASE> &phaseMode2PeakVolFraction,
                        arrayView1d<integer const> const &phaseTypes,
                        arrayView1d<integer const> const &phaseOrder,
                        arrayView1d<integer> const &mode,
                        arrayView3d<real64, relperm::USD_RELPERM> const &phaseTrappedVolFrac,
                        arrayView3d<real64, relperm::USD_RELPERM> const &phaseCapPressure,
                        arrayView4d<real64, relperm::USD_RELPERM_DS> const &dPhaseCapPressure_dPhaseVolFrac);

                //actual workers
                GEOS_HOST_DEVICE
                void computeBoundCapillaryPressure(TableFunction::KernelWrapper const &drainageCapPressureWrapper,
                                                   real64 const &phaseVolFraction,
                                                   real64 &phaseCapPressure,
                                                   real64 &dPhaseCapPressure_dPhaseVolFrac) const;

                GEOS_HOST_DEVICE
                void
                computeImbibitionWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &wettingKernelWapper,
                        const KilloughHysteresis::HysteresisCurve &wettingCurve,
                        const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMinHistoricalVolFraction,
                        const geos::real64 &phaseMaxHistoricalVolFraction,
                        const geos::real64 &phaseMode2PeakVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;

                //two phase flow overload
                GEOS_HOST_DEVICE
                void
                computeImbibitionWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &wettingKernelWapper,
                        const KilloughHysteresis::HysteresisCurve &wettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMinHistoricalVolFraction,
                        const geos::real64 &phaseMaxHistoricalVolFraction,
                        const geos::real64 &phaseMode2PeakVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;

                GEOS_HOST_DEVICE
                void
                computeImbibitionNonWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
                        const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                        const KilloughHysteresis::HysteresisCurve &wettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMaxHistoricalVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;

                //2phase flow overload
                GEOS_HOST_DEVICE
                void
                computeImbibitionNonWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
                        const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMaxHistoricalVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;



                //wrapper call wrt number of phase
                GEOS_HOST_DEVICE
                void computeTwoPhaseWetting(integer const ipWetting,
                                            integer const ipNonWetting,
                                            arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                            arraySlice1d<real64 const,
                                                    compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                            arraySlice1d<real64 const,
                                                    compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                            arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseTrappedVolFrac,
                                            arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseCapPressure,
                                            arraySlice2d<real64,
                                                    relperm::USD_RELPERM_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                            ModeIndexType &mode,
                                            arraySlice1d<real64, compflow::USD_PHASE - 1> &phaseMode2PeakVolFraction) const;


                GEOS_HOST_DEVICE
                void computeTwoPhaseNonWetting(integer const ipWetting,
                                               integer const ipNonWetting,
                                               arraySlice1d<real64 const,
                                                       compflow::USD_PHASE - 1> const &phaseVolFraction,
                                               arraySlice1d<real64 const,
                                                       compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                               arraySlice1d<real64 const,
                                                       compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                               arraySlice1d<real64,
                                                       relperm::USD_RELPERM - 2> const &phaseTrappedVolFrac,
                                               arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseCapPressure,
                                               arraySlice2d<real64, relperm::USD_RELPERM_DS -
                                                                    2> const &dPhaseCapPressure_dPhaseVolFrac,
                                               ModeIndexType &mode) const;

                GEOS_HOST_DEVICE
                void computeThreePhase(integer const ipWetting,
                                       integer const ipInter,
                                       integer const ipNonWetting,
                                       arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                       arraySlice1d<real64 const,
                                               compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                       arraySlice1d<real64 const,
                                               compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                       arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseTrappedVolFrac,
                                       arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseCapPressure,
                                       arraySlice2d<real64,
                                               relperm::USD_RELPERM_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                       ModeIndexType &mode,
                                       arraySlice1d<real64, compflow::USD_PHASE - 1> &phaseMode2PeakVolFraction) const;

                //uppermost call-wrappers
                // Standard 3-argument compute method for compatibility with InverseCapillaryPressure
                GEOS_HOST_DEVICE
                void compute(arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                             arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                             arraySlice2d<real64, cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac) const;

                GEOS_HOST_DEVICE
                virtual void compute(arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                     arraySlice1d<real64 const,
                                             compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                     arraySlice1d<real64 const,
                                             compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                     arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseTrappedVolFrac,
                                     arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                                     arraySlice2d<real64,
                                             cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                     ModeIndexType &mode,
                                     arraySlice1d<real64, compflow::USD_PHASE - 1> &phaseMode2PeakVolFraction) const;

                GEOS_HOST_DEVICE
                virtual void update(localIndex const k,
                                    localIndex const q,
                                    arraySlice1d<real64 const,
                                            compflow::USD_PHASE - 1> const &phaseVolFraction) const override;

                /**
                 * @brief Compute phase volume fraction from capillary pressure (inverse operation).
                 * @param phaseVolFraction [out] Computed phase volume fractions
                 * @param phaseMaxHistoricalVolFraction [in] Maximum historical phase volume fractions
                 * @param phaseMinHistoricalVolFraction [in] Minimum historical phase volume fractions
                 * @param phaseTrappedVolFrac [in] Trapped phase volume fractions
                 * @param phaseCapPressure [in] Target capillary pressures (input)
                 * @param dPhaseCapPressure_dPhaseVolFrac [out] Derivatives of capillary pressure w.r.t. phase volume fraction
                 * @param mode [in] Hysteresis mode (DRAINAGE, IMBIBITION, DRAINAGE_TO_IMBIBITION, IMBIBITION_TO_DRAINAGE)
                 * 
                 * Uses Newton-Raphson iteration to invert the capillary pressure function.
                 */
                GEOS_HOST_DEVICE
                void computeInv(arraySlice1d<real64, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMode2PeakVolFraction,
                                arraySlice1d<real64 const, cappres::USD_CAPPRES - 2> const &phaseTrappedVolFrac,
                                arraySlice1d<real64 const, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                                arraySlice2d<real64, cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                fields::cappres::ModeIndexType const &mode) const;

                /**
                 * @brief Evaluate the raw drainage or imbibition table at a given saturation,
                 *        bypassing the hysteresis mode-transition logic in compute().
                 * @param tableIdx Table index (0 = DRAINAGE, 1 = IMBIBITION)
                 * @param phaseVolFraction Phase volume fraction to evaluate at
                 * @param phaseCapPressure [out] Capillary pressure from the table
                 * @param dPhaseCapPressure_dPhaseVolFrac [out] Derivative of Pc w.r.t. S
                 */
                GEOS_HOST_DEVICE
                void computeRawTablePc( integer const tableIdx,
                                        real64 const & phaseVolFraction,
                                        real64 & phaseCapPressure,
                                        real64 & dPhaseCapPressure_dPhaseVolFrac ) const
                {
                    computeBoundCapillaryPressure(
                        m_wettingNonWettingCapillaryPressureKernelWrappers[tableIdx],
                        phaseVolFraction, phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac );
                }

                /**
                 * @brief Compute the actual Pc range [Pc_lo, Pc_hi] for a scanning curve.
                 *
                 * For Mode 2 (DRAINAGE_TO_IMBIBITION):
                 *   Pc_hi = Pc_drainage(Shy)          (departure point)
                 *   Pc_lo = Pc_imbibition(Swma)       (endpoint, where F=1)
                 * For Mode 3 (IMBIBITION_TO_DRAINAGE):
                 *   Pc_lo = Pc_imbibition(Shy)         (departure point)
                 *   Pc_hi = Pc_drainage(Scrt)          (endpoint, where F=1)
                 *
                 * @param phaseMinHistVolFrac  Minimum historical volume fraction (Shy for Mode 2)
                 * @param phaseMaxHistVolFrac  Maximum historical volume fraction (Shy for Mode 3)
                 * @param mode                 Hysteresis mode
                 * @param Pc_lo  [out] Lower bound of scanning curve Pc range
                 * @param Pc_hi  [out] Upper bound of scanning curve Pc range
                 */
                GEOS_HOST_DEVICE
                void computeScanningCurvePcRange( real64 const phaseMinHistVolFrac,
                                                  real64 const phaseMaxHistVolFrac,
                                                  ModeIndexType const mode,
                                                  real64 & Pc_lo,
                                                  real64 & Pc_hi ) const
                {
                    integer const ipWater = 0;
                    real64 dPc_dummy;

                    if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION )
                    {
                        // Shy = min historical saturation (departure from drainage)
                        real64 const Smin_curve = m_wettingCurve.oppositeBoundPhaseVolFraction;
                        real64 const Shy = LvArray::math::max( phaseMinHistVolFrac, Smin_curve );

                        // Compute Scrt (trapped saturation) from Land model
                        real64 Scrt = 0.0;
                        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(
                            m_wettingCurve, Shy, m_landParam[ipWater],
                            m_jerauldParam_a, m_jerauldParam_b, Scrt );

                        // For wetting phase, Swma = Scrt (= 1 - (1 - Scrt))
                        real64 const Swma = Scrt;

                        // Pc_hi = Pc_drainage(Shy) — departure point, highest Pc on scanning curve
                        computeRawTablePc( 0 /* drainage */, Shy, Pc_hi, dPc_dummy );

                        // Pc_lo = Pc_imbibition(Swma) — endpoint where F=1, lowest Pc on scanning curve
                        computeRawTablePc( 1 /* imbibition */, Swma, Pc_lo, dPc_dummy );
                    }
                    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE ||
                             mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
                    {
                        // Shy = max historical saturation (departure from imbibition)
                        real64 const Smax_curve = m_wettingCurve.drainageExtremaPhaseVolFraction;
                        real64 const Shy = LvArray::math::min( phaseMaxHistVolFrac, Smax_curve );

                        // Compute Scrt
                        real64 Scrt = 0.0;
                        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction(
                            m_wettingCurve, Shy, m_landParam[ipWater],
                            m_jerauldParam_a, m_jerauldParam_b, Scrt );

                        // Pc_lo = Pc_imbibition(Shy) — departure point, lowest Pc on scanning curve
                        computeRawTablePc( 1 /* imbibition */, Shy, Pc_lo, dPc_dummy );

                        // Pc_hi = Pc_drainage(Scrt) — endpoint where F=1, highest Pc on scanning curve
                        computeRawTablePc( 0 /* drainage */, Scrt, Pc_hi, dPc_dummy );
                    }
                    else
                    {
                        // Fallback: use raw table full range
                        computeRawTablePc( 0 /* drainage */, 0.0, Pc_hi, dPc_dummy );
                        computeRawTablePc( 1 /* imbibition */, 1.0, Pc_lo, dPc_dummy );
                    }

                    // Ensure Pc_lo <= Pc_hi
                    if( Pc_lo > Pc_hi )
                    {
                        real64 tmp = Pc_lo;
                        Pc_lo = Pc_hi;
                        Pc_hi = tmp;
                    }
                }

            private:

                /**
                 * @brief Helper function to compute Pc(S) for given S, mode, and historical values.
                 * @param S Phase volume fraction
                 * @param mode Hysteresis mode
                 * @param ipPhase Phase index
                 * @param phaseMinHistoricalVolFraction Minimum historical volume fraction
                 * @param phaseMaxHistoricalVolFraction Maximum historical volume fraction
                 * @param capPresKernelWrappers Capillary pressure kernel wrappers
                 * @param wettingCurve Wetting curve data
                 * @param nonWettingCurve Non-wetting curve data
                 * @param landParam Land parameter
                 * @param phaseIntermediateMinVolFraction Intermediate phase minimum volume fraction
                 * @param killoughCurvatureParam Killough curvature parameter
                 * @param jerauldParam_a Jerauld parameter a
                 * @param jerauldParam_b Jerauld parameter b
                 * @param isWettingPhase Whether this is the wetting phase
                 * @param precomputedScrt Precomputed trapped critical saturation (optional, use -1.0 to compute)
                 * @param precomputedDenomF Precomputed denominator for F calculation (optional, use 0.0 to compute)
                 * @param precomputedShy Precomputed historical saturation Shy (optional, use -1.0 to compute)
                 * @return Computed capillary pressure
                 * 
                 * Used for Newton-Raphson inversion in computeInv.
                 * If precomputed values are provided (precomputedScrt >= 0, precomputedDenomF != 0, precomputedShy >= 0),
                 * they will be used instead of recomputing them.
                 */
                GEOS_HOST_DEVICE
                real64 computeCapillaryPressureForSaturation(
                        real64 const S,
                        fields::cappres::ModeIndexType const &mode,
                        integer const ipPhase,
                        real64 const &phaseMinHistoricalVolFraction,
                        real64 const &phaseMaxHistoricalVolFraction,
                        real64 const &phaseMode2PeakVolFraction,
                        arrayView1d<TableFunction::KernelWrapper const> const &capPresKernelWrappers,
                        KilloughHysteresis::HysteresisCurve const &wettingCurve,
                        KilloughHysteresis::HysteresisCurve const &nonWettingCurve,
                        real64 const &landParam,
                        real64 const &phaseIntermediateMinVolFraction,
                        real64 const &killoughCurvatureParam,
                        real64 const &jerauldParam_a,
                        real64 const &jerauldParam_b,
                        bool const isWettingPhase,
                        real64 const precomputedScrt = -1.0,
                        real64 const precomputedDenomF = 0.0,
                        real64 const precomputedShy = -1.0) const;

                static constexpr real64 flowReversalBuffer = KilloughHysteresis::flowReversalBuffer;
//    ModeIndexType& m_mode;

                //2p
                arrayView1d<TableFunction::KernelWrapper const> const m_wettingNonWettingCapillaryPressureKernelWrappers;
                arrayView1d<TableFunction::KernelWrapper const> const m_inverseWettingNonWettingCapillaryPressureKernelWrappers;
                //3p
                arrayView1d<TableFunction::KernelWrapper const> const m_wettingIntermediateCapillaryPressureKernelWrappers;
                arrayView1d<TableFunction::KernelWrapper const> const m_inverseWettingIntermediateCapillaryPressureKernelWrappers;
                arrayView1d<TableFunction::KernelWrapper const> const m_nonWettingIntermediateCapillaryPressureKernelWrappers;
                arrayView1d<TableFunction::KernelWrapper const> const m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers;

                ///Land Coeff
                arrayView1d<integer const> m_phaseHasHysteresis;
                arrayView1d<real64 const> m_landParam;

                /// Parameter a introduced by Jerauld in the Land model
                const real64 m_jerauldParam_a;

                /// Parameter b introduced by Jerauld in the Land model
                const real64 m_jerauldParam_b;

                /// Curvature parameter in Killough wetting phase hysteresis (enpoints curvatures)
                const real64 m_killoughCurvatureParamCapPres;

                /// needed in 3p-wetting hysteresis as we need to get the max accessible pore space
                real64 const m_phaseIntermediateMinVolFraction;

                KilloughHysteresis::HysteresisCurve const m_wettingCurve;
                KilloughHysteresis::HysteresisCurve const m_nonWettingCurve;

                /// Minimum historical phase volume fraction for each phase
                arrayView2d<real64 const, compflow::USD_PHASE> m_phaseMinHistoricalVolFraction;

                /// Maximum historical phase volume fraction for each phase
                arrayView2d<real64 const, compflow::USD_PHASE> m_phaseMaxHistoricalVolFraction;

                /// Peak saturation reached during Mode 2 (DRAINAGE_TO_IMBIBITION) for each phase (mutable for updates)
                arrayView2d<real64, compflow::USD_PHASE> m_phaseMode2PeakVolFraction;

                // Drainage / Imbibition flags cellwise
                arrayView1d<ModeIndexType> m_mode;

            };

            /**
             * @brief Create an update kernel wrapper.
             * @return the wrapper
             */
            KernelWrapper createKernelWrapper();

            //might need it to be virtual one level higher --> from Killough/Hysteresis common class
            virtual void saveConvergedPhaseVolFractionState(
                    arrayView2d<real64 const, compflow::USD_PHASE> const &phaseVolFraction) const override;


            struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct {


                ///Land Coeff
                static constexpr char const *landParameterString() { return "landParameter"; }

                ///flag
                static constexpr char const *phaseHasHysteresisString() { return "phaseHasHysteresis"; }

                ///and packed curves data struct
                static constexpr char const *wettingCurveString() { return "wettingCurve"; };

                static constexpr char const *nonWettingCurveString() { return "nonWettingCurve"; };


                ///tables and assoc. wrappers
                //2phase
                static constexpr char const *
                drainageWettingNonWettingCapPresTableNameString() { return "drainageWettingNonWettingCapPressureTableName"; }

                static constexpr char const *
                imbibitionWettingNonWettingCapPresTableNameString() { return "imbibitionWettingNonWettingCapPressureTableName"; }

                //3phase
                static constexpr char const *
                drainageWettingIntermediateCapPresTableNameString() { return "drainageWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                drainageNonWettingIntermediateCapPresTableNameString() { return "drainageNonWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                imbibitionWettingIntermediateCapPresTableNameString() { return "imbibitionWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                imbibitionNonWettingIntermediateCapPresTableNameString() { return "imbibitionNonWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                wettingNonWettingCapillaryPressureKernelWrappersString() { return "wettingNonWettingCapillaryPressureKernelWrappers"; }

                static constexpr char const *
                wettingIntermediateCapillaryPressureKernelWrappersString() { return "wettingIntermediateCapillaryPressureKernelWrappers"; }

                static constexpr char const *
                nonWettingIntermediateCapillaryPressureKernelWrappersString() { return "nonWettingIntermediateCapillaryPressureKernelWrappers"; }

                //misc
                static constexpr char const *
                phaseIntermediateMinVolFractionString() { return "phaseIntermediateMinVolFraction"; }
                //to decide wheter drainage/drainage to imbibition or imbibition/imbibition to drainage
            };


        private:
            virtual void postProcessInput();

            virtual void initializePreSubGroups() override;

            void resizeFields(localIndex const size,
                              localIndex const numPts) override;


            /**
             * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
             */
            void createAllTableKernelWrappers();

            /**
             * @brief Compute the Land coefficient for the wetting and non-wetting phases
             */
            void computeLandCoefficient();

            ///data members


            //TODO impl
//  array1d< integer >  m_tCurveOption;

            KilloughHysteresis::HysteresisCurve m_wettingCurve;
            KilloughHysteresis::HysteresisCurve m_nonWettingCurve;

            ///tables
            //2p
            string m_drainageWettingNonWettingCapPresTableName;
            string m_imbibitionWettingNonWettingCapPresTableName;
            //3p
            string m_drainageWettingIntermediateCapPresTableName;
            string m_drainageNonWettingIntermediateCapPresTableName;
            string m_imbibitionWettingIntermediateCapPresTableName;
            string m_imbibitionNonWettingIntermediateCapPresTableName;
            // kernel wrappers
            /// Imbibition kernel wrappers for relative permeabilities in the following order:
            /// 0- drainage
            /// 1- imbibition (cf. struct ModeIndexType)
            //2p
            array1d<TableFunction::KernelWrapper> m_wettingNonWettingCapillaryPressureKernelWrappers;
            array1d<TableFunction::KernelWrapper> m_inverseWettingNonWettingCapillaryPressureKernelWrappers;
            //3p
            array1d<TableFunction::KernelWrapper> m_wettingIntermediateCapillaryPressureKernelWrappers;
            array1d<TableFunction::KernelWrapper> m_inverseWettingIntermediateCapillaryPressureKernelWrappers;
            array1d<TableFunction::KernelWrapper> m_nonWettingIntermediateCapillaryPressureKernelWrappers;
            array1d<TableFunction::KernelWrapper> m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers;
            
            // Store inverse tables to keep them alive
            std::vector<std::shared_ptr<TableFunction>> m_inverseTables;


            /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
            array1d<integer> m_phaseHasHysteresis;

            /// Trapping parameter from the Land model (typically called C)
            array1d<real64> m_landParam;

            /// Parameter a introduced by Jerauld in the Land model
            real64 m_jerauldParam_a;

            /// Parameter b introduced by Jerauld in the Land model
            real64 m_jerauldParam_b;

            /// Curvature parameter in Killough wetting phase hysteresis (Scanning curves curvatures)
            real64 m_killoughCurvatureParamCapPres;

            /// Cell-wise status imbibition, imbibitioon_to_drainage, ... etc
            array1d<integer> m_mode;

            // Max historical saturations
            /// Minimum historical phase volume fraction for each phase
            array2d<real64, compflow::LAYOUT_PHASE> m_phaseMinHistoricalVolFraction;

            /// Maximum historical phase volume fraction for each phase
            array2d<real64, compflow::LAYOUT_PHASE> m_phaseMaxHistoricalVolFraction;

            /// Peak saturation reached during Mode 2 (DRAINAGE_TO_IMBIBITION) for each phase
            array2d<real64, compflow::LAYOUT_PHASE> m_phaseMode2PeakVolFraction;

            //needed in hysteresis of wetting phase
            real64 m_phaseIntermediateMinVolFraction;

        };

        // Standard 3-argument compute method for compatibility with InverseCapillaryPressure
        // Uses zero/default values for historical fractions (no hysteresis in inverse computation)
        GEOS_HOST_DEVICE
        inline void TableCapillaryPressureHysteresis::KernelWrapper::compute(
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                arraySlice2d<real64, cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac
        ) const {
            // Create temporary arrays for historical fractions and trapped fraction
            // Initialize to zero/default values since inverse computation doesn't use hysteresis
            constexpr integer MAX_NUM_PHASES = CapillaryPressureBase::MAX_NUM_PHASES;
            real64 phaseMaxHistoricalVolFraction[MAX_NUM_PHASES]{};
            real64 phaseMinHistoricalVolFraction[MAX_NUM_PHASES]{};
            real64 phaseTrappedVolFrac[MAX_NUM_PHASES]{};
            real64 phaseMode2PeakVolFraction[MAX_NUM_PHASES]{};
            ModeIndexType mode = ModeIndexType::DRAINAGE; // Default to drainage mode
            
            // Create ArrayView from stack arrays, then get slices
            integer const numPhases = LvArray::integerConversion< integer >( phaseVolFraction.size() );
            localIndex dims[1] = { numPhases };
            localIndex strides[1] = { 1 };
            
            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseMaxHistSlice(
                phaseMaxHistoricalVolFraction, dims, strides );
            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseMinHistSlice(
                phaseMinHistoricalVolFraction, dims, strides );
            arraySlice1d< real64, cappres::USD_CAPPRES - 2 > phaseTrappedSlice(
                phaseTrappedVolFrac, dims, strides );
            arraySlice1d< real64, compflow::USD_PHASE - 1 > phaseMode2PeakSlice(
                phaseMode2PeakVolFraction, dims, strides );
            
            // Call the full compute method
            compute( phaseVolFraction,
                     phaseMaxHistSlice,
                     phaseMinHistSlice,
                     phaseTrappedSlice,
                     phaseCapPressure,
                     dPhaseCapPressure_dPhaseVolFrac,
                     mode,
                     phaseMode2PeakSlice );
        }

        GEOS_HOST_DEVICE
        inline void TableCapillaryPressureHysteresis::KernelWrapper::compute(
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseTrappedVolFrac,
                arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                arraySlice2d<real64, cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                ModeIndexType &mode,
                arraySlice1d<real64, compflow::USD_PHASE - 1> &phaseMode2PeakVolFraction
        ) const {
            // Early return if m_phaseOrder is empty or input arrays are empty
            if( m_phaseOrder.size() == 0 || 
                phaseVolFraction.size() == 0 ||
                phaseCapPressure.size() == 0 )
            {
                return;
            }

            LvArray::forValuesInSlice(dPhaseCapPressure_dPhaseVolFrac, [](real64 &val) { val = 0.0; });

            using PT = CapillaryPressureBase::PhaseType;
            // Check bounds before accessing m_phaseOrder
            integer const ipWater = ( PT::WATER < m_phaseOrder.size() ) ? m_phaseOrder[PT::WATER] : -1;
            integer const ipOil = ( PT::OIL < m_phaseOrder.size() ) ? m_phaseOrder[PT::OIL] : -1;
            integer const ipGas = ( PT::GAS < m_phaseOrder.size() ) ? m_phaseOrder[PT::GAS] : -1;

            if (ipWater >= 0 && ipOil >= 0 && ipGas >= 0) {
                computeThreePhase(ipWater, // wetting
                                  ipOil,   // intermediate
                                  ipGas,   // non-wetting
                                  phaseVolFraction,
                                  phaseMaxHistoricalVolFraction,
                                  phaseMinHistoricalVolFraction,
                                  phaseTrappedVolFrac,
                                  phaseCapPressure,
                                  dPhaseCapPressure_dPhaseVolFrac,
                                  mode,
                                  phaseMode2PeakVolFraction);

            } else if (ipWater < 0) {
                computeTwoPhaseNonWetting(ipOil, // leading
                                          ipGas, // deduced
                                          phaseVolFraction,
                                          phaseMaxHistoricalVolFraction,
                                          phaseMinHistoricalVolFraction,
                                          phaseTrappedVolFrac,
                                          phaseCapPressure,
                                          dPhaseCapPressure_dPhaseVolFrac,
                                          mode);
            } else if (ipOil < 0) {
                computeTwoPhaseWetting(ipWater, // leading
                                       ipGas,   // deduced
                                       phaseVolFraction,
                                       phaseMaxHistoricalVolFraction,
                                       phaseMinHistoricalVolFraction,
                                       phaseTrappedVolFrac,
                                       phaseCapPressure,
                                       dPhaseCapPressure_dPhaseVolFrac,
                                       mode,
                                       phaseMode2PeakVolFraction);
            } else if (ipGas < 0) {
                computeTwoPhaseWetting(ipWater, //leading
                                       ipOil,   //deduced
                                       phaseVolFraction,
                                       phaseMaxHistoricalVolFraction,
                                       phaseMinHistoricalVolFraction,
                                       phaseTrappedVolFrac,
                                       phaseCapPressure,
                                       dPhaseCapPressure_dPhaseVolFrac,
                                       mode,
                                       phaseMode2PeakVolFraction);
            }


        }

        GEOS_HOST_DEVICE
        inline void TableCapillaryPressureHysteresis::KernelWrapper::update(const geos::localIndex k,
                                                                            const geos::localIndex q,
                                                                            const arraySlice1d<const geos::real64,
                                                                                    compflow::USD_PHASE
                                                                                    - 1> &phaseVolFraction) const {
            // Create a reference to the Mode 2 peak slice for this element
            // m_phaseMode2PeakVolFraction uses compflow::LAYOUT_PHASE, so use compflow::USD_PHASE - 1
            arraySlice1d<real64, compflow::USD_PHASE - 1> phaseMode2PeakSlice = m_phaseMode2PeakVolFraction[k];
            compute(phaseVolFraction,
                    m_phaseMaxHistoricalVolFraction[k],
                    m_phaseMinHistoricalVolFraction[k],
                    m_phaseTrappedVolFrac[k][q],
                    m_phaseCapPressure[k][q],
                    m_dPhaseCapPressure_dPhaseVolFrac[k][q],
                    m_mode[k],
                    phaseMode2PeakSlice);
        }


    } //constitutive
} // geos

#endif //GEOS_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
