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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_DIETERICH_SEISMICITY_RATE_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_DIETERICH_SEISMICITY_RATE_HPP_

#include "physicsSolvers/inducedSeismicity/SeismicityRateBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

using namespace fields;


/**
 * @class DieterichSeismicityRate
 *
 * @brief Solving the ODE for seismicity rate from Dieterich, 1994
 *
 * @details This solver finds a solution R(x, t) - the seismicity rate - to the ordinary differential equation (ODE)
 * formulated by Dieterich, 1994 given a certain stressing history. The stressing history can consist
 * of mechanical stresses and pore pressure. The solver class includes a member variable
 * pointing to the stress solver that is specified in the XML file. SolverStep for the
 * stress solver is then called in the SolverStep function for the seismicity rate, to take
 * the updated stress history as the input.
 *
 * Solving the ODE is currently implemented by computing the closed-form integral solution
 * to the ODE which involves numerical calculation of an integral of a stress functional.
 * We initially solve for the log of the seismicity rate in order to avoid overflow that
 * typically occurs in the exponential of the stress history.
 */
class DieterichSeismicityRate : public SeismicityRateBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  DieterichSeismicityRate() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  DieterichSeismicityRate( const string & name,
                           Group * const parent );

  /// Destructor
  virtual ~DieterichSeismicityRate() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "DieterichSeismicityRate"; }

   /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  /**
   * @brief single step advance in computing the seismicity rate based on
   *  stress history according to closed form integral solution (Heimisson & Segall, 2018)
   *  to the ODE formulated by Dieterich, 1994
   * @param time_n time at previous converged step
   * @param dt time step size
   * @param subRegion ElementSubRegionBase to compute the solution in
   */
  void integralSolverStep( real64 const & time_n,
                           real64 const & dt,
                           ElementSubRegionBase & subRegion );

  /**@}*/

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * directEffectString() { return "directEffect"; }
    static constexpr char const * backgroundStressingRateString() { return "backgroundStressingRate"; }
  };

  virtual void initializePreSubGroups() override;

private:

  /**
   * @brief Checks stress rate that is argument of exponential in stress functional,
   *  called in integralSolverStep. If too large as to cause overflow, checks various
   *  conditions to give informative output message to user.
   */
  void checkExpArgument( real64 arg );

  real64 m_directEffect;
  real64 m_backgroundStressingRate;

  struct solverHelper
  {
    // Constructor
    solverHelper( ElementSubRegionBase & subRegion )
    {
      // Retrieve field variables
      R = subRegion.getField< inducedSeismicity::seismicityRate >();
      logDenom = subRegion.getField< inducedSeismicity::logDenom >();

      sigViews[0] = subRegion.getField< inducedSeismicity::initialProjectedNormalTraction >();
      sigViews[1] = subRegion.getField< inducedSeismicity::projectedNormalTraction_n >();
      sigViews[2] = subRegion.getField< inducedSeismicity::projectedNormalTraction >();

      tauViews[0] = subRegion.getField< inducedSeismicity::initialProjectedShearTraction >();
      tauViews[1] = subRegion.getField< inducedSeismicity::projectedShearTraction_n >();
      tauViews[2] = subRegion.getField< inducedSeismicity::projectedShearTraction >();

      if( subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::fluidNamesString() ) )
      {
        flowExists = true;

        pViews[0] = subRegion.getField< flow::initialPressure >();
        pViews[1] = subRegion.getField< flow::pressure_n >();
        pViews[2] = subRegion.getField< flow::pressure >();
      }
    }

    void computeSeismicityRate( ElementSubRegionBase & subRegion,
                                real64 const & time_n, real64 const & dt,
                                real64 const directEffectValue, real64 backgroundStressingRateValue )
    {
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        // arguments of stress exponential at current and previous time step
        real64 g = (tauViews[2][k] + backgroundStressingRateValue*(time_n+dt))/(directEffectValue*getEffectiveNormalTraction( 2, k ))
                   - tauViews[0][k]/(directEffectValue*getEffectiveNormalTraction( 0, k ));
        real64 g_n = (tauViews[1][k] + backgroundStressingRateValue*time_n)/(directEffectValue*getEffectiveNormalTraction( 1, k ))
                     - tauViews[0][k]/(directEffectValue*getEffectiveNormalTraction( 0, k ));

        // checkExpArgument();

        // Compute the difference of the log of the denominator of closed for integral solution.
        // This avoids directly computing the exponential of the current stress state which is more prone to overflow.
        logDenom[k] += std::log( 1 + dt/(2*(directEffectValue*getEffectiveNormalTraction( 0, k )/backgroundStressingRateValue))
                                 *(std::exp( g - logDenom[k] ) + std::exp( g_n - logDenom[k] ) ));

        // Convert log seismicity rate to raw value
        R[k] = LvArray::math::exp( g - logDenom[k] );
      } );
    }

    real64 getEffectiveNormalTraction( size_t const nIndex, localIndex const k )
    {
      if( flowExists )
      {
        return -sigViews[nIndex][k]-pViews[nIndex][k];
      }
      else
      {
        return -sigViews[nIndex][k];
      }
    }

    arrayView1d< real64 > R;
    arrayView1d< real64 > logDenom;

    arrayView1d< real64 const > sigViews[ 3 ];
    arrayView1d< real64 const > tauViews[ 3 ];
    arrayView1d< real64 const > pViews[ 3 ];

    bool flowExists = false;
  };

};

}/* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_DIETERICH_SEISMICITY_RATE_HPP_ */
