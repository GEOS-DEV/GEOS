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

namespace geos
{

//START_SPHINX_INCLUDE_BEGINCLASS
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

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override; 

//END_SPHINX_INCLUDE_BEGINCLASS

//START_SPHINX_INCLUDE_SOLVERINTERFACE
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  // void odeSolverStep( real64 const & time_n,
  //                     real64 const & dt,
  //                     integer const cycleNumber,
  //                     DomainPartition & domain );

  void integralSolverStep( real64 const & time_n,
                           real64 const & dt,
                           ElementSubRegionBase & subRegion );

  void updateMeanSolidStress( ElementSubRegionBase & subRegion );

//END_SPHINX_INCLUDE_SOLVERINTERFACE

  /**@}*/

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * directEffectString() { return "directEffect"; }
    static constexpr char const * backgroundStressingRateString() { return "backgroundStressingRate"; }
  };

  virtual void initializePreSubGroups() override;

protected:

  real64 m_directEffect;
  real64 m_backgroundStressingRate;
  
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_DIETERICH_SEISMICITY_RATE_HPP_ */
