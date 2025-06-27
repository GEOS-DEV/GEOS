/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GravitySolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYSOLVERBASE_HPP_

#include "mesh/MeshFields.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"


namespace geos
{

class GravitySolverBase : public PhysicsSolverBase
{
public:

  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  GravitySolverBase() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  GravitySolverBase( const std::string & name,
                     Group * const parent );

  virtual ~GravitySolverBase() override;

  GravitySolverBase( GravitySolverBase const & ) = delete;
  GravitySolverBase( GravitySolverBase && ) = default;

  GravitySolverBase & operator=( GravitySolverBase const & ) = delete;
  GravitySolverBase & operator=( GravitySolverBase && ) = delete;

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;
  virtual void initializePostInitialConditionsPreSubGroups() override;


  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;


  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;


  /// Bind data between input XML file and source code
  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * modeString() { return "mode"; }
    static constexpr char const * stationCoordinatesString() { return "stationCoordinates"; }
    static constexpr char const * outputPropertyLogString() { return "outputPropertyLog"; }
    static constexpr char const * outputGzLogString() { return "outputGzLog"; }
    static constexpr char const * outputGzString() { return "outputGz"; }
    static constexpr char const * outputGzBasenameString() { return "outputGzBasename"; }
    static constexpr char const * outputAdjointLogString() { return "outputAdjointLog"; }
    static constexpr char const * residueString() { return "residue"; }
    static constexpr char const * gzAtStationsString() { return "gzAtStations"; }
  };

  void reinit() override final;


protected:

  void saveGz( real64 const & time_n,
               integer const cycleNumber,
               string const basename,
               arrayView1d< real64 > const & gzAtStations );

  virtual real64 explicitStepModeling( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition & domain ) = 0;

  virtual real64 explicitStepAdjoint( real64 const & time_n,
                                      real64 const & dt,
                                      integer const cycleNumber,
                                      DomainPartition & domain ) = 0;


  enum class GravityMode { Modeling, Adjoint };
  static const std::unordered_map< std::string, GravityMode > modeMap;


  GravityMode parseMode( const std::string & modeStr )
  {
    auto it = modeMap.find( modeStr );
    if( it != modeMap.end())
      return it->second;
    throw std::invalid_argument( "Invalid mode string: " + modeStr );
  }



  std::string m_modeString;

  GravityMode m_mode;


  /// Coordinates of the gravimeter stations
  array2d< real64 > m_stationCoordinates;

  /// Dump vertical component Gz to disk
  localIndex m_outputGz;
  string m_outputGzBasename;

  /// Residue observed at station (only for adjoint computation)
  array1d< real64 > m_residue;

  /// Gz component recorded at stations
  array1d< real64 > m_gzAtStations;
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYSOLVERBASE_HPP_
