/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_



#include <string>
#include <limits>

#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "../dataRepository/ManagedGroup.hpp"
#include "../dataRepository/ExecutableGroup.hpp"
#include "common/DataTypes.hpp"
#include "mesh/MeshBody.hpp"


namespace geosx
{

class DomainPartition;

namespace systemSolverInterface
{
class EpetraBlockSystem;
class LinearSolverWrapper;
}
class SystemSolverParameters;


namespace dataRepository
{
namespace keys
{
string const courant = "courant";
string const maxDt   = "maxDt";
}
}

class SolverBase : public ExecutableGroup
{
public:

  explicit SolverBase( std::string const & name,
                       ManagedGroup * const parent );

  virtual ~SolverBase() override;

  static string CatalogName() { return "SolverBase"; }

  SolverBase() = default;
  SolverBase( SolverBase const & ) = default;
  SolverBase( SolverBase &&) = default;
  SolverBase& operator=( SolverBase const & ) = default;
  SolverBase& operator=( SolverBase&& ) = default;


//  virtual void Registration( dataRepository::WrapperCollection& domain );


  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::ManagedGroup * domain );


  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        int const cycleNumber,
                        dataRepository::ManagedGroup * domain ) override;


  /**
   * @brief function to perform setup for implicit timestep
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   */
  virtual void ImplicitStepSetup( real64 const& time_n,
                                  real64 const& dt,
                                  DomainPartition * const domain );


  virtual real64 NonlinearImplicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition * const domain );

  virtual real64 LinearImplicitStep(real64 const & time_n,
                                    real64 const & dt,
                                    integer const cycleNumber,
                                    DomainPartition * const domain );

  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param domain the domain partition
   * @param blockSystem the entire block system
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @return the residual for convergence evaluation
   */
  virtual real64 Assemble ( DomainPartition * const domain,
                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            real64 const time,
                            real64 const dt );

  virtual void AllocateTempVars();

  virtual void DeallocateTempVars();

  virtual void ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                    real64 const scalingFactor,
                                    localIndex const dofOffset,
                                    DomainPartition * const nodeManager );

  /**
   * @brief function to perform cleanup for implicit timestep
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   */
  virtual void ImplicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition * const domain );




  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  using CatalogInterface = cxx_utilities::CatalogInterface< SolverBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  struct viewKeyStruct
  {
    constexpr static auto verboseLevelString = "verboseLevel";
    constexpr static auto gravityVectorString = "gravityVector";

  } viewKeys;

  struct groupKeyStruct
  {
    dataRepository::GroupKey systemSolverParameters = { "SystemSolverParameters" };
  } groupKeys;


  /**
   * accessor for the system solver parameters.
   * @return
   */
  SystemSolverParameters * getSystemSolverParameters();

  R1Tensor const & getGravityVector() const { return m_gravityVector; }
  R1Tensor       & getGravityVector()       { return m_gravityVector; }
  R1Tensor const * globalGravityVector() const;
  integer verboseLevel() const { return m_verboseLevel; }

protected:
  /// This is a wrapper for the linear solver package
  systemSolverInterface::LinearSolverWrapper * m_linearSolverWrapper;

  /// this is a block structured linear system object used to hold the system
  systemSolverInterface::EpetraBlockSystem * m_linearSystem;

private:
  integer m_verboseLevel = 0;
  R1Tensor m_gravityVector;



};



} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
