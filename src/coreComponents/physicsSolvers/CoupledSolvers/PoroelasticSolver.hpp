/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file PoroelasticSolver.hpp
 *
 */

#ifndef POROELASTICSOLVER_HPP_
#define POROELASTICSOLVER_HPP_

#include "../SolverBase.hpp"

namespace geosx
{

class PoroelasticSolver : public SolverBase
{
public:
  PoroelasticSolver( const std::string& name,
                     Group * const parent );
  ~PoroelasticSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "Poroelastic"; }

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override final;

  virtual void
  ImplicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition * const domain ) override final;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition * const domain ) override;

  void UpdateDeformationForCoupling( DomainPartition * const domain );

  real64 SplitOperatorStep( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition * const domain);


  enum class couplingTypeOption : int
  {
    FixedStress,
    TightlyCoupled
  };



  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto couplingTypeOptionString = "couplingTypeOptionEnum";
    constexpr static auto couplingTypeOptionStringString = "couplingTypeOption";

    constexpr static auto totalMeanStressString = "totalMeanStress";
    constexpr static auto oldTotalMeanStressString = "oldTotalMeanStress";

    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto fluidSolverNameString = "fluidSolverName";
  } poroElasticSolverViewKeys;


protected:
  virtual void PostProcessInput() override final;

  virtual void InitializePostInitialConditions_PreSubGroups(dataRepository::Group * const problemManager) override final;


private:

  string m_solidSolverName;
  string m_flowSolverName;
  string m_couplingTypeOptionString;
  couplingTypeOption m_couplingTypeOption;

};

} /* namespace geosx */

#endif /* POROELASTICSOLVER_HPP_ */
