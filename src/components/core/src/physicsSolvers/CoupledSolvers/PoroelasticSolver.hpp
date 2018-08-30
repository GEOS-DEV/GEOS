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
                     ManagedGroup * const parent );
  ~PoroelasticSolver();

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "Poroelastic"; }

  virtual void FillDocumentationNode() override final;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup ) override final;

  virtual void ImplicitStepSetup( real64 const& time_n,
                                  real64 const& dt,
                                  DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem) override final;

  virtual void ImplicitStepComplete( real64 const& time_n,
                                     real64 const& dt,
                                     DomainPartition * const domain) override final;

  virtual void FinalInitialization( dataRepository::ManagedGroup * const problemManager ) override final;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * domain ) override;

//  virtual real64 ExplicitStep( real64 const & time_n,
//                               real64 const & dt,
//                               integer const cycleNumber,
//                               DomainPartition * const domain );
//
//  virtual real64 NonlinearImplicitStep( real64 const & time_n,
//                                        real64 const & dt,
//                                        integer const cycleNumber,
//                                        DomainPartition * const domain,
//                                        systemSolverInterface::EpetraBlockSystem * const blockSystem );
//
//  virtual real64 LinearImplicitStep(real64 const & time_n,
//                                    real64 const & dt,
//                                    integer const cycleNumber,
//                                    DomainPartition * const domain,
//                                    systemSolverInterface::EpetraBlockSystem * const blockSystem );
//
//  virtual void ImplicitStepSetup( real64 const& time_n,
//                                  real64 const& dt,
//                                  DomainPartition * const domain,
//                                  systemSolverInterface::EpetraBlockSystem * const blockSystem);
//
//  virtual void AssembleSystem( DomainPartition * const domain,
//                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
//                               real64 const time,
//                               real64 const dt );
//
//  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
//                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
//                                        real64 const time,
//                                        real64 const dt );
//
//  virtual real64
//  CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
//                         DomainPartition * const domain );
//
//
//  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
//                            SystemSolverParameters const * const params );
//
//  virtual void
//  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
//                       real64 const scalingFactor,
//                       DomainPartition * const domain );
//
//  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain );
//
//
//  virtual void ImplicitStepComplete( real64 const & time,
//                                     real64 const & dt,
//                                     DomainPartition * const domain );

  virtual void ReadXML_PostProcess() override final;

  void UpdateDeformationForCoupling( DomainPartition * const domain );

  real64 SplitOperatorStep( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition * const domain);


  enum class couplingTypeOption
  {
    FixedStress,
    FullyImplicit
  };

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto couplingTypeOptionString = "couplingTypeOption";
    constexpr static auto biotCoefficienttring = "BiotCoe";

    constexpr static auto effectiveStressString = "effectiveStress";
    constexpr static auto deltaEffectiveStressString = "deltaEffectiveStress";
    constexpr static auto deltaVolumetricStrainString = "deltaVolumetricStrain";

    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto fluidSolverNameString = "fluidSolverName";

    constexpr static auto deltaPorosityString = "deltaPorosity";
    constexpr static auto dPorosity_dPressureString = "dPorosity_dPressure";
  } viewKeys;




private:
  string m_solidSolverName;
  string m_flowSolverName;
  couplingTypeOption m_couplingTypeOption;
  real64 m_biotCoef;

};

} /* namespace geosx */

#endif /* POROELASTICSOLVER_HPP_ */
