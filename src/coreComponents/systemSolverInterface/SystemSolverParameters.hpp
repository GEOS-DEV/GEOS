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

/*
 * SystemSolverParameters.hpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_SYSTEMSOLVERPARAMETERS_HPP_
#define SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_SYSTEMSOLVERPARAMETERS_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class SystemSolverParameters : public dataRepository::ManagedGroup
{
public:
  SystemSolverParameters() = delete;

  SystemSolverParameters( std::string const & name,
                          ManagedGroup * const parent );

  SystemSolverParameters(SystemSolverParameters &&) = default;

  virtual ~SystemSolverParameters() override = default;

  static string CatalogName() { return "SystemSolverParameters"; }

  struct viewKeysStruct
  {
    static constexpr auto verbosityString           = "verbosityFlag";
    static constexpr auto solverTypeString          = "solverType";
    static constexpr auto krylovTolString           = "krylovTol";
    static constexpr auto numKrylovIterString       = "numKrylovIter";
    static constexpr auto kspaceString              = "kspace";
    static constexpr auto ilut_fillString           = "ilut_fill";
    static constexpr auto ilut_dropString           = "ilut_drop";
    static constexpr auto useMLPrecondString        = "useMLPrecond";
    static constexpr auto useInnerSolverString      = "useInnerSolver";
    static constexpr auto scalingOptionString       = "scalingOption";
    static constexpr auto useBicgstabString         = "useBicgstab";
    static constexpr auto useDirectSolverString     = "useDirectSolver";
    static constexpr auto KrylovResidualInitString  = "KrylovResidualInit";
    static constexpr auto KrylovResidualFinalString = "KrylovResidualFinal";
    static constexpr auto useNewtonSolveString      = "useNewtonSolve";
    static constexpr auto newtonTolString           = "newtonTol";
    static constexpr auto maxIterNewtonString       = "maxIterNewton";
    static constexpr auto numNewtonIterationsString = "numberOfNewtonIterations";
    static constexpr auto maxTimeStepCutsString     = "maxTimeStepCuts";
    static constexpr auto timeStepCutFactorString   = "timestepCutFactor";
    static constexpr auto maxLineSearchCutsString   = "maxLineSearchCuts";
    static constexpr auto lineSearchCutFactorString = "lineSearchCutFactor";
    static constexpr auto allowNonConvergedString   = "allowNonConverged";
    static constexpr auto doLineSearchString   = "doLineSearch";    

  } viewKeys;

  struct groupKeysStruct
  {} groupKeys;

  integer  verbose() const                    { return m_verbose; }
  string  solverType() const                  { return m_solverType; }
  real64 krylovTol() const                    { return m_krylovTol; }
  integer  numKrylovIter() const              { return m_numKrylovIter; }
  integer  kspace() const                     { return m_kspace; }
  real64 ilut_fill() const                    { return m_ilut_fill; }
  real64 ilut_drop() const                    { return m_ilut_drop; }
  integer   useMLPrecond() const              { return m_useMLPrecond; }
  integer   useInnerSolver() const            { return m_useInnerSolver; }
  integer  scalingOption() const              { return m_scalingOption; }
  integer   useBicgstab() const               { return m_useBicgstab; }
  integer   useDirectSolver() const           { return m_useDirectSolver; }
  real64 KrylovResidualInit() const           { return m_KrylovResidualInit; }
  real64 KrylovResidualFinal() const          { return m_KrylovResidualFinal; }
  integer   useNewtonSolve() const            { return m_useNewtonSolve; }
  real64 newtonTol() const                    { return m_newtonTol; }
  integer  maxIterNewton() const              { return m_maxIterNewton; }
  integer const & numNewtonIterations() const { return m_numNewtonIterations; }
  integer & numNewtonIterations()             { return m_numNewtonIterations; }

  integer maxTimeStepCuts() const             { return m_maxTimeStepCuts; }
  real64  timeStepCutFactor() const           { return m_timeStepCutFactor; }
  integer maxLineSearchCuts() const           { return m_maxLineSearchCuts; }
  real64  lineSearchCutFactor() const         { return m_lineSearchCutFactor; }
  integer allowNonConverged() const           { return m_allowNonConverged; }

  integer doLineSearch() const           { return m_doLineSearch; }  



  integer m_verbose;
  string  m_solverType;
  real64  m_krylovTol;
  integer m_numKrylovIter;
  integer m_kspace;
  real64  m_ilut_fill;
  real64  m_ilut_drop;
  integer m_useMLPrecond;
  integer m_useInnerSolver;
  integer m_scalingOption;
  integer m_useBicgstab;
  integer m_useDirectSolver;
  real64  m_KrylovResidualInit;
  real64  m_KrylovResidualFinal;
  integer m_useNewtonSolve;
  real64  m_newtonTol;
  integer m_maxIterNewton;
  integer m_numNewtonIterations;

  integer m_maxTimeStepCuts;
  real64  m_timeStepCutFactor;
  integer m_maxLineSearchCuts;
  real64  m_lineSearchCutFactor;
  integer m_allowNonConverged;

  integer m_doLineSearch;  


};

} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_SYSTEMSOLVERPARAMETERS_HPP_
        */
