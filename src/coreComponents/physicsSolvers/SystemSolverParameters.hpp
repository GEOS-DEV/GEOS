/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SystemSolverParameters.hpp
 */

#ifndef GEOSX_SYSTEMSOLVERPARAMETERS_HPP_
#define GEOSX_SYSTEMSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{

class SystemSolverParameters : public dataRepository::Group
{
public:
  SystemSolverParameters() = delete;

  SystemSolverParameters( std::string const & name,
                          Group * const parent );

  SystemSolverParameters(SystemSolverParameters &&) = default;

  virtual ~SystemSolverParameters() override = default;

  static string CatalogName() { return "SystemSolverParameters"; }

  virtual void PostProcessInput() override;

  struct viewKeysStruct
  {
    static constexpr auto solverTypeString          = "solverType";
    static constexpr auto krylovTolString           = "krylovTol";
    static constexpr auto useAdaptiveKrylovString   = "useAdaptiveKrylovTol";
    static constexpr auto numKrylovIterString       = "numKrylovIter";
    static constexpr auto kspaceString              = "kspace";
    static constexpr auto ilut_fillString           = "ilut_fill";
    static constexpr auto ilut_dropString           = "ilut_drop";
    static constexpr auto useMLPrecondString        = "useMLPrecond";
    static constexpr auto useInnerSolverString      = "useInnerSolver";
    static constexpr auto scalingOptionString       = "scalingOption";
    static constexpr auto useBicgstabString         = "useBicgstab";
    static constexpr auto useDirectSolverString     = "useDirectSolver";
    static constexpr auto krylovResidualInitString  = "krylovResidualInit";
    static constexpr auto krylovResidualFinalString = "krylovResidualFinal";
  } viewKeys;

  struct groupKeysStruct
  {} groupKeys;

  string  solverType() const                  { return m_solverType; }
  real64 krylovTol() const                    { return m_krylovTol; }
  integer useAdaptiveKrylovTol() const        { return m_useAdaptiveKrylovTol; }
  integer  numKrylovIter() const              { return m_numKrylovIter; }
  integer  kspace() const                     { return m_kspace; }
  real64 ilut_fill() const                    { return m_ilut_fill; }
  real64 ilut_drop() const                    { return m_ilut_drop; }
  integer   useMLPrecond() const              { return m_useMLPrecond; }
  integer   useInnerSolver() const            { return m_useInnerSolver; }
  integer  scalingOption() const              { return m_scalingOption; }
  integer   useBicgstab() const               { return m_useBicgstab; }
  integer   useDirectSolver() const           { return m_useDirectSolver; }
  real64 krylovResidualInit() const           { return m_krylovResidualInit; }
  real64 krylovResidualFinal() const          { return m_krylovResidualFinal; }

  string  m_solverType;
  real64  m_krylovTol;
  integer m_useAdaptiveKrylovTol;
  integer m_numKrylovIter;
  integer m_kspace;
  real64  m_ilut_fill;
  real64  m_ilut_drop;
  integer m_useMLPrecond;
  integer m_useInnerSolver;
  integer m_scalingOption;
  integer m_useBicgstab;
  integer m_useDirectSolver;
  real64  m_krylovResidualInit;
  real64  m_krylovResidualFinal;
  real64  m_krylovAuxTime;
  real64  m_krylovSetupTime;
  real64  m_krylovSolveTime;
  integer m_maxIters = 1000;


};

} /* namespace geosx */

#endif /*GEOSX_SYSTEMSOLVERPARAMETERS_HPP_*/
