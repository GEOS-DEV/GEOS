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
  virtual ~SystemSolverParameters();

  void FillDocumentationNode( dataRepository::ManagedGroup * const  );

  struct keyStruct
  {
    dataRepository::ViewKey verbosity           = { "verbosityFlag" };
    dataRepository::ViewKey krylovTol           = { "krylovTol" };
    dataRepository::ViewKey numKrylovIter       = { "numKrylovIter" };
    dataRepository::ViewKey kspace              = { "kspace" };
    dataRepository::ViewKey ilut_fill           = { "ilut_fill" };
    dataRepository::ViewKey ilut_drop           = { "ilut_drop" };
    dataRepository::ViewKey useMLPrecond        = { "useMLPrecond" };
    dataRepository::ViewKey useInnerSolver      = { "useInnerSolver" };
    dataRepository::ViewKey scalingOption       = { "scalingOption" };
    dataRepository::ViewKey useBicgstab         = { "useBicgstab" };
    dataRepository::ViewKey useDirectSolver     = { "useDirectSolver" };
    dataRepository::ViewKey KrylovResidualInit  = { "KrylovResidualInit" };
    dataRepository::ViewKey KrylovResidualFinal = { "KrylovResidualFinal" };
    dataRepository::ViewKey useNewtonSolve      = { "useNewtonSolve" };
    dataRepository::ViewKey newtonTol           = { "newtonTol" };
    dataRepository::ViewKey maxIterNewton       = { "maxIterNewton" };
  }keys;

  int32  verbose() const              { return *(this->getData<int32>( keys.verbosity )); }
  real64 krylovTol() const            { return *(this->getData<int32>( keys.krylovTol )); }
  int32  numKrylovIter() const        { return *(this->getData<int32>( keys.numKrylovIter )); }
  int32  kspace() const               { return *(this->getData<int32>( keys.kspace )); }
  real64 ilut_fill() const            { return *(this->getData<int32>( keys.ilut_fill )); }
  real64 ilut_drop() const            { return *(this->getData<int32>( keys.ilut_drop )); }
  bool   useMLPrecond() const         { return *(this->getData<int32>( keys.useMLPrecond )); }
  bool   useInnerSolver() const       { return *(this->getData<int32>( keys.useInnerSolver )); }
  int32  scalingOption() const        { return *(this->getData<int32>( keys.scalingOption )); }
  bool   useBicgstab() const          { return *(this->getData<int32>( keys.useBicgstab )); }
  bool   useDirectSolver() const      { return *(this->getData<int32>( keys.useDirectSolver )); }
  real64 KrylovResidualInit() const   { return *(this->getData<int32>( keys.KrylovResidualInit )); }
  real64 KrylovResidualFinal() const  { return *(this->getData<int32>( keys.KrylovResidualFinal )); }
  bool   useNewtonSolve() const       { return *(this->getData<int32>( keys.useNewtonSolve )); }
  real64 newtonTol() const            { return *(this->getData<int32>( keys.newtonTol )); }
  int32  maxIterNewton() const        { return *(this->getData<int32>( keys.maxIterNewton )); }


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_SYSTEMSOLVERPARAMETERS_HPP_ */
