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

  struct viewKeyStruct
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
  }viewKeys;

  int32  verbose() const              { return *(this->getData<int32>( viewKeys.verbosity )); }
  real64 krylovTol() const            { return *(this->getData<real64>( viewKeys.krylovTol )); }
  int32  numKrylovIter() const        { return *(this->getData<int32>( viewKeys.numKrylovIter )); }
  int32  kspace() const               { return *(this->getData<int32>( viewKeys.kspace )); }
  real64 ilut_fill() const            { return *(this->getData<real64>( viewKeys.ilut_fill )); }
  real64 ilut_drop() const            { return *(this->getData<int32>( viewKeys.ilut_drop )); }
  int32   useMLPrecond() const         { return *(this->getData<int32>( viewKeys.useMLPrecond )); }
  int32   useInnerSolver() const       { return *(this->getData<int32>( viewKeys.useInnerSolver )); }
  int32  scalingOption() const        { return *(this->getData<int32>( viewKeys.scalingOption )); }
  int32   useBicgstab() const          { return *(this->getData<int32>( viewKeys.useBicgstab )); }
  int32   useDirectSolver() const      { return *(this->getData<int32>( viewKeys.useDirectSolver )); }
  real64 KrylovResidualInit() const   { return *(this->getData<real64>( viewKeys.KrylovResidualInit )); }
  real64 KrylovResidualFinal() const  { return *(this->getData<real64>( viewKeys.KrylovResidualFinal )); }
  int32   useNewtonSolve() const       { return *(this->getData<int32>( viewKeys.useNewtonSolve )); }
  real64 newtonTol() const            { return *(this->getData<real64>( viewKeys.newtonTol )); }
  int32  maxIterNewton() const        { return *(this->getData<int32>( viewKeys.maxIterNewton )); }


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_SYSTEMSOLVERPARAMETERS_HPP_ */
