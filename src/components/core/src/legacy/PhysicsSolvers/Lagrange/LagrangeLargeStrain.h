/*
 * LagrangeLargeStrain.h
 *
 *  Created on: Nov 8, 2012
 *      Author: settgast1
 */

#ifndef LAGRANGELARGESTRAIN_H_
#define LAGRANGELARGESTRAIN_H_

#include "LagrangeSolverBase.h"

class LagrangeLargeStrain : public LagrangeSolverBase
{
public:
  LagrangeLargeStrain(  const std::string& name,
                        ProblemManagerT* const pm );
  virtual ~LagrangeLargeStrain();


  static const char* SolverName()
  {
    return "LagrangeLargeStrain";
  };

  void ApplyThermalStress( ElementRegionT& elemRegion,
                           NodeManager& nodeManager,
                           const localIndex& elementID,
                           Epetra_SerialDenseVector * rhs);

private:
  LagrangeLargeStrain();


  void ProcessElementRegion( NodeManager& nodeManager,
                             ElementRegionT& elemRegion,
                             const realT dt );


  void CalculateNodalForceFromStress(NodeManager& nodeManager,
                                     ElementRegionT& elemRegion,
                                     const localIndex& elementID,
                                     R2SymTensor& stress );

  virtual realT CalculateElementResidualAndDerivative( const MaterialBaseParameterData& matParams,
                                                       const FiniteElementBase& fe,
                                                       const Array2dT<R1Tensor>& dNdX,
                                                       const realT* const detJ,
                                                       R2SymTensor const * const refStress,
                                                       array<R1Tensor> const & u,
                                                       array<R1Tensor> const & uhat,
                                                       array<R1Tensor> const & uhattilde,
                                                       array<R1Tensor> const & vtilde,
                                                       realT const dt,
                                                       Epetra_SerialDenseMatrix& dRdU,
                                                       Epetra_SerialDenseVector& R ){ return 0;}



};

#endif /* LAGRANGELARGESTRAIN_H_ */
