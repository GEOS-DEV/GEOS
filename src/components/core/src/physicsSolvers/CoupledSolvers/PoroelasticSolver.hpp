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

  virtual real64 SolverStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::ManagedGroup * domain ) override;

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

  real64 SplitOperatorStep( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition * const domain);


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto fluidSolverNameString = "fluidSolverName";
  } viewKeys;




private:
  string m_solidSolverName;
  string m_flowSolverName;

};

} /* namespace geosx */

#endif /* POROELASTICSOLVER_HPP_ */
