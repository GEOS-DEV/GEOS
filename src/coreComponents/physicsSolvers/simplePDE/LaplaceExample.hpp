#pragma once

#include "physicsSolvers/PhysicsSolverBase.hpp"

namespace geos
{

class LaplaceExample final : public PhysicsSolverBase 
{
public:
  LaplaceExample() = delete;

  LaplaceExample( const string & name,
                  Group * const parent );

  virtual ~LaplaceExample() override = default;


  static string catalogName() { return "LaplaceExample"; }

  string getCatalogName() const override { return catalogName(); }


  virtual void 
  registerDataOnMesh( Group & meshBodies ) override;


  virtual real64 
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  // virtual void 
  // setSparsityPattern( DomainPartition & domain,
  //                     DofManager & dofManager,
  //                     CRSMatrix< real64, globalIndex > & localMatrix,
  //                     SparsityPattern< globalIndex > & pattern ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;


  struct viewKeyStruct : public PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * fieldVarName() { return "fieldName"; }
  };

  virtual void updateState( DomainPartition &  ) override {}

  virtual void implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            DomainPartition & GEOS_UNUSED_PARAM( domain ) ) override
  {}

private:
  string m_fieldName;


};


} // namespace geos
