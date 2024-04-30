/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * SolidMechanicsEmbeddedFractures.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_

#include "physicsSolvers/contact/ContactSolverBase.hpp"

namespace geos
{
using namespace constitutive;

class SolidMechanicsEmbeddedFractures : public ContactSolverBase
{
public:
  SolidMechanicsEmbeddedFractures( const string & name,
                                   Group * const parent );

  ~SolidMechanicsEmbeddedFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsEmbeddedFractures";
  }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void implicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override final;

  void updateState( DomainPartition & domain ) override final;

  void addCouplingNumNonzeros( DomainPartition & domain,
                               DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addCouplingSparsityPattern( DomainPartition const & domain,
                                   DofManager const & dofManager,
                                   SparsityPatternView< globalIndex > const & pattern ) const;

  void applyTractionBC( real64 const time_n,
                        real64 const dt,
                        DomainPartition & domain );

  virtual bool updateConfiguration( DomainPartition & domain ) override final;

  bool useStaticCondensation() const { return m_useStaticCondensation; }

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {
    constexpr static char const * useStaticCondensationString() { return "useStaticCondensation"; }
  };

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final;

  virtual void postProcessInput() override final;

private:

  void updateJump( DofManager const & dofManager,
                   real64 const dt,
                   DomainPartition & domain );

  /// decide whether to use static condensation or not
  integer m_useStaticCondensation;

};


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_ */
