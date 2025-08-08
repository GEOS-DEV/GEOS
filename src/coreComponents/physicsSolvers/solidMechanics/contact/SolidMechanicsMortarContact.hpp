/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * SolidMechanicsMortarContact.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_

#include "physicsSolvers/solidMechanics/contact/ContactSolverBase.hpp"

namespace geos
{

class SolidMechanicsMortarContact : public ContactSolverBase
{
public:
  SolidMechanicsMortarContact( const string & name,
                                            Group * const parent );

  ~SolidMechanicsMortarContact() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsMortarContact";
  }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
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
                            bool const setSparsity = true ) override final;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override final;

  virtual void implicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64 calculateResidualNorm( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs ) override;

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;

private:

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_ */
