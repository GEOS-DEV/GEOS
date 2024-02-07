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

/**
 * @file MatrixFreeSolidMechanicsFEM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MFSOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MFSOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_

#include "dataRepository/ExecutableGroup.hpp"
#include "physicsSolvers/simplePDE/LaplaceBaseH1.hpp"  // a base class shared by all Laplace solvers
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/common/LinearOperatorWithBC.hpp"

namespace geos
{



//START_SPHINX_INCLUDE_BEGINCLASS
class MatrixFreeSolidMechanicsFEM : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  MatrixFreeSolidMechanicsFEM() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  MatrixFreeSolidMechanicsFEM( const string & name,
                               Group * const parent );

  /// Destructor
  virtual ~MatrixFreeSolidMechanicsFEM() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "MatrixFreeSolidMechanicsFEM"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual void setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                          DofManager & dofManager ) const override;

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  void applyTractionBC( real64 const time,
                        DofManager const & dofManager,
                        DomainPartition & domain,
                        arrayView1d< real64 > const & localRhs );


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * kernelOptimizationOption() { return "kernelOptimizationOption"; }
  } solidMechanicsViewKeys;

protected:
  string m_fieldName;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;

  using FieldType = array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM >;

  int m_kernelOptimizationOption = 0;


};
} /* namespace geos */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MFSOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_ */
