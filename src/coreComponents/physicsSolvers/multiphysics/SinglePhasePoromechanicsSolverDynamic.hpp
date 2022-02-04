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
 * @file SinglePhasePoromechanicsSolverDynamic.hpp
 *
 *  Created on: Jan 29, 2022
 *      Author: ron
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSSOLVERDYNAMIC_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSSOLVERDYNAMIC_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsSolver.hpp"

namespace geosx {

class SinglePhasePoromechanicsSolverDynamic : public SinglePhasePoromechanicsSolver
{
public:
	SinglePhasePoromechanicsSolverDynamic( const string & name,
										   Group * const parent );

	~SinglePhasePoromechanicsSolverDynamic() override;

	  /**
	   * @brief name of the node manager in the object catalog
	   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
	   */
	  static string catalogName()
	  {
	    return "SinglePhasePoromechanicsDynamic";
	  }


	  virtual real64 solverStep( real64 const & time_n,
	                             real64 const & dt,
	                             int const cycleNumber,
	                             DomainPartition & domain ) override;

	  virtual void
	  explicitStepSetup( real64 const & time_n,
	                     real64 const & dt,
	                     DomainPartition & domain) override final;

	  virtual real64 explicitStep( real64 const & time_n,
	                               real64 const & dt,
	                               integer const cycleNumber,
	                               DomainPartition & domain ) override;

	  void updateDeformationForCoupling( DomainPartition & domain );

	  void applyPressureToFacesInExplicitSolver( DomainPartition & domain );

	  enum class CouplingTypeOption : integer
	  {
	    FIM,
		FEM_ExplicitTransient,
		FEM_ImplicitTransient
	  };

	  struct viewKeyStruct : SolverBase::viewKeyStruct
	  {
	    constexpr static char const * couplingTypeOptionString() { return "couplingTypeOptionEnum"; }

	    constexpr static char const * couplingTypeOptionStringString() { return "couplingTypeOption"; }

	  };

protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  CouplingTypeOption m_couplingTypeOption;

};

ENUM_STRINGS( SinglePhasePoromechanicsSolverDynamic::CouplingTypeOption,
              "FIM",
              "FEM_ExplicitTransient",
			  "FEM_ImplicitTransient" );


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSSOLVERDYNAMIC_HPP_ */
