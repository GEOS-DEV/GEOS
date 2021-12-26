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
 * @file SinglePhaseProppantBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_

#include "SinglePhaseBase.hpp"

namespace geosx
{

class SinglePhaseProppantBase : public SinglePhaseBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseProppantBase( const string & name,
                           Group * const parent );

  SinglePhaseProppantBase() = delete;

  /// deleted copy constructor
  SinglePhaseProppantBase( SinglePhaseProppantBase const & ) = delete;

  /// default move constructor
  SinglePhaseProppantBase( SinglePhaseProppantBase && ) = default;

  /// deleted assignment operator
  SinglePhaseProppantBase & operator=( SinglePhaseProppantBase const & ) = delete;

  /// deleted move operator
  SinglePhaseProppantBase & operator=( SinglePhaseProppantBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseProppantBase();

  virtual void updateFluidModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const override;

  virtual void updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion,
                                              localIndex const targetIndex ) const override;

protected:

  virtual void validateFluidModels( DomainPartition const & domain ) const override;

  virtual FluidPropViews getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const override;

private:


};
}
#endif /* GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_ */
