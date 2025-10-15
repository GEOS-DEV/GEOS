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

/**
 * @file SinglePhaseStressPathDrivenFVM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTRESSPATHDRIVENFVM_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTRESSPATHDRIVENFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "prescribedStressPathGeomechanics/OedometricStressPath.hpp"

namespace geos
{


/**
 * @class SinglePhaseStressPathDrivenFVM
 *
 * class to perform a single phase finite volume solve
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
class SinglePhaseStressPathDrivenFVM : public SinglePhaseFVM< SinglePhaseBase >
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseStressPathDrivenFVM( const string & name,
                                   dataRepository::Group * const parent );


  /// deleted default constructor
  SinglePhaseStressPathDrivenFVM() = delete;

  /// deleted copy constructor
  SinglePhaseStressPathDrivenFVM( SinglePhaseStressPathDrivenFVM const & ) = delete;

  /// default move constructor
  SinglePhaseStressPathDrivenFVM( SinglePhaseStressPathDrivenFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseStressPathDrivenFVM & operator=( SinglePhaseStressPathDrivenFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseStressPathDrivenFVM & operator=( SinglePhaseStressPathDrivenFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseStressPathDrivenFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SinglePhaseStressPathDrivenFVM";
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : SinglePhaseBase::viewKeyStruct
  {
    /// Name of the hydraulicApertureRelationName
    static constexpr char const * hydraulicApertureRelationNameString() { return "hydraulicApertureRelationName"; }
  };
  struct groupKeyStruct //: PhysicsSolverBase::groupKeyStruct
  {
    /// @return string for the PrescribedStressPathGeomechanicsBase wrapper
    static constexpr char const * oedometricStressPathString() { return "OedometricStressPath"; }
    static constexpr char const * prescribedStressPathGeomechanicsBaseString() { return "PrescribedStressPathGeomechanicsBase"; }
  };

  // The same as in PoromechanicsSolver.hpp
  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;

  virtual void updateState ( DomainPartition & domain ) override final;


private:

  void updateHydraulicAperture( DomainPartition & domain ) const;
  
  void updateFractureAperture( SurfaceElementSubRegion & subRegion ) const;

  /// Just like 'Nonlinear solver parameters' in PhysicsSolverBase.hpp
  OedometricStressPath m_oedometricStressPath;

  real64 m_apertureInSitu;

};

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTRESSPATHDRIVENFVM_HPP_
