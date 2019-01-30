/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * @file CompositionalMultiphaseWell.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP

#include "WellBase.hpp"

class Epetra_FECrsGraph;
class Epetra_FECrsMatrix;
class Epetra_FEVector;

namespace geosx
{

namespace systemSolverInterface
{
class EpetraBlockSystem;
class LinearSolverWrapper;
enum class BlockIDs;
}
  
namespace dataRepository
{
namespace keys
{
static constexpr auto compositionalMultiphaseWell = "CompositionalMultiphaseWell";
}
}

class DomainPartition;

class CompositionalMultiphaseWell : public WellBase
{
public:
  explicit CompositionalMultiphaseWell( string const & name, dataRepository::ManagedGroup * const parent );
  ~CompositionalMultiphaseWell() override;

  CompositionalMultiphaseWell() = delete;
  CompositionalMultiphaseWell( CompositionalMultiphaseWell const &) = delete;
  CompositionalMultiphaseWell( CompositionalMultiphaseWell && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return dataRepository::keys::compositionalMultiphaseWell; }

  virtual void InitializePostSubGroups( ManagedGroup * const rootGroup ) override;

  // Initialize all variables in the well
  void InitializeState( DomainPartition * const domain );

  // Apply well solution after the linear solve
  void ApplySolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                      real64 const scalingFactor,
                      DomainPartition * const domain );
  
  // Compute the rate for all the perforations and form the control equation
  void AssembleWellTerms( DomainPartition * const domain,
                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                          real64 const time_n,
                          real64 const dt );

  // Check if the well control needs to be switched
  void CheckControlSwitch();

  // Reset the well to its initial state at the beginning of the time step
  void ResetStateToBeginningOfStep( DomainPartition * const domain );

  // Report the total flow rate
  real64 GetTotalFlowRate();

  struct viewKeyStruct : public WellBase::viewKeyStruct
  {

    static constexpr auto pressureString               = "pressure";
    static constexpr auto deltaPressureString          = "deltaPressure";
    static constexpr auto globalCompDensityString      = "globalCompDensity";
    static constexpr auto deltaGlobalCompDensityString = "deltaGlobalCompDensity";
    static constexpr auto mixtureVelocityString        = "mixtureVelocity";
    static constexpr auto deltaMixtureVelocityString   = "deltaMixtureVelocity";

    static constexpr auto flowRateString = "flowRate";
    static constexpr auto bhpString      = "bhp";

    static constexpr auto globalComponentFracString        = "globalComponentFraction";
    static constexpr auto dGlobalComponentFrac_dPresString = "dGlobalComponentFraction_dPres";
    static constexpr auto dGlobalComponentFrac_dCompString = "dGlobalComponentFraction_dComp";
    

    using ViewKey = dataRepository::ViewKey;
    
    // primary solution field
    ViewKey pressure               = { pressureString };
    ViewKey deltaPressure          = { deltaPressureString };
    ViewKey globalCompDensity      = { globalCompDensityString };
    ViewKey deltaGlobalCompDensity = { deltaGlobalCompDensityString };
    ViewKey mixtureVelocity        = { mixtureVelocityString };
    ViewKey deltaMixtureVelovity   = { deltaMixtureVelocityString };
    
    // well controls
    ViewKey flowRate = { flowRateString };
    ViewKey bhp      = { bhpString };

    // global composition to input injection stream
    ViewKey globalComponentFrac        = { globalComponentFracString };
    ViewKey dGlobalComponentFrac_dPres = { dGlobalComponentFrac_dPresString };
    ViewKey dGlobalComponentFrac_dComp = { dGlobalComponentFrac_dCompString };    
    
  } viewKeysCompositionalMultiphaseWell;

  struct groupKeyStruct : public WellBase::groupKeyStruct
  {

  } groupKeysCompositionalMultiphaseWell;

private:

  // update each connection pressure from bhp and hydrostatic head
  void StateUpdate( DomainPartition const * domain, localIndex fluidIndex );

  // form the well control equation based on the type of well
  void FormControlEquation( DomainPartition const * const domain,
                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            real64 const time_n,
                            real64 const dt );
  
  array1d<real64> m_bhp;

  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;
  
};

} //namespace geosx


#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP
