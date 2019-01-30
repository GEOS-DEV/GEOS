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

    static constexpr auto pressureString      = "pressure";
    static constexpr auto deltaPressureString = "deltaPressure";
    
    static constexpr auto flowRateString = "flowRate";
    static constexpr auto bhpString      = "bhp";

    static constexpr auto avgDensityString = "avgDensity";

    static constexpr auto phaseComponentFracString        = "phaseComponentFraction";
    static constexpr auto dPhaseComponentFrac_dPresString = "dPhaseComponentFraction_dPres";
    static constexpr auto dPhaseComponentFrac_dCompString = "dPhaseComponentFraction_dComp";
    static constexpr auto phaseDensityString              = "phaseDensity";
    static constexpr auto dPhaseDensity_dPresString       = "dPhaseDensity_dPres";
    static constexpr auto dPhaseDensity_dCompString       = "dPhaseDensity_dComp";
    static constexpr auto phaseViscosityString            = "phaseViscosity";
    static constexpr auto dPhaseViscosity_dPresString     = "dPhaseViscosity_dPres";
    static constexpr auto dPhaseViscosity_dCompString     = "dPhaseViscosity_dComp";
    static constexpr auto phaseMobilityString             = "phaseMobility";
    static constexpr auto dPhaseMobility_dPresString      = "dPhaseMobility_dPres";
    static constexpr auto dPhaseMobility_dCompString      = "dPhaseMobility_dComp";

    using ViewKey = dataRepository::ViewKey;
    
    // primary solution field
    ViewKey pressure      = { pressureString };
    ViewKey deltaPressure = { deltaPressureString };

    // well controls
    ViewKey flowRate = { flowRateString };
    ViewKey bhp      = { bhpString };

    // average density to compute hydrostatic head
    ViewKey avgDensity = { avgDensityString };

    // intermediate values used to compute perforation rates (+ derivatives w.r.t well and res vars)
    ViewKey phaseComponentFrac        = { phaseComponentFracString };
    ViewKey dPhaseComponentFrac_dPres = { dPhaseComponentFrac_dPresString };
    ViewKey dPhaseComponentFrac_dComp = { dPhaseComponentFrac_dCompString };
    ViewKey phaseDensity              = { phaseDensityString };
    ViewKey dPhaseDensity_dPres       = { dPhaseDensity_dPresString };
    ViewKey dPhaseDensity_dComp       = { dPhaseDensity_dCompString };
    ViewKey phaseViscosity            = { phaseViscosityString };
    ViewKey dPhaseViscosity_dPres     = { dPhaseViscosity_dPresString };
    ViewKey dPhaseViscosity_dComp     = { dPhaseViscosity_dCompString };
    ViewKey phaseMobility             = { phaseMobilityString };
    ViewKey dPhaseMobility_dPres      = { dPhaseMobility_dPresString };
    ViewKey dPhaseMobility_dComp      = { dPhaseMobility_dCompString };
    
  } viewKeysCompositionalMultiphaseWell;

  struct groupKeyStruct : public WellBase::groupKeyStruct
  {

  } groupKeysCompositionalMultiphaseWell;

private:

  // update each connection pressure from bhp and hydrostatic head
  void StateUpdate( DomainPartition const * domain, localIndex fluidIndex );

  // form the well control equation based on the type of well
  void FormControlEquation( DomainPartition const * const domain,
                            Epetra_FECrsMatrix * const jacobian,
                            Epetra_FEVector * const residual,
                            real64 const time_n,
                            real64 const dt );
  
  array1d<real64> m_bhp;

};

} //namespace geosx


#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP
