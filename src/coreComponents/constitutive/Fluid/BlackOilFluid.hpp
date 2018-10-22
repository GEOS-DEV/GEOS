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

/**
  * @file BlackOilFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BLACKOILFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BLACKOILFLUID_HPP_

#include "constitutive/Fluid/MultiFluidBase.hpp"

namespace PVTPackage
{
class BlackOilMultiphaseSystem;
}

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const blackOilFluid = "BlackOilFluid";
}
}

namespace constitutive
{

class BlackOilFluid : public MultiFluidBase
{
public:

  BlackOilFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~BlackOilFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::blackOilFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override final {}

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

  virtual void StateUpdatePointMultiphaseFluid(real64 const & pres,
                                               real64 const & temp,
                                               real64 const * composition,
                                               localIndex const k,
                                               localIndex const q) override;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr auto surfaceDensitiesString = "surfaceDensities";
    static constexpr auto molarWeightsString = "molarWeights";
    static constexpr auto tableFilesString = "tableFiles";
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey surfaceDensities = { surfaceDensitiesString };
    ViewKey molarWeights     = { molarWeightsString };
    ViewKey tableFiles       = { tableFilesString };

  } viewKeys;

private:

  void createFluid();

  // Black-oil phase/component description
  array1d<real64> m_surfaceDensities;
  array1d<real64> m_molarWeights;

  // Black-oil table filenames
  string_array m_tableFiles;

  // PVTPackage fluid object
  PVTPackage::BlackOilMultiphaseSystem * m_fluid;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BLACKOILFLUID_HPP_
