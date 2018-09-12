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
  * @file CompositionalMultiphaseFluid.cpp
  */

#include "CompositionalMultiphaseFluid.hpp"

#ifdef GEOSX_USE_PVT_PACKAGE
#include "MultiphaseSystem/CompositionalMultiphaseSystem.hpp"
#endif

using namespace PVTPackage;

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid(std::string const & name, ManagedGroup * const parent)
  : ConstitutiveBase(name, parent),
    m_fluid(nullptr)
{

}

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid()
{
  delete m_fluid;
}

std::unique_ptr<ConstitutiveBase>
CompositionalMultiphaseFluid::DeliverClone(string const & name, ManagedGroup * const parent) const
{
  auto clone = std::make_unique<CompositionalMultiphaseFluid>( name, parent );
  // TODO actually clone
#ifdef GEOSX_USE_PVT_PACKAGE
  //clone->m_fluid = new CompositionalMultiphaseSystem(...);
#endif
  return clone;
}

void CompositionalMultiphaseFluid::AllocateConstitutiveData(dataRepository::ManagedGroup * const parent,
                                                            localIndex const numConstitutivePointsPerParentIndex)
{
  // TODO
}

void CompositionalMultiphaseFluid::FillDocumentationNode()
{
  // TODO
}

void CompositionalMultiphaseFluid::ReadXML_PostProcess()
{
  // TODO read input

#ifdef GEOSX_USE_PVT_PACKAGE
  std::vector<std::string> Labels = { "N2","C10","C20","H2O" };
  auto Pc = { 34e5,25.3e5,14.6e5,220.5e5 };
  auto Tc = { 126.2,622.0,782.0,647.0 };
  auto Omega = { 0.04,0.443,0.816,0.344 };
  auto Mw = { 28e-3,134e-3,275e-3,18e-3 };
  auto nbc = Pc.size();
  const ComponentProperties CompProps(nbc, Labels, Mw, Tc, Pc, Omega);
  m_fluid = new CompositionalMultiphaseSystem({ PHASE_TYPE::OIL, PHASE_TYPE::GAS, PHASE_TYPE::LIQUID_WATER_RICH },
                                              { EOS_TYPE::PENG_ROBINSON, EOS_TYPE::PENG_ROBINSON, EOS_TYPE::PENG_ROBINSON },
                                              COMPOSITIONAL_FLASH_TYPE::TRIVIAL,
                                              CompProps);
#else
  GEOS_ERROR("Cannot use compositional fluid model without PVTPackage. Rebuild with ENABLE_PVT_PACKAGE=ON");
#endif
}

void CompositionalMultiphaseFluid::FinalInitialization(ManagedGroup * const parent)
{

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx