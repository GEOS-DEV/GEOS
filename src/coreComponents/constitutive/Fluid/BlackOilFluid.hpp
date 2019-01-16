/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "constitutive/Fluid/MultiFluidPVTPackageWrapper.hpp"

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

class BlackOilFluid : public MultiFluidPVTPackageWrapper
{
public:

  BlackOilFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~BlackOilFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::blackOilFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }


  struct viewKeyStruct : MultiFluidPVTPackageWrapper::viewKeyStruct
  {
    static constexpr auto surfaceDensitiesString = "surfaceDensities";
    static constexpr auto tableFilesString = "tableFiles";
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey surfaceDensities = { surfaceDensitiesString };
    ViewKey tableFiles       = { tableFilesString };

  } viewKeysBlackOilFluid;

protected:
  virtual void PostProcessInput() override;

private:

  void createFluid() override;

  // Black-oil phase/component description
  array1d<real64> m_surfaceDensities;

  // Black-oil table filenames
  string_array m_tableFiles;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BLACKOILFLUID_HPP_
