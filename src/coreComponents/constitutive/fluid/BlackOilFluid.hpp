/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BlackOilFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "constitutive/fluid/MultiFluidPVTPackageWrapper.hpp"

namespace geosx
{

namespace constitutive
{

class BlackOilFluid : public MultiFluidPVTPackageWrapper
{
public:

  using exec_policy = serialPolicy;

  BlackOilFluid( string const & name, Group * const parent );

  virtual ~BlackOilFluid() override;

  std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "BlackOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }


  struct viewKeyStruct : MultiFluidPVTPackageWrapper::viewKeyStruct
  {
    static constexpr char const * surfaceDensitiesString() { return "surfaceDensities"; }
    static constexpr char const * tableFilesString() { return "tableFiles"; }
  };

protected:
  virtual void postProcessInput() override;

private:

  void createFluid() override;

  // Black-oil phase/component description
  array1d< real64 > m_surfaceDensities;

  // Black-oil table filenames
  path_array m_tableFiles;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
