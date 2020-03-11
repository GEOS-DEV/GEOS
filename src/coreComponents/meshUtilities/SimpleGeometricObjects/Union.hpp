/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file Union.hpp
 *
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_UNION_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_UNION_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class Union : public SimpleGeometricObjectBase
{
public:
  Union( const std::string& name,
       Group * const parent );

  virtual ~Union() override;

  static string CatalogName() { return "Union"; }

  bool IsCoordInObject( const R1Tensor& coord ) const override final;


private:

  string_array m_objects;
  
  struct viewKeyStruct
  {
    static constexpr auto objects = "objects";
  } viewKeys;


};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_UNION_HPP_
        */
