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

/*
 * TaperLayer.cpp
 *
 */

#include "TaperLayer.hpp"

namespace geos
{
using namespace dataRepository;

TaperLayer::TaperLayer( string const & name, Group * const parent ):
  FieldSpecificationBase( name, parent )
{
  registerWrapper( viewKeyStruct::xMinString(), &m_xMin ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { -LvArray::NumericLimits< real32 >::max,
                            -LvArray::NumericLimits< real32 >::max,
                            -LvArray::NumericLimits< real32 >::max } ).
    setDescription( "Minimum (x,y,z) coordinates of the inner PML boundaries" );

  registerWrapper( viewKeyStruct::xMaxString(), &m_xMax ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { LvArray::NumericLimits< real32 >::max,
                            LvArray::NumericLimits< real32 >::max,
                            LvArray::NumericLimits< real32 >::max } ).
    setDescription( "Maximum (x,y,z) coordinates of the inner PML boundaries" );

  registerWrapper( viewKeyStruct::taperCoeffString(), &m_taperCoeff ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.001 ).
    setDescription( "Desired coefficient of the taper region, used to compute the damping profile" );

  registerWrapper( viewKeyStruct::thicknessMinXYZString(), &m_thicknessMinXYZ ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {-1.0, -1.0, -1.0} ).
    setDescription( "Thickness of the taper region, at left, front, and top sides, used to compute the damping profile" );

  registerWrapper( viewKeyStruct::thicknessMaxXYZString(), &m_thicknessMaxXYZ ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {-1.0, -1.0, -1.0} ).
    setDescription( "Thickness of the taper region, at right, back, and bottom sides, used to compute the damping profile" );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::initialConditionString() ).
    setInputFlag( InputFlags::FALSE );
  initialCondition( false ); // to make sure this is not called by applyInitialConditions

}

void TaperLayer::postProcessInput()
{
  GEOS_THROW_IF( (m_xMax[0]<m_xMin[0] || m_xMax[1]<m_xMin[1] || m_xMax[2]<m_xMin[2]),
                 getCatalogName() << " " << getDataContext() << " "
                                  << viewKeyStruct::xMinString()
                                  << " must be smaller than "
                                  << viewKeyStruct::xMaxString(),
                 InputError );

  GEOS_THROW_IF( (m_taperCoeff<=0 || m_taperCoeff>1),
                 getCatalogName() << " " << getDataContext() << " "
                                  << viewKeyStruct::taperCoeffString()
                                  << " must satisfy 0 < taperCoeff <= 1",
                 InputError );

  GEOS_LOG_RANK_0_IF( (m_xMin[0]<smallestXMin || m_xMin[1]<smallestXMin || m_xMin[2]<smallestXMin),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::xMinString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_xMax[0]>largestXMax || m_xMax[1]>largestXMax || m_xMax[2]>largestXMax),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::xMaxString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_thicknessMinXYZ[0]<0 || m_thicknessMinXYZ[1]<0 || m_thicknessMinXYZ[2]<0),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::thicknessMinXYZString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_thicknessMaxXYZ[0]<0 || m_thicknessMaxXYZ[1]<0 || m_thicknessMaxXYZ[2]<0),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::thicknessMaxXYZString()
                                       << " will be computed internally" );

}


REGISTER_CATALOG_ENTRY( FieldSpecificationBase, TaperLayer, string const &, Group * const )

} /* namespace geos */
