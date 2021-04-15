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
 * @file ContactRelationBase.cpp
 */

#include "ContactRelationBase.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"



namespace geosx
{

using namespace dataRepository;


namespace constitutive
{


ContactRelationBase::ContactRelationBase( string const & name,
                                          Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_penaltyStiffness( 0.0 ),
  m_apertureFunction( nullptr ),
  m_apertureTolerance( 1.0e-99 )
{
  registerWrapper( viewKeyStruct::penaltyStiffnessString(), &m_penaltyStiffness ).
    //setInputFlag( InputFlags::REQUIRED ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the penetration penalty stiffness. Units of Pressure/length" );

  registerWrapper( viewKeyStruct::apertureToleranceString(), &m_apertureTolerance ).
    setApplyDefaultValue( 1.0e-9 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to be used to avoid floating point errors in expressions involving aperture. "
                    "For example in the case of dividing by the actual aperture (not the effective aperture "
                    "that results from the aperture function) this value may be used to avoid 1/0 errors. "
                    "Note that this value may have some physical significance in its usage, as it may be used "
                    "to smooth out highly nonlinear behavior associated with 1/0 in addition to avoiding the "
                    "1/0 error." );

}

ContactRelationBase::~ContactRelationBase()
{}

real64 ContactRelationBase::limitTangentialTractionNorm( real64 const GEOSX_UNUSED_PARAM( normalTraction ) ) const
{
  GEOSX_ERROR( "ContactRelationBase::limitTangentialTractionNorm called!. Should be overridden." );
  return 0;
}

real64 ContactRelationBase::dLimitTangentialTractionNorm_dNormalTraction( real64 const GEOSX_UNUSED_PARAM( normalTraction ) ) const
{
  GEOSX_ERROR( "ContactRelationBase::dLimitTangentialTractionNorm_dNormalTraction called!. Should be overridden." );
  return 0;
}

Group *
ContactRelationBase::createChild( string const & catalogKey, string const & childName )
{
  FunctionBase::CatalogInterface::CatalogType const & functionCatalog = FunctionBase::getCatalog();
  GEOSX_ERROR_IF( !functionCatalog.count( catalogKey ), catalogKey << " is an invalid key ContactRelationBase child group." );

  m_apertureFunction = &(FunctionManager::getInstance().registerGroup( childName, FunctionBase::CatalogInterface::factory( catalogKey, childName, this ) ) );

  return m_apertureFunction;
}


void ContactRelationBase::setSchemaDeviations( xmlWrapper::xmlNode,
                                               xmlWrapper::xmlNode schemaParent,
                                               integer )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "1";

    xmlWrapper::xmlNode tableFunctionNode = targetChoiceNode.prepend_child( "xsd:element" );
    tableFunctionNode.append_attribute( "name" ) = "TableFunction";
    tableFunctionNode.append_attribute( "type" ) = "TableFunctionType";
  }
}



void ContactRelationBase::initializePreSubGroups()
{
  TableFunction * const apertureTable = dynamic_cast< TableFunction * >(m_apertureFunction);
  if( apertureTable!=nullptr )
  {
    ArrayOfArraysView< real64 > xvals0 = apertureTable->getCoordinates();
    array1d< real64 > & yvals = apertureTable->getValues();

    GEOSX_ERROR_IF( xvals0.size() > 1,
                    "Aperture limiter table cannot be greater than a 1d table." );

    arraySlice1d< real64 > xvals = xvals0[0];
    localIndex const size = xvals.size();

    GEOSX_ERROR_IF( xvals0( 0, size-1 ) > 0.0 || xvals0( 0, size-1 ) < 0.0,
                    "Invalid aperture limiter table. Last coordinate must be zero!!" );

    GEOSX_ERROR_IF( xvals.size() < 2,
                    "Invalid aperture limiter table. Must have more than two points specified" );

    localIndex n=xvals.size()-1;
    real64 const slope = (yvals[n]-yvals[n-1]) / (xvals[n]-xvals[n-1]);

    GEOSX_ERROR_IF( slope >= 1.0,
                    "Invalid aperture table. slope of last two points >= 1 is invalid." );

    real64 m_apertureTransition = (yvals[n] - slope * xvals[n] ) / ( 1.0 - slope );

    xvals0.emplaceBack( 0, m_apertureTransition );
    yvals.emplace_back( m_apertureTransition );

    xvals0.emplaceBack( 0, m_apertureTransition*10e9 );
    yvals.emplace_back( m_apertureTransition*10e9 );

    apertureTable->reInitializeFunction();
  }


}

void ContactRelationBase::computeTraction( arraySlice1d< real64 const > const & dispJump,
                                           arraySlice1d< real64 > const & tractionVector ) const
{
  tractionVector[0] = dispJump[0] >= 0 ? 0.0 : m_penaltyStiffness * dispJump[0];
  tractionVector[1] = 0.0;
  tractionVector[2] = 0.0;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ContactRelationBase, string const &, Group * const )

}
} /* namespace geosx */
