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
 * @file ContactRelationBase.cpp
 */

#include "ContactRelationBase.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "managers/Functions/TableFunction.hpp"



namespace geosx
{

using namespace dataRepository;


namespace constitutive
{


ContactRelationBase::ContactRelationBase( string const & name,
                                          ManagedGroup * const parent ):
  ConstitutiveBase(name, parent),
  m_penaltyStiffness(0.0),
  m_apertureFunction(nullptr)
{
  RegisterViewWrapper( viewKeyStruct::penaltyStiffnessString, &m_penaltyStiffness, 0 )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Value of the penetration penalty stiffness. Units of Pressure/length");


}

ContactRelationBase::~ContactRelationBase()
{
}




ManagedGroup *
ContactRelationBase::CreateChild( string const & catalogKey, string const & childName )
{
  ManagedGroup * rval = nullptr;

  FunctionBase::CatalogInterface::CatalogType const & functionCatalog = FunctionBase::GetCatalog();
  if( functionCatalog.count(catalogKey) )
  {
    m_apertureFunction = (FunctionBase::CatalogInterface::Factory( catalogKey, childName, this )).release();
    rval =  NewFunctionManager::Instance()->RegisterGroup( childName, m_apertureFunction, 0 );
  }
  else
  {
    GEOS_ERROR(catalogKey<<" is an invalid key ContactRelationBase child group.");
  }
  return rval;
}

void ContactRelationBase::InitializePreSubGroups( ManagedGroup * const )
{
  TableFunction * const apertureTable = dynamic_cast<TableFunction*>(m_apertureFunction);
  if( apertureTable!=nullptr )
  {
    array1d<array1d<real64>> & xvals0 = apertureTable->getCoordinates();
    array1d<real64> &          yvals = apertureTable->getValues();

    GEOS_ERROR_IF( xvals0.size() > 1,
                   "Aperture limiter table cannot be greater than a 1d table.");

    array1d<real64> & xvals = xvals0[0];

    GEOS_ERROR_IF( xvals.back() > 0.0 || xvals.back() < 0.0 ,
                   "Invalid aperture limiter table. Last coordinate must be zero!!" );

    GEOS_ERROR_IF( xvals.size() < 2,
                   "Invalid aperture limiter table. Must have more than two points specified");

    localIndex n=xvals.size()-1;
    real64 const slope = (yvals[n]-yvals[n-1]) / (xvals[n]-xvals[n-1]) ;

    GEOS_ERROR_IF( slope >= 1.0,
                   "Invalid aperture table. slope of last two points >= 1 is invalid.");

    real64 m_apertureTransition = (yvals[n] - slope * xvals[n] ) / ( 1.0 - slope );

    xvals.push_back(m_apertureTransition);
    yvals.push_back(m_apertureTransition);

    xvals.push_back(m_apertureTransition*10e9);
    yvals.push_back(m_apertureTransition*10e9);

    apertureTable->reInitializeFunction();
  }

//  for( int i=0 ; i<200 ; ++i )
//  {
//    real64 coord = 0.01*i-1.0;
//    std::cout<<coord<<" "<<apertureTable->Evaluate( &coord )<<std::endl;
//  }

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ContactRelationBase, string const &, ManagedGroup * const )

}
} /* namespace geosx */

