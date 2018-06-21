// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file ContactManagerT.cpp
 * @author Randolph Settgast
 * @date created on Sep 17, 2010
 */

#include "ContactManagerBaseT.h"
#include "ContactManagerT.h"

/**
 * @brief Constructor for the contacts / neighbor pairs
 * @author Scott Johnson
 * @date Jun 15, 2011
 * We are currently assuming that "contact" is the union of all "neighbor pairs"
 * and "inactive + active physical contacts"
 * Note that neighbor pairs that are also contacts would be the intersection of
 * these two sets
 * This data structure could be split into different structures later if it is
 * found that these are too cumbersome
 */
ContactManagerT::ContactManagerT():
  ContactManagerBaseT(ObjectDataStructureBaseT::ContactBaseManager),
  m_contactToIntersectionPolygonPointsMap(m_VariableOneToManyMaps["contactToIntersectionPolygonPointsMap"]),
  m_intersectionPolygonPoints(),
  ctr(0.0),
  dd(1e-4)
{
  this->AddKeylessDataField<R1Tensor>( "face1ParentSoln", true, true);
  this->AddKeylessDataField<R1Tensor>( "face2ParentSoln", true, true);
}

ContactManagerT::~ContactManagerT()
{}

localIndex ContactManagerT::NumberOfPolygons() const
{
  array<lArray1d>& tmp = m_contactToIntersectionPolygonPointsMap;
  localIndex num = 0;
  for( array<lArray1d>::size_type i = 0 ; i < tmp.size() ; ++i)
    num += tmp[i].size() > 0 ? 1 : 0;
  return num;
}

/**
 * @brief FOR DEFAULT VISUALIZATION - visit needs at least one polygon to
 * visualize a mesh or it gets angry
 * @author Scott Johnson
 * Calculates the center as the center of the domain and assigns the polygon a
 * dimension of 1e-4 the minimum dimension
 * @param[in] xmin Domain minima
 * @param[in] xmax Domain maxima
 */
void ContactManagerT::SetDefaultPolygonDimensions(const R1Tensor& xmin, const R1Tensor& xmax)
{
  R1Tensor dx = xmax;
  if(xmax.MaxVal() > 1e10 && xmin.MinVal() > 1e10)
  {
    this->ctr = 0.0;
    dx = 1e10;
  }
  else
  {
    this->ctr = xmax;
    this->ctr += xmin;
    this->ctr *= 0.5;
    dx -= xmin;
  }
  this->dd = dx.MinVal();
  this->dd *= 1e-4;
}
