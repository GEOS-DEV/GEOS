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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ElementTester.cpp
 * @author walsh24
 * @date Feb 18, 2011
 */

#include "SolverFactory.h"
#include "ElementTester.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"
#include "ArrayT/Array3dT.h"

// Finite element libraries
#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/GaussQuadrature.h"
#include "ElementLibrary/LagrangeBasis.h"

using namespace PS_STR;

#include <iostream>
#include <fstream>


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// Least-Squares Finite element advection diffusion reaction solver

ElementTester::ElementTester( const std::string& name,
                              ProblemManagerT* const pm ):
  SolverBase(name,pm),
  m_CommPtr(&(pm->m_epetraComm))
{}

ElementTester::~ElementTester()
{
  // TODO Auto-generated destructor stub
}

void ElementTester::ReadXML( TICPP::HierarchicalDataNode* const hdn ) {}


void ElementTester::RegisterFields( PhysicalDomainT * domain )
{

  domain->m_feNodeManager.AddKeyedDataField<FieldInfo::referencePosition>();
  domain->m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain->m_feNodeManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();

}

/**
 *
 *
 **/
double ElementTester::TimeStep( const realT& time,
                                const realT& dt,
                                const int cycleNumber,
                                PhysicalDomainT * domain,
                                const array<string>& namesOfSolverRegions,
                                SpatialPartition& partition,
                                FractunatorBase* const fractunator )
{
  realT dt_taken = dt;

  int myPID = m_CommPtr->MyPID();

  if(myPID == 0)
  {

    for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegionIter = domain->m_feElementManager.m_ElementRegions.begin() ;
         elementRegionIter != domain->m_feElementManager.m_ElementRegions.end() ;
         ++elementRegionIter )
    {
      ElementRegionT& elementRegion = elementRegionIter->second;
      // Store K and RHS of Kf=P
      TestRegion(time,dt,elementRegion,domain);
    }

  }

  return dt_taken;
}


/// returns K in row column value (r,c,v) format, and right hand side Kf = P
void ElementTester::TestRegion(const realT& time,
                               const realT& dt,
                               ElementRegionT& elementRegion,
                               PhysicalDomainT * domain){

  LagrangeBasis<3> feBasis(1);
  GaussQuadrature<3> feQuadrature(1);
  FiniteElement<3> feElement(feBasis, feQuadrature );
//  const unsigned int numQuadraturePoints = feElement.n_quadrature_points();

  const localIndex numElements = elementRegion.m_numElems;
  numNodesPerElem = elementRegion.m_numNodesPerElem;

  // Node properties
  nodeCoords.resize(numNodesPerElem);

  // loop over elements
  for( localIndex el = 0 ; el < numElements ; ++el)
  {

    // Rearrange nodes to match lexographic ordering
    nodeList = elementRegion.m_toNodesRelation[el];

    // Get coordinates
    for( int a=0 ; a<numNodesPerElem ; ++a )
    {
      const localIndex nd = nodeList[a];
      nodeCoords[a] = (*domain->m_feNodeManager.m_refposition)[nd];
      nodeCoords[a] += (*domain->m_feNodeManager.m_displacement)[nd];
    }

    // Update finite element geometry and record element data
    feElement.reinit(nodeCoords);

    std::cout << std::endl;
    std::cout << "Element " << el << std::endl;

    ReportCoordinates();

    CheckHandedness();

    CheckElement(feElement);

    std::cout << "**Rotating element 90 deg in x/y**" << std::endl;
    for( array<R1Tensor>::size_type i =0 ; i < nodeCoords.size() ; ++i)
    {
      realT x = nodeCoords[i][0];
      realT y = nodeCoords[i][1];
      nodeCoords[i][0] = -y;
      nodeCoords[i][1] =  x;
    }

    CheckElement(feElement);

    std::cout << "**Rotating element 90 deg in x/z**" << std::endl;
    for(array<R1Tensor>::size_type i =0 ; i < nodeCoords.size() ; ++i)
    {
      realT x = nodeCoords[i][0];
      realT z = nodeCoords[i][2];
      nodeCoords[i][0] =  z;
      nodeCoords[i][2] = -x;
    }

    CheckElement(feElement);

  } //loop over elements
}


void ElementTester::ReportCoordinates(){

  // Coordinates
  std::cout << "Node coordinates" << std::endl;
  for(int i =0 ; i < numNodesPerElem ; ++i)
  {
    std::cout << "    " << i << ": " << nodeCoords[i] << std::endl;
  }
}

void ElementTester::CheckHandedness(){
  R1Tensor dx,dy,dz;
  dx = nodeCoords[1] - nodeCoords[0];
  dy = nodeCoords[2] - nodeCoords[0];
  dz = nodeCoords[4] - nodeCoords[0];

  R1Tensor dxb,dyb,dzb;
  if(numNodesPerElem > 7)
  {
    dxb = nodeCoords[7] - nodeCoords[6];
    dyb = nodeCoords[7] - nodeCoords[5];
    dzb = nodeCoords[7] - nodeCoords[3];
  }

  std::cout << "Element handedness: " << ((dx*Cross(dy,dz) >0) ? 1 : (dx*Cross(dy,dz) <0) ? -1 : 0) << std::endl;
  if(numNodesPerElem > 7)
  {
    std::cout << "   Opposite Corner: " << ((dxb*Cross(dyb,dzb) >0) ? 1 : (dxb*Cross(dyb,dzb) <0) ? -1 : 0)  << std::endl;
  }
}

void ElementTester::CheckElement(FiniteElement<3>& feElement){
  feElement.reinit(nodeCoords);

  R1Tensor gradXYZ[3];
  gradXYZ[0] = 0;
  gradXYZ[1] = 0;
  gradXYZ[2] = 0;

  for(int i =0 ; i < 8 ; ++i)
  {
    gradXYZ[0] += feElement.gradient(i,0) * nodeCoords[i][0];
    gradXYZ[1] += feElement.gradient(i,0) * nodeCoords[i][1];
    gradXYZ[2] += feElement.gradient(i,0) * nodeCoords[i][2];
  }

  std::cout << "dNdX*X: Should equal Identity matrix" << std::endl;
  std::cout << "    " << gradXYZ[0] << std::endl;
  std::cout << "    " << gradXYZ[1] << std::endl;
  std::cout << "    " << gradXYZ[2] << std::endl;

  realT JQw = feElement.JxW(0);
  std::cout << "Jacobian x Quadrature weight (Should be positive): " << JQw << std::endl;
}

//////////////////////////////////////////
/// Register solver in the solver factory
REGISTER_SOLVER( ElementTester )
