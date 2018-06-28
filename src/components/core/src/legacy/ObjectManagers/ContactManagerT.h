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
 * @file ContactManagerT.h
 * @author Scott Johnson
 * @date Jun 15, 2011
 */
#ifndef CONTACTMANAGERT_H_
#define CONTACTMANAGERT_H_

#include "../../dataRepository/Group.hpp"
#include "ContactManagerBaseT.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

/**
 * @author Scott Johnson
 * @brief Class to manage collections of face-face contacts
 */
class ContactManagerT : public ContactManagerBaseT
{
public:
  /**
   * @brief Contact manager constructor
   * @author Scott Johnson
   */
  ContactManagerT();//ExternalFaceManagerT* fm);

  /**
   * @brief Contact manager destructor
   * @author Scott Johnson
   */
  virtual ~ContactManagerT();

  localIndex NumberOfPolygons() const;

public:
  void SetDefaultPolygonDimensions(const R1Tensor& xmin, const R1Tensor& xmax);

  inline realT DefaultPolygonDimension() const { return dd;}

  inline R1Tensor DefaultPolygonCenter() const { return ctr;}

  ///maps the discrete element to the indices in the _external_ face manager
  OrderedVariableOneToManyRelation& m_contactToIntersectionPolygonPointsMap;

  ///stateless cache of the points comprising the contact patches
  array<R1Tensor> m_intersectionPolygonPoints;

private:
  // FOR DEFAULT VISUALIZATION
  R1Tensor ctr;
  realT dd;
};

#endif /* CONTACTMANAGERT_H_ */
