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
 * @file InitialConditions.h
 * @author walsh24
 * @date Feb 28, 2011
 */

#ifndef InitialConditionFACTORY_H_
#define InitialConditionFACTORY_H_

#include "Common/Common.h"
#include "Common/typedefs.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "Utilities/StringUtilities.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ObjectManagers/PhysicalDomainT.h"

#include <map>
#include <string>
#include <vector>
#include "../../legacy/Common/GPException.h.old"


class ProblemManagerT;

/// Base class for all initial conditions
class InitialConditionBase
{
public:
  InitialConditionBase( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  InitialConditionBase( const PhysicalDomainT::ObjectDataStructureKeys& objectType, const std::string& regionName,
                        const std::string& fieldName,  const  std::vector<std::string>& setNames,  const FieldType& fieldType);
  virtual ~InitialConditionBase();

  virtual void RegisterFields( PhysicalDomainT& domain ) = 0;

  virtual void Apply( PhysicalDomainT& domain )=0;

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  bool IsThisInitialStress();

protected:
  PhysicalDomainT::ObjectDataStructureKeys objectType_;
  std::string regionName_;
  std::string fieldName_;
  std::vector<std::string> setNames_;
  FieldType fieldType_;
  bool m_additive;
  realT m_scale;
};


/// Set initial conditions for field based on data in an ascii file
class ReadInitialConditionFromFile : public InitialConditionBase
{
public:
  ReadInitialConditionFromFile( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~ReadInitialConditionFromFile(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  static const char* InitialConditionName(){return "ReadInitialConditionFromFile";};

protected:
  std::string filename_;
  bool isIndexedFile_;
};

/// Set initial condition for a field or subset thereof to a constant value
class ConstantInitialCondition : public InitialConditionBase
{
public:
  ConstantInitialCondition( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~ConstantInitialCondition(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn );

  static const char* InitialConditionName(){return "ConstantInitialCondition";};

protected:
  std::string valueStr_;
};

/// Set intial condition for a field or subset thereof based on a table
class InitialConditionTable : public InitialConditionBase
{
public:
  InitialConditionTable( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~InitialConditionTable(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn );

  static const char* InitialConditionName(){return "InitialConditionTable";};

protected:
  template <class TABLE, class ARRAY, class ARRAY2>
  inline static bool SetField(const TABLE* ptable, int component, const ARRAY* pos, const NodeManager& nodeManager, ElementRegionT& region,
                              const array1d<lSet>& sets, ARRAY2* pfield)
  {
    if((!pfield) || (!ptable))
      return false;
    const TABLE& table = *ptable;
    ARRAY2& field = *pfield;
    const bool isFE = pos == 0;
    if(isFE)
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
          {
            R1Tensor cpos(region.GetElementCenter(*si, nodeManager));
            field[*si][component] = table.Lookup(cpos);
          }
        }
      }
      else
      {
        for (localIndex i = 0 ; i < region.DataLengths() ; ++i)
        {
          R1Tensor cpos(region.GetElementCenter(i, nodeManager));
          field[i][component] = table.Lookup(cpos);
        }
      }
    }
    else
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
            field[*si][component] = table.Lookup((*pos)[*si]);
        }
      }
      else
      {
        localIndex i = 0;
        for (typename ARRAY::const_iterator it = pos->begin() ; it != pos->end() ; ++it, ++i)
          field[i][component] = table.Lookup(*it);
      }
    }
    return true;
  }

  template <class TABLE, class ARRAY, class ARRAY2>
  inline static bool SetField(const TABLE* ptable, const ARRAY* pos, const NodeManager& nodeManager, ElementRegionT& region, const array1d<lSet>& sets,
                              ARRAY2* pfield)
  {
    if((!pfield) || (!ptable))
      return false;
    const TABLE& table = *ptable;
    ARRAY2& field = *pfield;
    const bool isFE = pos == 0;
    if(isFE)
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
          {
            R1Tensor cpos(region.GetElementCenter(*si, nodeManager));
            field[*si] = table.Lookup(cpos);
          }
        }
      }
      else
      {
        for (localIndex i = 0 ; i < region.DataLengths() ; ++i)
        {
          R1Tensor cpos(region.GetElementCenter(i, nodeManager));
          field[i] = table.Lookup(cpos);
        }
      }
    }
    else
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
            field[*si] = table.Lookup((*pos)[*si]);
        }
      }
      else
      {
        localIndex i = 0;
        for (typename ARRAY::const_iterator it = pos->begin() ; it != pos->end() ; ++it, ++i)
          field[i] = table.Lookup(*it);
      }
    }
    return true;
  }

  template <class TABLE, class ARRAY, class ARRAY2>
  inline static bool AddToField(const TABLE* ptable, int component, const ARRAY* pos, const NodeManager& nodeManager, ElementRegionT& region,
                                const array1d<lSet>& sets, ARRAY2* pfield, realT& scale)
  {
    if((!pfield) || (!ptable))
      return false;
    const TABLE& table = *ptable;
    ARRAY2& field = *pfield;
    const bool isFE = pos == 0;
    if(isFE)
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
          {
            R1Tensor cpos(region.GetElementCenter(*si, nodeManager));
            field[*si][component] += table.Lookup(cpos) * scale;
          }
        }
      }
      else
      {
        for (localIndex i = 0 ; i < region.DataLengths() ; ++i)
        {
          R1Tensor cpos(region.GetElementCenter(i, nodeManager));
          field[i][component] += table.Lookup(cpos) * scale;
        }
      }
    }
    else
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
            field[*si][component] += table.Lookup((*pos)[*si]) * scale;
        }
      }
      else
      {
        localIndex i = 0;
        for (typename ARRAY::const_iterator it = pos->begin() ; it != pos->end() ; ++it, ++i)
          field[i][component] += table.Lookup(*it) * scale;
      }
    }
    return true;
  }

  template <class TABLE, class ARRAY, class ARRAY2>
  inline static bool AddToField(const TABLE* ptable, const ARRAY* pos, const NodeManager& nodeManager, ElementRegionT& region, const array1d<lSet>& sets,
                                ARRAY2* pfield, realT& scale)
  {
    if((!pfield) || (!ptable))
      return false;
    const TABLE& table = *ptable;
    ARRAY2& field = *pfield;
    const bool isFE = pos == 0;
    if(isFE)
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
          {
            R1Tensor cpos(region.GetElementCenter(*si, nodeManager));
            field[*si] += table.Lookup(cpos) * scale;
          }
        }
      }
      else
      {
        for (localIndex i = 0 ; i < region.DataLengths() ; ++i)
        {
          R1Tensor cpos(region.GetElementCenter(i, nodeManager));
          field[i] += table.Lookup(cpos) * scale;
        }
      }
    }
    else
    {
      if(sets.size() > 0)
      {
        for (array1d<lSet>::const_iterator iset = sets.begin() ; iset != sets.end() ; ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin() ; si != set.end() ; ++si)
            field[*si] += table.Lookup((*pos)[*si]) * scale;
        }
      }
      else
      {
        localIndex i = 0;
        for (typename ARRAY::const_iterator it = pos->begin() ; it != pos->end() ; ++it, ++i)
          field[i] += table.Lookup(*it) * scale;
      }
    }
    return true;
  }

  std::string m_tableName;
  string_array m_variableNames;
  array1d<FieldType> m_variableTypes;
  int m_component;
};

/// Set intial condition for a field or subset thereof based on a function of
// the other field variables
class InitialConditionFunction : public InitialConditionBase
{
public:
  InitialConditionFunction( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~InitialConditionFunction(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  static const char* InitialConditionName(){return "InitialConditionFunction";};

protected:
  std::string functionName_;
  string_array variableNames_;
  array1d<FieldType> variableTypes_;
  int component_;

};

///// Set intial condition for a field or subset thereof based on a 3D table
//class InitialConditionFrom3DTable: public InitialConditionBase
//{
//public:
//  InitialConditionFrom3DTable( const TICPP::HierarchicalDataNode* const
// InitialConditionNode, const ProblemManagerT* const problemManager);
//  virtual ~InitialConditionFrom3DTable(){};
//
//  virtual void RegisterFields( PhysicalDomainT& domain );
//
//  virtual void Apply( PhysicalDomainT& domain );
//
//  virtual void ReadXML( const TICPP::HierarchicalDataNode* const hdn );
//
//  static const char* InitialConditionName(){return
// "InitialConditionFrom3DTable";};
//
//protected:
//  std::string m_tableName;
//};
//


/// Calculate the center of each face from the node positions
class CalculateFaceCenters : public InitialConditionBase
{
public:
  CalculateFaceCenters( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateFaceCenters(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML(  TICPP::HierarchicalDataNode* hdn  ){};

  static const char* InitialConditionName(){return "CalculateFaceCenters";};

};

/// Calculate the center of each face from the node positions
class CalculateElementCenters : public InitialConditionBase
{
public:
  CalculateElementCenters( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateElementCenters(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML(  TICPP::HierarchicalDataNode* hdn  ){};

  static const char* InitialConditionName(){return "CalculateElementCenters";};

};

/// Calculate the normal of each face from the common plane contact
class CalculateFaceNormals : public InitialConditionBase
{
public:
  CalculateFaceNormals(  TICPP::HierarchicalDataNode*  InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateFaceNormals(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML(  TICPP::HierarchicalDataNode*  hdn  ){};

  static const char* InitialConditionName(){return "CalculateFaceNormals";};

};


/// Calculate an aperture of each face from the common plane contact
class CalculateAperture : public InitialConditionBase
{
public:
  CalculateAperture(  TICPP::HierarchicalDataNode*  InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateAperture(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn  ){};

  static const char* InitialConditionName(){return "CalculateAperture";};

};

/// Join sets
class LinkFractureFaces : public InitialConditionBase
{
public:
  LinkFractureFaces(  TICPP::HierarchicalDataNode*  InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~LinkFractureFaces(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn );

  static const char* InitialConditionName(){return "JoinFaces";};

protected:

};

//////////////////////////

// Initial Condition Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newInitialCondition"
//   * A base class to derive the functions to generate InitialCondition
// pointers: "InitialConditionInitializer"
//   * A String-to-InitialCondition-Intializer map hidden behind the
// getInitialConditionCatalogue function
//   * A template to create InitialCondition initializers:
// "InitialConditionRegistrator"
//   * A compiler directive to simplify autoregistration:
// "REGISTER_InitialCondition"
//
// Most initial conditions will only need to use one or two of the parts:
//   * To register a new InitialCondition in the factory:
// REGISTER_InitialCondition( InitialConditionClassName )
//   * To load a InitialCondition pointer from the factory:
//       InitialConditionBase* aInitialConditionPtr =
// newInitialCondition(InitialConditionString, args );

/// The InitialCondition Factory.
InitialConditionBase* newInitialCondition(const std::string& InitialConditionName, TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);

/// Base class to generate new InitialCondition pointers
class InitialConditionInitializer
{
public:
  virtual InitialConditionBase* initializeInitialCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm) = 0;
  virtual ~InitialConditionInitializer() = 0;
};
inline InitialConditionInitializer::~InitialConditionInitializer() { }

/// Interface to the InitialCondition name -> InitialCondition initializer map
std::map<std::string, InitialConditionInitializer*> & getInitialConditionCatalogue();

/// Return a list of supported InitialCondition names
void getInitialConditionNames( std::vector<std::string>& nameList);

/// Template for creating classes derived from InitialConditionInitializer
template< class InitialConditionType >
class InitialConditionRegistrator : public InitialConditionInitializer
{

public:
  InitialConditionRegistrator(void){
    std::string InitialConditionName = std::string(InitialConditionType::InitialConditionName() );
    getInitialConditionCatalogue() [InitialConditionName] = this;
  };

  InitialConditionBase* initializeInitialCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm  ){
    return new InitialConditionType( hdn, pm );
  };
};

/// Compiler directive to simplify autoregistration
#define REGISTER_InitialCondition( ClassName ) namespace { InitialConditionRegistrator<ClassName> reg_ ## ClassName; }

#endif
