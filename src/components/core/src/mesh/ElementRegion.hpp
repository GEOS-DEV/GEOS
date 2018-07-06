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
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTREGION_H
#define ELEMENTREGION_H

#include "managers/ObjectManagerBase.hpp"
#include "FaceManager.hpp"



class StableTimeStep;

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const defaultMaterial = "material";
//string const numNodesPerElement = "numNodesPerElement";
//string const nodeList = "nodeList";
//string const constitutiveMap = "constitutiveMap";
string const numericalMethod = "numericalMethod";
string const cellBlockSubRegions = "cellBlockSubRegions";
string const cellBlockSubRegionNames = "cellBlocks";
}
}

class CellBlockSubRegion;



/**
 * Class to manage the data stored at the element level.
 */
class ElementRegion : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "ElementRegion"; }

  virtual const string getCatalogName() const override final
  { return ElementRegion::CatalogName(); }


  ///@}


  ElementRegion() = delete;

  ElementRegion( string const & name, ManagedGroup * const parent );


  ElementRegion(const ElementRegion& init);
//  ElementRegion( ElementRegion&& init);


  virtual void FillDocumentationNode() override final;

  virtual void ReadXML_PostProcess() override;

  void SetConstitutiveMap( dataRepository::ManagedGroup const * domain,
                           map<string,localIndex> & counts );

  virtual ~ElementRegion() override;

//  Array2dT<integer> & m_toNodesRelation;


  virtual void InitializePreSubGroups( ManagedGroup * const group ) override final;

  subGroupMap & GetSubRegions()
  {
    return GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups();
  }

  subGroupMap const & GetSubRegions() const
  {
    return GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups();
  }

  CellBlockSubRegion const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(regionName);
  }
  CellBlockSubRegion * GetSubRegion( string const & regionName )
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(regionName);
  }

  CellBlockSubRegion const * GetSubRegion( localIndex const & index ) const
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(index);
  }
  CellBlockSubRegion * GetSubRegion( localIndex const & index )
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(index);
  }

  localIndex numSubRegions() const
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups().size();
  }
  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda )
  {
    ManagedGroup * cellBlockSubRegions = this->GetGroup(dataRepository::keys::cellBlockSubRegions);

    cellBlockSubRegions->forSubGroups<CellBlockSubRegion>( [&]( CellBlockSubRegion * subRegion ) -> void
      {
        lambda( subRegion );
      });

//    set<string> names;
//    for( auto & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
//    {
//      names.insert(iterCellBlocks.first);
//    }
//    for( auto & name : names )
//    {
//      lambda( cellBlockSubRegions->GetGroup<CellBlockSubRegion>(name) );
//    }

  }


  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda ) const
  {
    ManagedGroup const * cellBlockSubRegions = this->GetGroup(dataRepository::keys::cellBlockSubRegions);

    cellBlockSubRegions->forSubGroups<CellBlockSubRegion>( [&]( CellBlockSubRegion const * subRegion ) -> void
      {
        lambda( subRegion );
      });

//        set<string> names;
//    for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
//    {
//      names.insert(iterCellBlocks.first);
//    }
//
//    for( auto const & name : names )
//    {
//      lambda( cellBlockSubRegions->GetGroup<CellBlockSubRegion>(name) );
//    }

  }


  string const & getNumericalMethod() const
  {
    return this->getData<string>(dataRepository::keys::numericalMethod);
  }


private:
  ElementRegion& operator=(const ElementRegion& rhs);
//  string & m_elementType;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
