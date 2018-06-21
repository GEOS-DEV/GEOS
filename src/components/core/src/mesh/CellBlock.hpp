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
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTOBJECTT_H_
#define ELEMENTOBJECTT_H_

#include "managers/ObjectManagerBase.hpp"
#include "legacy/ObjectManagers/EnergyT.h"
//#include "common/InterObjectRelation.hpp"
//#include "legacy/ArrayT/bufvector.h"
#include "FaceManager.hpp"


class StableTimeStep;

namespace geosx
{

namespace dataRepository
{
namespace keys
{
//string const defaultMaterial = "material";
//string const numNodesPerElement = "numNodesPerElement";
//string const nodeList = "nodeList";
//string const constitutiveMap = "constitutiveMap";
}
}



/**
 * Class to manage the data stored at the element level.
 */
class CellBlock : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "CellBlock"; }

  virtual const string getCatalogName() const override final
  { return CellBlock::CatalogName(); }


  ///@}


  CellBlock() = delete;

  CellBlock( string const & name, ManagedGroup * const parent );


  CellBlock(const CellBlock& init);
//  ElementRegion( ElementRegion&& init);


  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

//  map<string,integer> SetConstitutiveMap( dataRepository::ManagedGroup const &
// domain );

  virtual ~CellBlock() override;

  void GetFaceNodes( const localIndex elementIndex,
                     const localIndex localFaceIndex,
                     localIndex_array& nodeIndicies) const;

  R1Tensor GetElementCenter(localIndex k, const NodeManager& nodeManager, const bool useReferencePos = true) const;



//  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;
//
//  virtual int PackUpDownMapsSize( localIndex_array const & packList ) const override;
//
//  virtual int PackUpDownMaps( buffer_unit_type * & buffer,
//                              localIndex_array const & packList ) const override;
//
//  virtual int UnpackUpDownMaps( buffer_unit_type const * & buffer,
//                                localIndex_array const & packList ) override;

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto numNodesPerElementString     = "numNodesPerElement";
    static constexpr auto nodeListString               = "nodeList";
    static constexpr auto numEdgesPerElementString     = "numEdgesPerElement";
    static constexpr auto edgeListString               = "edgeList";
    static constexpr auto numFacesPerElementString     = "numFacesPerElement";
    static constexpr auto faceListString               = "faceList";

    dataRepository::ViewKey numNodesPerElement = { numNodesPerElementString };
    dataRepository::ViewKey nodeList           = { nodeListString };
    dataRepository::ViewKey numEdgesPerElement = { numEdgesPerElementString };
    dataRepository::ViewKey edgeList           = { edgeListString };
    dataRepository::ViewKey numFacesPerElement = { numFacesPerElementString };
    dataRepository::ViewKey faceList           = { faceListString };
  } m_CellBlockViewKeys;

//  class groupKeyStruct
//  {
//  } groupKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockViewKeys; }

//  virtual groupKeyStruct & groupKeys() { return m_ObjectManagerBaseGroupKeys; }
//  virtual groupKeyStruct const & groupKeys() const { return m_ObjectManagerBaseGroupKeys; }





  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }
  localIndex       & numNodesPerElement()       { return m_numNodesPerElement; }
  localIndex const & numEdgesPerElement() const { return m_numEdgesPerElement; }
  localIndex       & numEdgesPerElement()       { return m_numEdgesPerElement; }
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }
  localIndex       & numFacesPerElement()       { return m_numFacesPerElement; }

  FixedOneToManyRelation & nodeList()                    { return m_toNodesRelation; }
  FixedOneToManyRelation const & nodeList() const        { return m_toNodesRelation; }

  FixedOneToManyRelation       & edgeList()       { return m_toEdgesRelation; }
  FixedOneToManyRelation const & edgeList() const { return m_toEdgesRelation; }

  FixedOneToManyRelation       & faceList()       { return m_toFacesRelation; }
  FixedOneToManyRelation const & faceList() const { return m_toFacesRelation; }


//protected:

private:
  localIndex m_numNodesPerElement;
  localIndex m_numEdgesPerElement;
  localIndex m_numFacesPerElement;
  FixedOneToManyRelation  m_toNodesRelation;
  FixedOneToManyRelation  m_toEdgesRelation;
  FixedOneToManyRelation  m_toFacesRelation;


  CellBlock& operator=(const CellBlock& rhs);
//  string & m_elementType;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
