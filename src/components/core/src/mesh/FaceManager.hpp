/**
 * @file FaceManager.h
 * @author settgast1
 */

#ifndef FACEMANAGER_H_
#define FACEMANAGER_H_

#include "common/InterObjectRelation.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class NodeManager;
class ElementRegionManager;
class CellBlockSubRegion;

class FaceManager : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName() 
  { return "FaceManager"; }

  virtual const string getCatalogName() const override final
  { return FaceManager::CatalogName(); }


  ///@}
  ///
  ///
  ///
  ///
  FaceManager( string const &, ManagedGroup * const parent );
  virtual ~FaceManager() override final;

//  void Initialize(  ){}

  virtual void FillDocumentationNode() override final;


  void BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elemManager );

  void  AddNewFace( localIndex const & kReg,
                    localIndex const & kSubReg,
                    localIndex const & ke,
                    localIndex const & kelf,
                    localIndex & numFaces,
                    array<localIndex_array>& facesByLowestNode,
                    localIndex_array& tempNodeList,
                    array<localIndex_array>& tempFaceToNodeMap,
                    CellBlockSubRegion const & elementRegion );



  void SortAllFaceNodes( NodeManager const & nodeManager,
                         ElementRegionManager const & elemManager);

  void SortFaceNodes( NodeManager const & nodeManager,
                      R1Tensor const & elementCenter,
                      const localIndex faceIndex );

  void SetDomainBoundaryObjects( NodeManager * const nodeManager );

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual int PackUpDownMapsSize( localIndex_array const & packList ) const override;
  virtual int PackUpDownMaps( buffer_unit_type * & buffer,
                              localIndex_array const & packList ) const override;

  virtual int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                localIndex_array const & packList ) override;


  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & nodeManager,
                                                   array<globalIndex_array>& faceToNodes ) override final;
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nodeListString              = "nodeList";
    static constexpr auto edgeListString              = "edgeList";
    static constexpr auto elementRegionListString     = "elemRegionList";
    static constexpr auto elementSubRegionListString  = "elemSubRegionList";
    static constexpr auto elementListString           = "elemList";

    dataRepository::ViewKey nodeList              = { nodeListString };
    dataRepository::ViewKey edgeList              = { edgeListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };
  } viewKeys;

  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;

  OrderedVariableOneToManyRelation const & nodeList() const        { return m_nodeList; }
  OrderedVariableOneToManyRelation & nodeList()                    { return m_nodeList; }
  Array2dT<localIndex> const & elementRegionList() const    { return this->getReference< Array2dT<localIndex> >(viewKeys.elementRegionList); }
  Array2dT<localIndex> & elementRegionList()                { return this->getReference< Array2dT<localIndex> >(viewKeys.elementRegionList); }
  Array2dT<localIndex> const & elementSubRegionList() const { return this->getReference< Array2dT<localIndex> >(viewKeys.elementSubRegionList); }
  Array2dT<localIndex> & elementSubRegionList()             { return this->getReference< Array2dT<localIndex> >(viewKeys.elementSubRegionList); }
  Array2dT<localIndex> const & elementList() const          { return this->getReference< Array2dT<localIndex> >(viewKeys.elementList); }
  Array2dT<localIndex> & elementList()                      { return this->getReference< Array2dT<localIndex> >(viewKeys.elementList); }



private:

  template<bool DOPACK>
  int PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const;


  OrderedVariableOneToManyRelation m_nodeList;
  Array2dT<localIndex>  m_elementRegionList;
  Array2dT<localIndex> m_elementSubRegionList;
  Array2dT<localIndex> m_elementList;

  FaceManager() = delete;
  FaceManager( FaceManager const &) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
