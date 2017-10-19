/**
 * @file FaceManager.h
 * @author settgast1
 */

#ifndef FACEMANAGER_H_
#define FACEMANAGER_H_

#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class NodeManager;
class ElementRegionManager;
class CellBlockSubRegion;

class FaceManager: public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static string CatalogName()
  {
    return "FaceManager";
  }

  string getCatalogName() const override final
  {
    return FaceManager::CatalogName();
  }


  ///@}
  ///
  ///
  ///
  ///
  FaceManager( string const & , ManagedGroup * const parent );
  virtual ~FaceManager();

//  void Initialize(  ){}

  void FillDocumentationNode( dataRepository::ManagedGroup * const  ) override final;


  void BuildFaces( NodeManager const * const nodeManager, ElementRegionManager * const elemManager );

  void  AddNewFace( localIndex const & kReg,
                    localIndex const & kSubReg,
                    localIndex const & ke,
                    localIndex const & kelf,
                    localIndex & numFaces,
                    Array1dT<lArray1d>& facesByLowestNode,
                    lArray1d& tempNodeList,
                    Array1dT<lArray1d>& tempFaceToNodeMap,
                    CellBlockSubRegion const & elementRegion );



  void SortAllFaceNodes( NodeManager const & nodeManager,
                         ElementRegionManager const & elemManager);

  void SortFaceNodes( NodeManager const & nodeManager,
                      R1Tensor const & elementCenter,
                      const localIndex faceIndex );

  struct viewKeysStruct
  {
    dataRepository::ViewKey nodeList              = { "nodeList" };
    dataRepository::ViewKey elementRegionList     = { "elemRegionList" };
    dataRepository::ViewKey elementSubRegionList  = { "elemSubRegionList" };
    dataRepository::ViewKey elementList           = { "elemList" };
  }viewKeys;

  struct groupKeysStruct
  {
  }groupKeys;

  array< int32_array > const & nodeList() const        { return this->getReference< array< int32_array > >(viewKeys.nodeList); }
  array< int32_array > & nodeList()                    { return this->getReference< array< int32_array > >(viewKeys.nodeList); }
  Array2dT<int32> const & elementRegionList() const    { return this->getReference< Array2dT<int32> >(viewKeys.elementRegionList); }
  Array2dT<int32> & elementRegionList()                { return this->getReference< Array2dT<int32> >(viewKeys.elementRegionList); }
  Array2dT<int32> const & elementSubRegionList() const { return this->getReference< Array2dT<int32> >(viewKeys.elementSubRegionList); }
  Array2dT<int32> & elementSubRegionList()             { return this->getReference< Array2dT<int32> >(viewKeys.elementSubRegionList); }
  Array2dT<int32> const & elementList() const          { return this->getReference< Array2dT<int32> >(viewKeys.elementList); }
  Array2dT<int32> & elementList()                      { return this->getReference< Array2dT<int32> >(viewKeys.elementList); }




private:

  FaceManager() = delete;
  FaceManager( FaceManager const &) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
