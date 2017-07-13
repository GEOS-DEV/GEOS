/**
 * @file FaceManager.h
 * @author settgast1
 */

#ifndef FACEMANAGER_H_
#define FACEMANAGER_H_

#include "managers/ObjectManagerBase.hpp"

namespace geosx
{


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


//  void BuildFaces( const NodeManager& nodeManager, const ElementManager& elemManager );


public:



private:

  FaceManager() = delete;
  FaceManager( FaceManager const &) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
