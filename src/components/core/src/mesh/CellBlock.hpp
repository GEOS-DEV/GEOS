//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    static constexpr auto numFacesPerElementString     = "numFacesPerElement";
    static constexpr auto faceListString               = "faceList";

    dataRepository::ViewKey numNodesPerElement = { numNodesPerElementString };
    dataRepository::ViewKey nodeList           = { nodeListString };
    dataRepository::ViewKey numFacesPerElement = { numFacesPerElementString };
    dataRepository::ViewKey faceList           = { faceListString };
  } viewKeys;

  class groupKeyStruct
  {
public:
  } groupKeys;





  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }
  localIndex       & numNodesPerElement()       { return m_numNodesPerElement; }
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }
  localIndex       & numFacesPerElement()       { return m_numFacesPerElement; }

//protected:

  FixedOneToManyRelation  m_toNodesRelation;
  FixedOneToManyRelation  m_toFacesRelation;
  localIndex m_numNodesPerElement;
  localIndex m_numFacesPerElement;

  CellBlock& operator=(const CellBlock& rhs);
//  string & m_elementType;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
