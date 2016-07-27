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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * ExternalFaceManagerT.cpp
 * @date Jun 15, 2011
 * @author Scott Johnson
 */

//-------------------------
//NOTE: there was a major trimming of this class after r590; revert to r590 to see the previous methods
//-------------------------

#include "FaceManagerT.h"
#include "ExternalFaceManagerT.h"
#include "PhysicalDomainT.h"
#include "ContactManagerT.h"
#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif
#include "Utilities/GeometryUtilities.h"
#include "IO/BinStream.h"
#include "ArrayT/Array3dT.h"
#include <iostream> // to be removed - dumps shape
#include <limits>
#include <cassert>

#ifdef SRC_INTERNAL
#include "PhysicsSolvers/BackgroundAMR.h"
#include "MPI_Communications/Bifroest/GeodynBifroest.h"
#endif

#include "Constitutive/Interface/InterfaceFactory.h"

/**
 * @brief Constructor to set internal pointers to the external node and face managers
 * @author Scott Johnson
 * @param fm External face manager pointer
 */
ExternalFaceManagerT::ExternalFaceManagerT(FaceManagerT* fm = 0, FaceManagerT* fmDiscreteElement = 0) :
ObjectDataStructureBaseT(ObjectDataStructureBaseT::ExternalFaceManager),
m_faceManager(NULL),
m_discreteElementFaceManager(NULL),
m_neighborList(m_VariableOneToManyMaps["neighborList"]),
m_neighborListInverse(m_UnorderedVariableOneToManyMaps["neighborListInverse"]),
m_contactActive(false),
m_contactSelf(false),
m_autoContact(true),
m_contact(NULL),
m_externalFaceToFaceMap(m_OneToOneMaps["externalFaceToFaceMap"]),
nfe(0),
nde(0),
m_sorted(false),
m_sorter(NULL)
{
  this->m_tol.maximumSeparation = 0.0;
  this->m_tol.searchRadiusFactor = 0.0;
  this->m_tol.searchRadiusVelocityFactor = 0.0;
  this->m_tol.cosMin = 0.7;
  this->m_tol.spatial = 1.0e-14;
  this->m_tol.area = 1.0e-6;
  this->m_tol.penetration = 0.9;
  this->m_tol.feParentSolution = 1.0e-12;

  this->m_faceManager = fm;
  this->m_discreteElementFaceManager = fmDiscreteElement;

  m_contact = InterfaceFactory::NewInterface("PenaltyCoulomb");
  m_contact->SetVariableParameters(false);
  m_contact->resize(0,1);

  //states for contact library interface
  this->AddKeylessDataField<realT> ( "normalApproach", true, true);
  this->AddKeylessDataField<realT> ( "normalApproachMax", true, true);
  this->AddKeylessDataField<realT> ( "area", true, true);
  this->AddKeylessDataField<R1Tensor>( "shearSlip", true, true);

  this->AddKeylessDataField<realT> ( "boundingRadiusLastSort" , true, false );
  this->AddKeylessDataField<R1Tensor> ( "faceCenterLastSort" , true, false );
  this->AddKeylessDataField<realT> ( "boundingRadius" , true, true );
  this->AddKeylessDataField<R1Tensor> ( "faceCenter" , true, true );
  this->AddKeylessDataField<R1Tensor> ( "faceNormal" , true, true );
  this->AddKeylessDataField<R1Tensor> ( "faceVelocity" , true, true );

  this->AddKeylessDataField<localIndex>( "parentElement" , true, false );
  this->AddKeylessDataField<int>( "parentElementRegion" , true, false );

  this->AddKeylessDataField<int>( "excludeFromContact" , true, true );

  //@annavarapusr1: Adding stabilization and weighting fields for Nitsche's method
  this->AddKeylessDataField<realT> ( "nitscheStab_n" , false, false );
  this->AddKeylessDataField<realT> ( "nitscheStab_t1" , false, false );
  this->AddKeylessDataField<realT> ( "nitscheStab_t2" , false, false );
  this->AddKeylessDataField<realT> ( "nitscheGamma" , false, false );

  m_sorter = SpatialSorting::SpatialSorterFactory::NewSpatialSorter("N2");
}

ExternalFaceManagerT::~ExternalFaceManagerT()
{
  if(m_contact)
    delete m_contact;

  if(m_sorter)
    delete m_sorter;
}

void ExternalFaceManagerT::SplitFace(const localIndex parentIndex,
                                     const localIndex newFaceIndex,
                                     const NodeManagerT& nodeManager)
{
  //THIS IS ONLY FOR INTERNAL FE FACES BEING SPLIT IN TWAIN
  const localIndex nfe1 = nfe + 1;

  //(1) Update base
  insert(nfe);
  insert(nfe1);

  //(2) Update external face to face map
  m_externalFaceToFaceMap[nfe] = parentIndex;
  m_externalFaceToFaceMap[nfe1] = newFaceIndex;

  //(3) Update the face to external face map and counts
  lArray1d& externalFaceIndex = m_faceManager->GetFieldData<localIndex>("externalFaceIndex");
  externalFaceIndex[parentIndex] = nfe;
  externalFaceIndex[newFaceIndex] = nfe1;

  //(4) Update DE indices to reflect change in face index offset (due to change in face manager size)
  for(localIndex i = nfe+2; i < m_DataLengths; ++i)
    ++m_externalFaceToFaceMap[i]; //i.e., one _new_ face has been added to the face manager, so the face index offset is increased by 1
  // FIXME: Fu note: This is a bug.  In 2D, two new faces are added.  Need to take care of this.

  //(5) Update neighbor list and inverse
  SpatialSorting::SpatialSorterBase::Add(nfe, nfe1, m_neighborList, m_neighborListInverse);

  //(6) Update counts
  nfe += 2;

  //(7) Resolve ghosting issues
  PostSortUpdate(nodeManager); //TODO: check this logic!!!!
}


void ExternalFaceManagerT::SplitFace(const localIndex parentIndex,
                                     const localIndex newFaceIndex1,
                                     const localIndex newFaceIndex2,
                                     const NodeManagerT& nodeManager)
{
  //THIS IS ONLY FOR INTERNAL FE FACES BEING SPLIT IN TWAIN
  const localIndex nfe1 = nfe + 1;
  const localIndex nfe2 = nfe + 2;

  //(1) Update base
  insert(nfe);
  insert(nfe1);
  insert(nfe2);

  //(2) Update external face to face map
  m_externalFaceToFaceMap[nfe] = parentIndex;
  m_externalFaceToFaceMap[nfe1] = newFaceIndex1;
  m_externalFaceToFaceMap[nfe2] = newFaceIndex2;

  //(3) Update the face to external face map and counts
  lArray1d& externalFaceIndex = m_faceManager->GetFieldData<localIndex>("externalFaceIndex");
  externalFaceIndex[parentIndex] = nfe;
  externalFaceIndex[newFaceIndex1] = nfe1;
  externalFaceIndex[newFaceIndex2] = nfe2;

  //(4) Update DE indices to reflect change in face index offset (due to change in face manager size)
  for(localIndex i = nfe+3; i < m_DataLengths; ++i)
  {
    m_externalFaceToFaceMap[i] += 2; //i.e., one _new_ face has been added to the face manager, so the face index offset is increased by 1
  // FIXME: Fu note: This is a bug.  In 2D, two new faces are added.  Need to take care of this.
  }
  //(5) Update neighbor list and inverse
  SpatialSorting::SpatialSorterBase::Add(nfe, nfe1, m_neighborList, m_neighborListInverse);
  SpatialSorting::SpatialSorterBase::Add(nfe, nfe2, m_neighborList, m_neighborListInverse);

  //(6) Update counts
  nfe += 3;

  //(7) Resolve ghosting issues
  PostSortUpdate(nodeManager); //TODO: check this logic!!!!
}

globalIndex ExternalFaceManagerT::insert(const localIndex i, const bool assignGlobals )
{

  globalIndex gi = ObjectDataStructureBaseT::insert( i, assignGlobals );
  m_contact->insert(i);
  return gi;
}
void ExternalFaceManagerT::erase( const localIndex i )
{
  ObjectDataStructureBaseT::erase(i);
  m_contact->erase(i);
}
globalIndex ExternalFaceManagerT::resize( const localIndex size, const bool assignGlobals )
{
  globalIndex gi = ObjectDataStructureBaseT::resize( size, assignGlobals );
  m_contact->resize(size);
  return gi;
}

void ExternalFaceManagerT::DeserializeObjectField(const std::string& name, const rArray1d& field)
{
  if(m_DataLengths == 0)
    return;
  m_contact->SetValues(name, field);
}

void ExternalFaceManagerT::DeserializeObjectFields(const sArray1d& names, const Array1dT<rArray1d>& fields)
{
  if(m_DataLengths == 0)
    return;
  m_contact->SetValues(names, fields);
}

void ExternalFaceManagerT::ReadXML(TICPP::HierarchicalDataNode* ContactNode)
{
  m_contactActive = ContactNode->GetAttributeOrDefault<bool> ("active", true);
  if(!m_contactActive)
    return;

  m_smoothedContact = ContactNode->GetAttributeOrDefault<bool> ("smoothed", false);
  m_contactSelf = ContactNode->GetAttributeOrDefault<bool> ("allowSelfContact", false);
  m_autoContact = ContactNode->GetAttributeOrDefault<bool> ("auto", true);

  m_tol.maximumSeparation =
      ContactNode->GetAttributeOrDefault<realT> ("maximumSeparation", 0);
  m_tol.searchRadiusFactor =
      ContactNode->GetAttributeOrDefault<realT> ("searchRadiusFactor", 0);
  m_tol.searchRadiusVelocityFactor =
      ContactNode->GetAttributeOrDefault<realT> ("searchRadiusVelocityFactor", 0);
  m_tol.feParentSolution =
      ContactNode->GetAttributeOrDefault<realT> ("feParentSolnTol", 1e-12);

  //Contact parameters
  this->m_contact->ReadXML(*ContactNode->Next(true));

  //Tolerances
  //    apertureFactor(0.0)  = % of face radius by which to increase the radial search distance
  //                           (i.e., the greater this value the larger the atomic search cost but the less frequent)
  //    cosMin(0.2)          = minimum cos of the angle between two opposing faces below which contact is rejected
  //                           (i.e., if faces are too far away from parallel, contact is no longer assessed; 0 is orthogonal)
  //    spatial(1e-12)       = ratio of distance to minimum dimension of the face used as a threshold for several values
  //                           (i.e., very small distances are neglected)
  //    area(1e-6)           = ratio of common plane area to face area, below which contact is rejected
  //                           (i.e., very small contact areas are rejected)
  //    penetration(0.7)     = ratio of penetration distance to face minimum dimension, above which contact is rejected
  //                           (i.e., overpenetration causes loss of contact)
  //                           the penetration distance is relative to the common plane, so if you have a cube
  //                           you can assure that you avoid opposing sides "seeing" each other with a value < 0.5/sqrt(2)
  //    maximumAperture(0.0) = additional search distance in the "negative"
  //                           (i.e., out-of-contact normal) direction to use to assess contact
  //
  //    NOTE:
  //    radius_of_search = (1 + apertureFactor) * radius + maximumAperture
  m_tol.cosMin      = ContactNode->GetAttributeOrDefault<realT> ("cosMinTol", 0.2);
  m_tol.spatial     = ContactNode->GetAttributeOrDefault<realT> ("spatialTol", 1e-12);
  m_tol.area        = ContactNode->GetAttributeOrDefault<realT> ("areaTol", 1e-6);
  m_tol.penetration = ContactNode->GetAttributeOrDefault<realT> ("penetrationTol", 0.7);
}

void ExternalFaceManagerT::BuildExternalFaces()
{
  this->nfe = 0;

  const iArray1d& isExternal = m_faceManager->m_isExternal;

  // get number of faces from FEM
  //set the external face index for those faces that are external
  {
    lArray1d& externalFaceIndex = m_faceManager->GetFieldData<localIndex>("externalFaceIndex");
    localIndex a = 0;
    for( iArray1d::const_iterator i=isExternal.begin() ; i!=isExternal.end() ; ++i, ++a)
    {
      externalFaceIndex[a] = (*i == 1) ? (this->nfe++) : std::numeric_limits<localIndex>::max();
    }
  }

  //all DE faces are external
  this->nde = this->m_discreteElementFaceManager->DataLengths();

  const localIndex numExternalFaces = this->nfe + this->nde;

  // assign the external face manager size to be that of all external faces (i.e., DE and FE)
  this->resize( numExternalFaces );

  // assign FE face values to first nfe entries of m_externalFaceToFaceMap
  // assign values of m_externalFaceToFaceMap
  lArray1d::iterator iter_xfc=this->m_externalFaceToFaceMap.begin();
  {
    iArray1d::const_iterator iter_fcfe=isExternal.begin() ;
    for(localIndex ifcfe = 0 ; iter_fcfe!=isExternal.end() ; ++iter_fcfe, ++ifcfe )
    {
      if( *iter_fcfe == 1 )
      {
        *iter_xfc = ifcfe;
        ++iter_xfc;
      }
    }
  }

  // externalFaceToFaceMap now contains indices of the FE faces in the first part
  // and DE face indices offset by the TOTAL number of FE faces (m_faceManager->DataLengths()) in the second part
  // see FaceIndex function for the back conversion call
  localIndex ixfc = this->DiscreteElementFaceIndexOffset();
  for( ; iter_xfc!=this->m_externalFaceToFaceMap.end() ; ++iter_xfc, ++ixfc)
  {
    *iter_xfc = ixfc;
  }

  // need to make sure the element properties are set before attempting this
  // in the interim, just use a user-defined bulk modulus set in the Contact tag
  //SetBulkModuli();
}

void ExternalFaceManagerT::BuildExternalFacesDiscreteElements()
{
  //all DE faces are external
  this->nde = this->m_discreteElementFaceManager->DataLengths();

  const localIndex numExternalFaces = this->nfe + this->nde;

  // assign the external face manager size to be that of all external faces (i.e., DE and FE)
  this->resize( numExternalFaces );

  localIndex ixfc = this->DiscreteElementFaceIndexOffset();
  for(localIndex i = this->nfe; i < numExternalFaces; ++i, ++ixfc)
  {
    m_externalFaceToFaceMap[i] = ixfc;
  }
}

/**
 * @brief Form a list of faces associated with each face in the list
 * @author Scott Johnson
 * Populates the m_neighborList structure and flags whether a repopulation was performed
 * Note: lastSort fields and current bounding sphere fields (and face velocity) are all filled by the end of the call
 * @param[in] nodeManager Node collection
 * @param[in] discreteElementNodeManager Node collection for the discrete elements
 * @param[in] m_discreteElementManager Collection of discrete elements
 * @param[in] dt Timestep
 */
bool ExternalFaceManagerT::RecalculateNeighborList(
    const NodeManagerT& nodeManager,
    const NodeManagerT& discreteElementNodeManager,
    const DiscreteElementManagerT& m_discreteElementManager,
    const realT dt, const bool forceRecalculate, const bool useReferencePosition)
{
  bool sort = false;
  if(!this->m_faceManager && !this->m_discreteElementFaceManager)
    return sort;
  if(!(m_autoContact || forceRecalculate))
    return sort;

  const localIndex num = this->DataLengths();
  if(num == 0)
    return sort;

  localIndex faceIndex = 0;
  bool is_fe = false;

  rArray1d& radii = this->GetFieldData<realT>("boundingRadius");
  Array1dT<R1Tensor>& centers = this->GetFieldData<R1Tensor>("faceCenter");
  Array1dT<R1Tensor>& normals = this->GetFieldData<R1Tensor>("faceNormal");
  Array1dT<R1Tensor>& vel = this->GetFieldData<R1Tensor>("faceVelocity");

  const iArray1d& excludeFromContact = this->GetFieldData<int>("excludeFromContact");

  //-----------------------------------------------------------------------
  // DETERMINE WHETHER ANYTHING NEEDS TO BE RESORTED
  //-----------------------------------------------------------------------
  lSet toResort;
  if(this->m_sorted)
  {
    SetExcludeFromContact( nodeManager, false );

    const rArray1d& lradii = this->GetFieldData<realT>("boundingRadiusLastSort");
    const Array1dT<R1Tensor>& lcenters = this->GetFieldData<R1Tensor>("faceCenterLastSort");
    R1Tensor dx;

    //iterate through external faces
    for (localIndex ixfc1 = 0; ixfc1 < this->DataLengths(); ++ixfc1)
    {
      if(excludeFromContact[ixfc1]>0)
        continue;

      //get the index in the appropriate faceManager instance
      faceIndex = this->FaceIndex(ixfc1, is_fe);
      if (is_fe)
      {
        this->m_faceManager->FaceBoundingSphere(nodeManager,
                                                faceIndex, dt,
                                                centers[ixfc1],
                                                vel[ixfc1],
                                                radii[ixfc1],
                                                normals[ixfc1],
                                                useReferencePosition);
      }
      else
      {
        this->m_discreteElementFaceManager->FaceBoundingSphere(discreteElementNodeManager,
                                                               faceIndex, dt,
                                                               centers[ixfc1],
                                                               vel[ixfc1],
                                                               radii[ixfc1],
                                                               normals[ixfc1],
                                                               useReferencePosition);
      }

      //determine whether my bounding sphere is at or crossing the last radius boundary
      {
        const realT lradius = lradii[ixfc1] - 2.0*radii[ixfc1];
        if(!SpatialSorting::SpatialSorterBase::Close(radii[ixfc1],
                                                     centers[ixfc1],
                                                     lradius,
                                                     lcenters[ixfc1]))
        {
          toResort.insert(ixfc1);
        }
      }
    }
    sort = toResort.size() > 0;
  }
  else
  {
    SetExcludeFromContact( nodeManager, true );

    sort = true;
    this->m_sorted = true;

    //iterate through
    const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = this->m_faceManager->m_toElementsRelation;
    const Array1dT<lArray1d>& fton = this->m_discreteElementFaceManager->m_toNodesRelation;
    const lArray1d& ntod = discreteElementNodeManager.GetFieldData<FieldInfo::demIndex>();//s.m_nodeToDiscreteElementMap;

    iArray1d faceAttachedToALocalNode;
    SetGhostRank(nodeManager, discreteElementNodeManager, excludeFromContact, faceAttachedToALocalNode);

    lArray1d& parentElement = this->GetFieldData<localIndex>("parentElement");
    iArray1d& parentElementRegion = this->GetFieldData<int>("parentElementRegion");

    for (localIndex ixfc1 = 0; ixfc1 < this->DataLengths(); ++ixfc1)
    {
      if(excludeFromContact[ixfc1]>0)
        continue;
      faceIndex = this->FaceIndex(ixfc1, is_fe);
      if (is_fe)
      {
        if(faceAttachedToALocalNode[ixfc1]==0)
        {
          parentElement[ixfc1] = std::numeric_limits<localIndex>::max();
          parentElementRegion[ixfc1] = std::numeric_limits<int>::max();
        }
        else
        {
          parentElement[ixfc1] = ftoe[faceIndex][0].second;
          parentElementRegion[ixfc1] = ftoe[faceIndex][0].first->m_regionNumber;
        }

        this->m_faceManager->FaceBoundingSphere(nodeManager, faceIndex, dt,
                                                centers[ixfc1], vel[ixfc1],
                                                radii[ixfc1], normals[ixfc1],
                                                useReferencePosition);
      }
      else
      {
        parentElement[ixfc1] = faceAttachedToALocalNode[ixfc1]==0 ? std::numeric_limits<localIndex>::max() :
            ntod[fton[faceIndex][0]];

        parentElementRegion[ixfc1] = faceAttachedToALocalNode[ixfc1]==0 ? std::numeric_limits<int>::max() : 0;

        this->m_discreteElementFaceManager->FaceBoundingSphere(discreteElementNodeManager,
                                                               faceIndex, dt,
                                                               centers[ixfc1], vel[ixfc1],
                                                               radii[ixfc1], normals[ixfc1],
                                                               useReferencePosition);
      }
    }
  }

  //-----------------------------------------------------------------------
  // LEAVE IF NO SORTING IS REQUIRED
  //-----------------------------------------------------------------------
  if(!sort)
    return false;


  //update the last radii and last centers
  {
    rArray1d& lradii = this->GetFieldData<realT>("boundingRadiusLastSort");
    Array1dT<R1Tensor>& lcenters = this->GetFieldData<R1Tensor>("faceCenterLastSort");

    const realT rvelFct = 0.5 * dt * this->m_tol.searchRadiusVelocityFactor;

    //reset the values of the last sort <-- remove to do partial updates!!!!
    toResort.clear();
    if (toResort.size() == 0)
    {
      //-----------------------------------------------------------------------
      // NOTHING TO RESORT BUT SORT TRIGGERED? THEN DO TOTAL RESORT
      //-----------------------------------------------------------------------
      std::copy(centers.begin(), centers.end(), lcenters.begin());
      std::copy(radii.begin(), radii.end(), lradii.begin());

      //now, increment the search radius by the adjustments
      for (localIndex ixfc1 = 0; ixfc1 < this->DataLengths(); ++ixfc1)
      {
        R1Tensor rvelTmp(vel[ixfc1]);
        rvelTmp *= rvelFct;
        IncrementSearchRadius(rvelTmp.L2_Norm(), lradii[ixfc1]);
        lcenters[ixfc1] += rvelTmp;
      }
      m_sorter->Sort(lradii, lcenters,
                     this->m_neighborList,
                     this->m_neighborListInverse,
                     excludeFromContact.data());
    }
    else
    {
      //-----------------------------------------------------------------------
      // PARTIAL RESORT
      //-----------------------------------------------------------------------
      for (std::set<localIndex>::const_iterator it = toResort.begin();
          it != toResort.end(); ++it)
      {
        R1Tensor rvelTmp(vel[*it]);
        rvelTmp *= rvelFct;
        lcenters[*it] = centers[*it];
        lcenters[*it] += rvelTmp;
        lradii[*it] = radii[*it];
        IncrementSearchRadius(rvelTmp.L2_Norm(), lradii[*it]);
      }
      m_sorter->Update(lradii, lcenters, toResort,
                       this->m_neighborList,
                       this->m_neighborListInverse);
    }
  }

  //-----------------------------------------------------------------------
  // DEAL WITH GHOSTS AND SELF-CONTACT
  //-----------------------------------------------------------------------
  PostSortUpdate(nodeManager, discreteElementNodeManager, toResort);
  return true;
}

void ExternalFaceManagerT::PostSortUpdate( const NodeManagerT& nodeManager,
                                           const NodeManagerT& discreteElementNodeManager,
                                           const lSet& toResort)
{
  //-----------------------------------------------------------------------
  // DEAL WITH GHOSTS AND SELF-CONTACT
  //-----------------------------------------------------------------------
  const iArray1d& excludeFromContact = this->GetFieldData<int>("excludeFromContact");
  iArray1d faceAttachedToALocalNode;
  SetGhostRank(nodeManager, discreteElementNodeManager, excludeFromContact,
               faceAttachedToALocalNode);
  RemoveInvalidPairsFromNeighborList(faceAttachedToALocalNode, toResort);
}

void ExternalFaceManagerT::PostSortUpdate( const NodeManagerT& nodeManager)
{
  //-----------------------------------------------------------------------
  // DEAL WITH GHOSTS AND SELF-CONTACT
  //-----------------------------------------------------------------------
  lSet toResort;
  const iArray1d& excludeFromContact = this->GetFieldData<int>("excludeFromContact");
  iArray1d faceAttachedToALocalNode;
  SetGhostRank(nodeManager, excludeFromContact,
               faceAttachedToALocalNode);
  RemoveInvalidPairsFromNeighborList(faceAttachedToALocalNode, toResort);
}

void ExternalFaceManagerT::SetGhostRank(    const NodeManagerT& nodeManager,
                                            const NodeManagerT& discreteElementNodeManager,
                                            const iArray1d& excludeFromContact,
                                            iArray1d& faceAttachedToALocalNode)
{
  faceAttachedToALocalNode.resize(this->DataLengths());

  //create a temporary ghost flag array and an array to cache whether a face is attached to a local node
  iArray1d& ghostRank = this->GetFieldData<FieldInfo::ghostRank>();
  {
    const Array1dT<lArray1d>& faceToNodesFE = this->m_faceManager->m_toNodesRelation;
    const Array1dT<lArray1d>& faceToNodesDE = this->m_discreteElementFaceManager->m_toNodesRelation;
    const iArray1d& nodeGhostFE = nodeManager.GetFieldData<FieldInfo::ghostRank>();
    const iArray1d& nodeGhostDE = discreteElementNodeManager.GetFieldData<FieldInfo::ghostRank>();
    const iArray1d& ghostRankFE = this->m_faceManager->GetFieldData<FieldInfo::ghostRank>();
    const iArray1d& ghostRankDE = this->m_discreteElementFaceManager->GetFieldData<FieldInfo::ghostRank>();


    for (localIndex ixfc1 = 0; ixfc1 < this->DataLengths(); ++ixfc1)
    {
      faceAttachedToALocalNode[ixfc1] = 0;
      if(excludeFromContact[ixfc1]>0)
        continue;

      bool is_fe;
      const localIndex faceIndex = this->FaceIndex(ixfc1, is_fe);
      const lArray1d& faceToNodes = is_fe ? faceToNodesFE[faceIndex] : faceToNodesDE[faceIndex];
      const iArray1d& nodeGhost = is_fe ? nodeGhostFE : nodeGhostDE;
      const iArray1d& ghostRankT = is_fe ? ghostRankFE : ghostRankDE;
      for(lArray1d::const_iterator it = faceToNodes.begin(); it != faceToNodes.end(); ++it)
      {
        if(nodeGhost[*it] < 0)
        {
          faceAttachedToALocalNode[ixfc1] = 1;
          break;
        }
      }
      ghostRank[ixfc1] = ghostRankT[faceIndex];
    }
  }
}


void ExternalFaceManagerT::SetGhostRank(    const NodeManagerT& nodeManager,
                                            const iArray1d& excludeFromContact,
                                            iArray1d& faceAttachedToALocalNode)
{
  faceAttachedToALocalNode.resize(this->DataLengths());

  //create a temporary ghost flag array and an array to cache whether a face is attached to a local node
  iArray1d& ghostRank = this->GetFieldData<FieldInfo::ghostRank>();
  {
    const Array1dT<lArray1d>& faceToNodesFE = this->m_faceManager->m_toNodesRelation;
    const iArray1d& nodeGhostFE = nodeManager.GetFieldData<FieldInfo::ghostRank>();
    const iArray1d& ghostRankFE = this->m_faceManager->GetFieldData<FieldInfo::ghostRank>();
    for (localIndex ixfc1 = 0; ixfc1 < nfe; ++ixfc1)
    {
      faceAttachedToALocalNode[ixfc1] = 0;
      if(excludeFromContact[ixfc1]>0)
        continue;
      bool is_fe = true;
      const localIndex faceIndex = this->FaceIndex(ixfc1, is_fe);
      const lArray1d& faceToNodes = faceToNodesFE[faceIndex];
      for(lArray1d::const_iterator it = faceToNodes.begin(); it != faceToNodes.end(); ++it)
      {
        if(nodeGhostFE[*it] < 0)
        {
          faceAttachedToALocalNode[ixfc1] = 1;
          break;
        }
      }
      ghostRank[ixfc1] = ghostRankFE[faceIndex];
    }
  }
}


void ExternalFaceManagerT::RemoveInvalidPairsFromNeighborListSub(const lArray1d& parentElement,
                                                                 const iArray1d& parentElementRegion,
                                                                 const iArray1d& faceAttachedToALocalNode,
                                                                 const localIndex a,
                                                                 lArray1d& current)
{
  lArray1d::iterator it1 = current.begin();
  const bool is_fe0 = IsFiniteElement(a);
  while(it1 != current.end())
  {
    const localIndex b = *it1;
    const bool is_fe1 = IsFiniteElement(b);

    //remove if (1) neither face is attached to a local node and, hence, both are ghosts
    //          (2) they are part of the same element, are not mixed contacts (DEM-FEM), and do not allow self contact
    //              (either by flag for FE or by nature of the element for DE)
    const bool remove = (faceAttachedToALocalNode[a] == 0 && faceAttachedToALocalNode[b] == 0) ||
        ((parentElement[a] == parentElement[b]) &&
            (parentElementRegion[a] == parentElementRegion[b]) &&
            (is_fe0 == is_fe1) &&
            (!this->m_contactSelf || is_fe0==false));
    if(remove)
    {
      lSet::iterator itset = this->m_neighborListInverse[b].find(a);
      if(itset != this->m_neighborListInverse[b].end())
      {
        this->m_neighborListInverse[b].erase(itset);
      }
      it1 = current.erase(it1);
    }
    else
    {
      ++it1;
    }
  }
}

/**
 * @brief Remove any ghost-on-ghost pairs or invalid self-contact from the neighbor list
 * @author Scott Johnson
 */
void ExternalFaceManagerT::RemoveInvalidPairsFromNeighborList(const iArray1d& faceAttachedToALocalNode,
                                                              const lSet& toResort)
{
  const lArray1d& parentElement = this->GetFieldData<localIndex>("parentElement");
  const iArray1d& parentElementRegion = this->GetFieldData<int>("parentElementRegion");
  localIndex a = 0;
  if(toResort.size() == 0)
  {
    //go through the entire neighbor list ... this is a total sort
    for(Array1dT<lArray1d>::iterator it0 = this->m_neighborList.begin(); it0 != this->m_neighborList.end(); ++it0, ++a)
    {
      lArray1d& current = *it0;
      RemoveInvalidPairsFromNeighborListSub(parentElement, parentElementRegion, faceAttachedToALocalNode, a, current);
    }
  }
  else
  {
    //for everything in the resort list, you need to also check anything
    //that is proximate, so both the neighbor list and the inverse list
    lSet check(toResort);
    for(lSet::const_iterator it0 = toResort.begin(); it0 != toResort.end(); ++it0)
    {
      const lArray1d& current = this->m_neighborList[*it0];
      for(lArray1d::const_iterator it = current.begin(); it != current.end(); ++it)
        check.insert(*it);
      const lSet& currentSet = this->m_neighborListInverse[*it0];
      for(lSet::const_iterator it = currentSet.begin(); it != currentSet.end(); ++it)
        check.insert(*it);
    }

    for(lSet::const_iterator it0 = check.begin(); it0 != check.end(); ++it0)
    {
      lArray1d& current = this->m_neighborList[*it0];
      RemoveInvalidPairsFromNeighborListSub(parentElement, parentElementRegion, faceAttachedToALocalNode, *it0, current);
    }
  }
}

/**
 * @brief Update the geometric properties of each contact
 * @author Scott Johnson
 * This updates the contact geometry without updating any material states or stresses
 * In addition to the variables below, this also fills the common plane points as well as
 * the default polygon's dimensions
 * @param[in] dt Timestep duration
 * @param[in,out] domain Domain
 * @param[out] xs Nodal positions on the faces
 * @param[out] areaCommonPlanesC The area of each common plane
 * @param[out] areaFacesC The area of each face that composes a common plane
 */
void ExternalFaceManagerT::UpdateGeometricContactPropertiesSub(const realT dt,
                                                               PhysicalDomainT& domain,
                                                               const Array1dT<Array1dT<R1Tensor> >& xs  )
{
  //check the tolerance for cosMin ... important!
  if(this->m_tol.cosMin <= 0)
    throw GPException("Cannot handle orthogonal or non-convex surface contact currently ... set tol.cosMin > 0");

  //-----------------------------
  //TO SET -->
  //-----------------------------
  iArray1d& activeC = domain.m_contactManager.GetFieldData<int>( "active");// activeC = 0.0;
  rArray1d& areaC = domain.m_contactManager.GetFieldData<realT>( "area"); areaC = 0.0;
  rArray1d& area = this->GetFieldData<realT>("area"); area = 0.0; //area of a face that also belongs to a common plane
  Array1dT<R1Tensor>& xpolyPts = domain.m_contactManager.m_intersectionPolygonPoints;   xpolyPts.clear();
  Array1dT<lArray1d>& contactToXPoly = domain.m_contactManager.m_contactToIntersectionPolygonPointsMap;
  Array1dT<R1Tensor>& normalC   = domain.m_contactManager.GetFieldData<R1Tensor>( "normal"); normalC = 0.0;
  Array1dT<R1Tensor>& applicationPointC = domain.m_contactManager.GetFieldData<R1Tensor>( "applicationPoint");
  //  std::unique_ptr<InterfaceBase>& contactC   = domain.m_contactManager.m_contact;

  //-----------------------------
  //ALREADY SET -->
  //-----------------------------
  const Array1dT<R1Tensor>& xfc     = this->GetFieldData<R1Tensor> ( "faceCenter");
  //  const Array1dT<R1Tensor>& dxfc    = this->GetFieldData<R1Tensor> ( "faceVelocity");
  const Array1dT<R1Tensor>& nfc     = this->GetFieldData<R1Tensor> ( "faceNormal");

  //-----------------------------
  //-----------------------------
  //-----------------------------

  //set extrema
  R1Tensor tmp, xmin, xmax;
  {
    xmin = std::numeric_limits<realT>::max();
    xmax = xmin;
    xmax *= -1.0;
  }

  //face properties
  bool is_fe1 = false, is_fe2 = false;

  //Common plane properties
  R1Tensor centerCommonPlane, normalCommonPlane;
  Array1dT<R1Tensor> pointsCommonPlane;

  if(m_tol.cosMin <= 0)
    throw GPException("Cannot handle non-convex contact with the common plane");

  //  const realT big = std::numeric_limits<realT>::max(); //big Suitably large number used for initializing mins and maxes
  //  const realT small = 1e-6; //small Small number for tolerance

  //-------ITERATE THROUGH POSSIBLE CONTACTS---------

  localIndex index = 0;

  for(Array1dT<lArray1d>::size_type ixfc1 = 0; ixfc1 < this->m_neighborList.size(); ++ixfc1) {

    ExternalFaceStruct efs1;

    //Retrieve additional non-cached face properties for face #1
    localIndex faceIndex1 = this->FaceIndex(ixfc1, is_fe1);
    if (is_fe1)
      this->m_faceManager->FaceProperties(domain.m_feNodeManager, faceIndex1, dt, efs1);
    else
      this->m_discreteElementFaceManager->FaceProperties(domain.m_discreteElementSurfaceNodes,
                                                         faceIndex1, dt, efs1);
    xmax.SetMax(efs1.xmax);
    xmin.SetMin(efs1.xmin);

    //Iterate through secondary faces
    for(lArray1d::size_type it = 0; it < this->m_neighborList[ixfc1].size(); ++it, ++index) {

      const localIndex ixfc2 = this->m_neighborList[ixfc1][it];

      ExternalFaceStruct efs2;

      //Retrieve additional non-cached face properties for face #2
      localIndex faceIndex2 = this->FaceIndex(ixfc2, is_fe2);
      if (is_fe2)
        this->m_faceManager->FaceProperties(domain.m_feNodeManager, faceIndex2, dt, efs2);
      else
        this->m_discreteElementFaceManager->FaceProperties(domain.m_discreteElementSurfaceNodes,
                                                           faceIndex2, dt, efs2);
      xmax.SetMax(efs2.xmax);
      xmin.SetMin(efs2.xmin);

      //***CONTACT CONDITION***
      //Determine the contact condition for the two external faces
      if(index >= domain.m_contactManager.DataLengths())
        throw GPException("index incremented past the end of the contact manager fields!");

      contactToXPoly[index].clear();
      {
        if( domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension==2 )
        {
          R1Tensor p1, p2;
          const int newActive = CommonEdgeInterferenceGeometry(xfc[ixfc1], nfc[ixfc1], efs1,
                                                               xfc[ixfc2], nfc[ixfc2], efs2,
                                                               this->m_tol,
                                                               centerCommonPlane,
                                                               normalC[index],
                                                               areaC[index], pointsCommonPlane);
          applicationPointC[index] = centerCommonPlane;
          activeC[index] = UpdateContactFlag(activeC[index], newActive);
        }
        else
        {
          const int newActive = CommonPlaneInterferenceGeometry(xfc[ixfc1], nfc[ixfc1], efs1,
                                                                xfc[ixfc2], nfc[ixfc2], efs2,
                                                                this->m_tol,
                                                                centerCommonPlane,
                                                                applicationPointC[index],
                                                                normalC[index],
                                                                areaC[index], pointsCommonPlane);
          activeC[index] = UpdateContactFlag(activeC[index], newActive);
        }
      }

      //  If contact has been positively identified, then go on to update
      //  the neighbor list properties
      if(IsAnyContact(activeC[index]))
      {
        if(IsMechanicalContact(activeC[index]))
        {
          area[ixfc1] += areaC[index];
          area[ixfc2] += areaC[index];
        }

        //set polygon points for future visualization
        contactToXPoly[index].reserve(pointsCommonPlane.size());
        for(Array1dT<R1Tensor>::size_type ii = 0; ii < pointsCommonPlane.size(); ++ii)
        {
          contactToXPoly[index].push_back(xpolyPts.size());
          xpolyPts.push_back(pointsCommonPlane[ii]);
        }
      }
    }//ixfc2
  }//ixfc1

  //FOR DEFAULT VISUALIZATION: SET POLYGON EXTENTS; POLYGON APPENDED IN THE WRITESILO FUNCTION (IF NECESSARY)
  domain.m_contactManager.SetDefaultPolygonDimensions(xmin, xmax);
}


/**
 * @brief Get the face position and velocity
 * @author Scott Johnson
 * @param[in] normal average normal over the face
 * @param[in] applicationPoint point to query
 * @param[in] xs list of nodal positions
 * @param[in] vs list of nodal velocities
 * @param[in] is_fe flag indicating whether this is a finite element face with 4-nodes
 * @param[out] faceParentSolution parent solution to cache in the case of finite elements
 * @param[out] point position of the face
 * @param[out] velocity velocity of the face
 */
void ExternalFaceManagerT::FacePositionAndVelocityFE(const R1Tensor& normal,
                                                     const R1Tensor& applicationPoint,
                                                     const Array1dT<R1Tensor>& xs,
                                                     const Array1dT<R1Tensor>& vs,
                                                     const realT tolParentSolution,
                                                     R1Tensor& faceParentSolution,
                                                     R1Tensor& point,
                                                     R1Tensor& velocity)
{
  point = 0.0;
  velocity = 0.0;

  realT N[4];
  GeometryUtilities::FindProjectionInParentSpace(applicationPoint, normal, xs, faceParentSolution,
                                                 N, true, tolParentSolution);
  for (localIndex a = 0; a < 4; ++a)
  {
    R1Tensor temp(xs[a]);
    temp *= N[a];
    point += temp;

    temp = vs[a];
    temp *= N[a];
    velocity += temp;
  }
}

void ExternalFaceManagerT::FacePositionAndVelocityFE(const R1Tensor& normal,
                                                     const R1Tensor& applicationPoint,
                                                     const Array1dT<R1Tensor>& xs,
                                                     const Array1dT<R1Tensor>& vs,
                                                     R1Tensor& point,
                                                     R1Tensor& velocity)
{
  point = 0.0;
  velocity = 0.0;

  //get parent solution
  realT N[] = {0,0};
  {
    R1Tensor temp(xs[1]);
    temp -= xs[0];
    const realT L2 = Dot(temp, temp);

    temp = applicationPoint;
    temp -= xs[0];
    const realT L2a = Dot(temp, temp);

    N[0] = sqrt(L2a / L2);
    N[1] = 1.0 - N[0];
  }

  for (localIndex a = 0; a < 2; ++a)
  {
    R1Tensor temp(xs[a]);
    temp *= N[a];
    point += temp;

    temp = vs[a];
    temp *= N[a];
    velocity += temp;
  }
}

void ExternalFaceManagerT::FacePositionAndVelocityFE(const Array1dT<R1Tensor>& xs,
                                                     const Array1dT<R1Tensor>& vs,
                                                     R1Tensor& point,
                                                     R1Tensor& velocity)
{
  point = 0.0;
  velocity = 0.0;
  Array1dT<R1Tensor>::const_iterator itx = xs.begin();
  Array1dT<R1Tensor>::const_iterator itv = vs.begin();
  const realT fct = xs.size() > 0 ? 1.0/xs.size() : 0.0;
  for (;itx != xs.end(); ++itx, ++itv)
  {
    point += *itx;
    velocity += *itv;
  }
  point *= fct;
  velocity *= fct;
}

realT ExternalFaceManagerT::ApplyStress(const R1Tensor& normal,
                                        const R1Tensor& applicationPoint,
                                        const Array1dT<R1Tensor>& xs,
                                        const lArray1d& faceToNodes,
                                        const Array1dT<lSet>&  nodeToFaces,
                                        const rArray1d& masses,
                                        const bool is_fe,
                                        const R1Tensor& stress,
                                        const realT area,
                                        const realT tolParentSolution,
                                        R1Tensor& faceParentSolution,
                                        Array1dT<R1Tensor>& forces,
                                        Array1dT<R1Tensor>& contactForces)
{
  realT massFace = 0.0;
  if (is_fe)
  {
    realT N[4];
    GeometryUtilities::FindProjectionInParentSpace(applicationPoint, normal, xs, faceParentSolution,
                                                   N, true, tolParentSolution);
    R1Tensor temp1;
    for (lArray1d::size_type a = 0; a < faceToNodes.size(); ++a)
    {
      const localIndex inode = faceToNodes[a];
      temp1 = stress;
      temp1 *= area;
      temp1 *= -N[a];
      forces[inode] += temp1;
      contactForces[inode] += temp1;
      massFace += masses[inode] / nodeToFaces[inode].size();
    }
  }
  else
  {
    R1Tensor temp1;
    for (lArray1d::size_type a = 0; a < faceToNodes.size(); ++a)
    {
      const localIndex inode = faceToNodes[a];
      temp1 = stress;
      temp1 *= area;
      temp1 *= -1.0 / faceToNodes.size();
      forces[inode] += temp1;
      contactForces[inode] += temp1;
      massFace += masses[inode] / nodeToFaces[inode].size();
    }
  }
  return massFace;
}

/**
 * @brief Update the forces on the nodes (FE) or DE's due to the contact forces
 * @author Scott Johnson
 * @param dt Timestep duration
 * @param domain Domain object to update
 */
void ExternalFaceManagerT::UpdateAndApplyContactStresses(StableTimeStep& maxdt, const realT dt,
                                                         PhysicalDomainT& domain,
                                                         const Array1dT<Array1dT<R1Tensor> >& xs)
{
  //-----------------------------
  //TO SET -->
  //-----------------------------
  Array1dT<R1Tensor>& forcesFE = domain.m_feNodeManager.GetFieldData<FieldInfo::force>();
  Array1dT<R1Tensor>& contactForcesFE =
      domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();

  Array1dT<R1Tensor>& forcesDE =
      domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::force>();
  Array1dT<R1Tensor>& contactForcesDE = domain.m_discreteElementSurfaceNodes.GetFieldData<
      FieldInfo::contactForce>();

  const rArray1d& massDE = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::mass>();
  const rArray1d& massFE = domain.m_feNodeManager.GetFieldData<FieldInfo::mass>();

  //-----------------------------
  //SET BY UPDATEGEOMETRICPROPERTIES -->
  //-----------------------------
  const Array1dT<R1Tensor>& normalC = domain.m_contactManager.GetFieldData<R1Tensor>("normal");
  const Array1dT<R1Tensor>& velocityC = domain.m_contactManager.GetFieldData<R1Tensor>("velocity");
  const iArray1d& activeC = domain.m_contactManager.GetFieldData<int>("active");
  const rArray1d& areaC = domain.m_contactManager.GetFieldData<realT>("area");
  const rArray1d& normalApproachC = domain.m_contactManager.GetFieldData<realT>("normalApproach");
  //const rArray1d& normalApproachMaxC = domain.m_contactManager.GetFieldData<realT>("normalApproachMax");
  const Array1dT<R1Tensor>& applicationPointC = domain.m_contactManager.GetFieldData<R1Tensor>(
      "applicationPoint");
  Array1dT<R1Tensor>& shearSlipC = domain.m_contactManager.GetFieldData<R1Tensor>("shearSlip");
  Array1dT<R1Tensor>& face1ParentSoln = domain.m_contactManager.GetFieldData<R1Tensor>(
      "face1ParentSoln"); //set before, but overwritten by ApplyStress
  Array1dT<R1Tensor>& face2ParentSoln = domain.m_contactManager.GetFieldData<R1Tensor>(
      "face2ParentSoln"); //set before, but overwritten by ApplyStress

  InterfaceBase*& contactC = domain.m_contactManager.m_contact;

  Array1dT<R1Tensor>& shearSlip = this->GetFieldData<R1Tensor>(
      "shearSlip");
  const rArray1d& area = this->GetFieldData<realT>("area");

  const bool is2D = domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension==2;

  //-----------------------------
  //-----------------------------
  //-----------------------------
  //-----------------------------

  //-----ITERATE NEIGHBORLIST AND UPDATE FORCES AND STRESSES-----
  const OrderedVariableOneToManyRelation& deFaceToNodes = domain.m_discreteElementSurfaceFaces.m_toNodesRelation;
  localIndex index = 0;
  for (Array1dT<lArray1d>::size_type ixfc1 = 0; ixfc1 < this->m_neighborList.size(); ++ixfc1)
  {
    bool is_fe1;
    localIndex faceIndex1 = this->FaceIndex(ixfc1, is_fe1);
    for (lArray1d::size_type it = 0; it < this->m_neighborList[ixfc1].size(); ++it, ++index)
    {
      bool is_fe2;
      const localIndex ixfc2 = this->m_neighborList[ixfc1][it];
      localIndex faceIndex2 = this->FaceIndex(ixfc2, is_fe2);
      if (index >= domain.m_contactManager.DataLengths())
        throw GPException("index incremented past the end of the contact manager fields!");

      //FIRST, UPDATE THE CONTACT STATES THAT WILL BE USED FOR ANY UPDATE CASE
      InterfaceBaseStateData& contactStateC = *contactC->StateData(index, 0);
      const int iactiveC = activeC[index];

      //-----------------------------------------------------------
      //CONTACT 1: ONLY FOR NON-MECHANICAL CONTACT
      //           zero stresses and states
      if (!IsMechanicalContact(iactiveC)) //this is NOT a new or old mechanical contact
      {
        //zero out the stress
        contactStateC *= 0.0;
        shearSlipC[index] = 0.0;

        //update contact state cache
        contactStateC.dxndt = Dot(velocityC[index], normalC[index]);
        contactStateC.dt = dt;
        contactStateC.normalApproach = normalApproachC[index];

        //--contact 1b--
        //try to predict the timestep associated with an eminent contact
        const realT dvn = Dot(velocityC[index], normalC[index]);
        if(IsOpenContact(iactiveC) && dvn > 0)
        {
          realT massFace1 = 0.0;
          {
            const lArray1d& faceToNodes = is_fe1 ? domain.m_feFaceManager.m_toNodesRelation(faceIndex1) :
                deFaceToNodes(faceIndex1);
            const rArray1d& mass = is_fe1 ? massFE : massDE;
            const Array1dT<lSet>&  nodeToFaces = is_fe1 ? domain.m_feNodeManager.m_nodeToFaceMap :
                domain.m_discreteElementSurfaceNodes.m_nodeToFaceMap;
            for(localIndex a = 0; a < faceToNodes.size(); ++a)
              massFace1 += mass[faceToNodes[a]] / nodeToFaces[faceToNodes[a]].size();
          }
          realT massFace2 = 0.0;
          {
            const lArray1d& faceToNodes = is_fe2 ? domain.m_feFaceManager.m_toNodesRelation(faceIndex2) :
                deFaceToNodes(faceIndex2);
            const rArray1d& mass = is_fe2 ? massFE : massDE;
            const Array1dT<lSet>&  nodeToFaces = is_fe2 ? domain.m_feNodeManager.m_nodeToFaceMap :
                domain.m_discreteElementSurfaceNodes.m_nodeToFaceMap;
            for(localIndex a = 0; a < faceToNodes.size(); ++a)
              massFace2 += mass[faceToNodes[a]] / nodeToFaces[faceToNodes[a]].size();
          }

          //----GET THE STABLE TIMESTEP ESTIMATE ----
          maxdt.SetIfSmaller(StableTimestep(massFace1, massFace2, areaC[index],
                                            contactC->StiffnessProjected(
                                                index)));
          //-----------------------------------------
        }

        continue;
      }

      //-----------------------------------------------------------
      //CONTACT 2: ONLY FOR MECHANICAL CONTACT
      //           increment stresses
      //           note: convention is state1-state0 (stress returned is that on face 1)
      //           e.g., Dot(ncp, normal0) > 0 etc

      //IF THIS IS NEW, POPULATE FROM FACE STATES
      if (IsNewMechanicalContact(iactiveC))
      {
        //TODO: need to make sure "index" applies to the same one as previously, since a neighbor list update can invalidate the mapping; if fixed do NOT set activeC = 0.0 above

        //this is a new contact, so get initial state from constituent faces!
        //also set the accumulated plasticEnergy from the faces
        const realT fct1 = area[ixfc1];
        const realT fct2 = area[ixfc2];

        contactC->MapFromRegion(fct1, fct2,
                                *m_contact->StateData(ixfc1,0),
                                *m_contact->StateData(ixfc2,0),
                                index, 0);

        GeometryUtilities::MapFromRegion(shearSlip[ixfc1], shearSlip[ixfc2], fct1, fct2, shearSlipC[index], true);
      }

      R1Tensor dxs;
      contactStateC.UpdateOrientation(normalApproachC[index], dt, normalC[index], velocityC[index], dxs, shearSlipC[index]);

      //-----------------------------------------------------------
      //CONTACT 4: ONLY FOR FE-DE OR FE-FE CONTACT
      //           resolve the stresses
      contactC->StrainDrivenUpdate(index);
      contactC->SetPermeabilityTerm(index);

      //store the components of stress in the total stress tensor
      R1Tensor stress(normalC[index]);
      stress *= contactStateC.stress;
      stress += contactStateC.stressShearVector;

      const realT mass1 = ApplyStress(normalC[index],
                                      applicationPointC[index],
                                      xs[ixfc1],
                                      is_fe1 ? domain.m_feFaceManager.m_toNodesRelation(faceIndex1) :
                                          deFaceToNodes(faceIndex1),
                                      is_fe1 ? domain.m_feNodeManager.m_nodeToFaceMap :
                                          domain.m_discreteElementSurfaceNodes.m_nodeToFaceMap,
                                       is_fe1 ? massFE : massDE,
                                       is_fe1 && !is2D,
                                       stress,
                                       areaC[index],
                                       this->m_tol.feParentSolution,
                                       face1ParentSoln[index],
                                       is_fe1 ? forcesFE : forcesDE,
                                       is_fe1 ? contactForcesFE : contactForcesDE);

      R1Tensor tstress = stress;
      tstress *= -1.0;

      const realT mass2 = ApplyStress(normalC[index],
                                      applicationPointC[index],
                                      xs[ixfc2],
                                      is_fe2 ? domain.m_feFaceManager.m_toNodesRelation(faceIndex2) :
                                          deFaceToNodes(faceIndex2),
                                      is_fe2 ? domain.m_feNodeManager.m_nodeToFaceMap :
                                          domain.m_discreteElementSurfaceNodes.m_nodeToFaceMap,
                                      is_fe2 ? massFE : massDE,
                                      is_fe2 && !is2D,
                                      tstress,
                                      areaC[index],
                                      this->m_tol.feParentSolution,
                                      face2ParentSoln[index],
                                      is_fe2 ? forcesFE : forcesDE,
                                      is_fe2 ? contactForcesFE : contactForcesDE);

      //----GET THE STABLE TIMESTEP ESTIMATE ----
      maxdt.SetIfSmaller(StableTimestep(mass1, mass2, areaC[index],
                                        contactC->StiffnessProjected(index)));
      //-----------------------------------------

      //UPDATE THE SHEAR SLIP
      shearSlipC[index] += dxs;
    } //end foreach ixfc2
  } //end foreach ixfc1

  //-----------------------------------------------------------
  //CONTACT 5: redistribute the common plane states to the faces for evolution to the next step
  m_contact->ZeroStates();
  shearSlip = 0.0;

  index = 0;
  for (Array1dT<lArray1d>::size_type ixfc1 = 0; ixfc1 < this->m_neighborList.size(); ++ixfc1)
  {
    for (lArray1d::size_type it = 0; it < this->m_neighborList[ixfc1].size(); ++it, ++index)
    {
      const localIndex ixfc2 = this->m_neighborList[ixfc1][it];
      if (index >= domain.m_contactManager.DataLengths())
        throw GPException("index incremented past the end of the contact manager fields!");
      if (IsMechanicalContact(activeC[index]))
      {
        const realT fctNormal = areaC[index];
        const realT fct1 = area[ixfc1];
        const realT fct2 = area[ixfc2];

        contactC->MapToRegion(fctNormal, fct1, fct2, index, 0,
                              *m_contact->StateData(ixfc1,0),
                              *m_contact->StateData(ixfc2,0));

        GeometryUtilities::MapToRegion(fctNormal, fct1, fct2,
                                       shearSlipC[index],
                                       shearSlip[ixfc1],
                                       shearSlip[ixfc2], true);
      }
    }
  }
}

/**
 * @brief Update the forces on the nodes (FE) or DE's due to the contact forces
 * @author Scott Johnson
 * @param dt Timestep duration
 * @param domain Domain object to update
 */
void ExternalFaceManagerT::UpdateGeometricContactProperties(const realT dt,
                                                            PhysicalDomainT& domain,
                                                            Array1dT<Array1dT<R1Tensor> >& xs,
                                                            const bool updateFESoln)
{
  //-----------------------------
  // TO SET -->
  //-----------------------------
  //set by UpdateGeometricContactProperties sub function
  const rArray1d& area = this->GetFieldData<realT>("area");
  const rArray1d& areaC = domain.m_contactManager.GetFieldData<realT>("area");
  iArray1d& activeC = domain.m_contactManager.GetFieldData<int>("active"); //may be reset for DE-DE contact ... only reason it's not const
  const Array1dT<R1Tensor>& xpolyPts = domain.m_contactManager.m_intersectionPolygonPoints;
  const Array1dT<lArray1d>& contactToXPoly = domain.m_contactManager.m_contactToIntersectionPolygonPointsMap;
  const Array1dT<R1Tensor>& normalC = domain.m_contactManager.GetFieldData<R1Tensor>("normal");

  //set in this function
  rArray1d& normalApproach = this->GetFieldData<realT>("normalApproach");
  rArray1d& normalApproachMax = this->GetFieldData<realT>("normalApproachMax");
  rArray1d& normalApproachC = domain.m_contactManager.GetFieldData<realT>("normalApproach");
  rArray1d& normalApproachMaxC = domain.m_contactManager.GetFieldData<realT>("normalApproachMax");
  Array1dT<R1Tensor>& applicationPointC = domain.m_contactManager.GetFieldData<R1Tensor>("applicationPoint");
  Array1dT<R1Tensor>& face1ParentSoln = domain.m_contactManager.GetFieldData<R1Tensor>("face1ParentSoln");
  Array1dT<R1Tensor>& face2ParentSoln = domain.m_contactManager.GetFieldData<R1Tensor>("face2ParentSoln");
  Array1dT<R1Tensor>& velocityC = domain.m_contactManager.GetFieldData<R1Tensor>("velocity"); //velocityC = 0.0;

  //-----------------------------
  // ALREADY SET -->
  //-----------------------------
  const Array1dT<R1Tensor>& pos_de = domain.m_discreteElementManager.GetFieldData<
      FieldInfo::currentPosition>();
  const OrderedVariableOneToManyRelation& deFaceToNodes = domain.m_discreteElementSurfaceFaces.m_toNodesRelation;

  //-----------------------------
  //-----------------------------
  //-----------------------------

  Array1dT<Array1dT<R1Tensor> > nodalVelocities;
  const Array1dT<R1Tensor>& vel_fe = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();
  const Array1dT<R1Tensor>& vel_de = domain.m_discreteElementSurfaceNodes.GetFieldData<
      FieldInfo::velocity>();
  {
    nodalVelocities.resize(this->DataLengths());
    bool is_fe;
    for (localIndex i = 0; i < this->DataLengths(); ++i)
    {
      const localIndex faceIndex = this->FaceIndex(i, is_fe);
      if (is_fe)
        domain.m_feFaceManager.NodalPositions(domain.m_feNodeManager, faceIndex, xs[i]);
      else
        domain.m_discreteElementSurfaceFaces.NodalPositions(domain.m_discreteElementSurfaceNodes,
                                                            faceIndex, xs[i]);

      const lArray1d& faceToNodeMap =
          is_fe ? domain.m_feFaceManager.m_toNodesRelation[faceIndex] :
              deFaceToNodes[faceIndex];
      nodalVelocities[i].resize(faceToNodeMap.size());
      CopyGlobalToLocal(faceToNodeMap, is_fe ? vel_fe : vel_de, nodalVelocities[i]);
    }
  }

  //GET THE GEOMETRIC PROPERTIES
  UpdateGeometricContactPropertiesSub(dt, domain, xs);

  //---------------STEP 2: PHYSICS OF CONTACT----------------------
  bool is_fe1 = false, is_fe2 = false;

  //-----ITERATE NEIGHBORLIST AND UPDATE-----
  const Array1dT<R1Tensor>& nfc     = this->GetFieldData<R1Tensor> ( "faceNormal");
  localIndex index = 0;

  domain.m_discreteElementManager.m_discreteElementToDiscreteElementContactsMap.clear();

  const lArray1d& nodeToDE = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::demIndex>();
  for (Array1dT<lArray1d>::size_type ixfc1 = 0; ixfc1 < this->m_neighborList.size(); ++ixfc1)
  {
    const localIndex faceIndex1 = this->FaceIndex(ixfc1, is_fe1);
    for (lArray1d::size_type it = 0; it < this->m_neighborList[ixfc1].size(); ++it, ++index)
    {
      const localIndex ixfc2 = this->m_neighborList[ixfc1][it];
      const localIndex faceIndex2 = this->FaceIndex(ixfc2, is_fe2);
      if (index >= domain.m_contactManager.DataLengths())
        throw GPException("index incremented past the end of the contact manager fields!");
      if (IsAnyContact(activeC[index]))
      {
        //-----------------------------------------------------------
        //CONTACT 1: get the geometry for either mechanical or flow contact

        //REMOVED TO ACCOMMODATE DISTRIBUTED CONTACT ... SEE CommonPlaneInterferenceGeometry FOR DEFINITION
        if(!this->m_smoothedContact)
          GeometryUtilities::Centroid_3DPolygon(contactToXPoly[index], xpolyPts, applicationPointC[index]);
#if 0
        R1Tensor tmpPt;
        static int tmpIndex = 0;
        GeometryUtilities::Centroid_3DPolygon(contactToXPoly[index], xpolyPts, tmpPt);
        std::cout << (tmpIndex++) << "\napp_pt: " << applicationPointC[index](0) << " " << applicationPointC[index](1) << " " << applicationPointC[index](2) << "\n";
        std::cout << "tmp_pt: " << tmpPt(0) << " " << tmpPt(1) << " " << tmpPt(2) << "\n\n";
#endif

        // the applicationPoint is the centroid of the common plane
        R1Tensor cpFacePosition1 = static_cast<R1Tensor>(0.0);
        R1Tensor cpFaceVelocity1 = static_cast<R1Tensor>(0.0);
        if(is_fe1)
        {
          if(updateFESoln)
          {
            if(domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension==2)
            {
              //TODO: add in linear interpolation of velocity
              FacePositionAndVelocityFE(normalC[index], applicationPointC[index], xs[ixfc1],
                                        nodalVelocities[ixfc1], cpFacePosition1,
                                        cpFaceVelocity1);

            }
            else
            {
              FacePositionAndVelocityFE(normalC[index], applicationPointC[index], xs[ixfc1],
                                        nodalVelocities[ixfc1], this->m_tol.feParentSolution,
                                        face1ParentSoln[index], cpFacePosition1,
                                        cpFaceVelocity1);
            }
          }
          else
            FacePositionAndVelocityFE(xs[ixfc1], nodalVelocities[ixfc1],cpFacePosition1, cpFaceVelocity1);
        }
        else
        {
          localIndex ide1 =
              nodeToDE[deFaceToNodes[faceIndex1][0]];
          GeometryUtilities::ProjectPointToPlaneAlongUnitVector(applicationPointC[index], xs[ixfc1][0], nfc[ixfc1],
                                                                normalC[index], cpFacePosition1);
          domain.m_discreteElementManager.VelocityAtPoint(ide1, applicationPointC[index], cpFaceVelocity1);
        }

        R1Tensor cpFacePosition2 = static_cast<R1Tensor>(0.0);
        R1Tensor cpFaceVelocity2 = static_cast<R1Tensor>(0.0);
        if(is_fe2)
        {
          if(updateFESoln)
          {
            if(domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension==2)
            {
              //TODO: add in linear interpolation of velocity
              FacePositionAndVelocityFE(normalC[index], applicationPointC[index], xs[ixfc2],
                                        nodalVelocities[ixfc2], cpFacePosition2,
                                        cpFaceVelocity2);

            }
            else
            {
              FacePositionAndVelocityFE(normalC[index], applicationPointC[index], xs[ixfc2],
                                      nodalVelocities[ixfc2], this->m_tol.feParentSolution,
                                      face2ParentSoln[index], cpFacePosition2,
                                      cpFaceVelocity2);
            }
          }
          else
            FacePositionAndVelocityFE(xs[ixfc2], nodalVelocities[ixfc2],cpFacePosition2, cpFaceVelocity2);
        }
        else
        {
          localIndex ide2 =
              nodeToDE[deFaceToNodes[faceIndex2][0]];
          GeometryUtilities::ProjectPointToPlaneAlongUnitVector(applicationPointC[index], xs[ixfc2][0], nfc[ixfc2],
                                                                normalC[index], cpFacePosition2);
          domain.m_discreteElementManager.VelocityAtPoint(ide2, applicationPointC[index], cpFaceVelocity2);
        }

        R1Tensor cpDistance = cpFacePosition1;
        cpDistance -= cpFacePosition2;
        normalApproachC[index] = Dot(cpDistance, normalC[index]);
        if(normalApproachC[index] < 0 && IsMechanicalContact(activeC[index]))
        {
          //TODO: WE NEED TO AMEND THE CONTACT; THIS IS A PATHOLOGICAL CASE
          SetCommonPlaneGeometryAsOverlap(index, ixfc1, ixfc2, domain, xs);
        }
        if(normalApproachC[index] > normalApproachMaxC[index])
          normalApproachMaxC[index] = normalApproachC[index];

        if(!IsMechanicalContact(activeC[index]))
          continue;

        velocityC[index] = cpFaceVelocity2;
        velocityC[index] -= cpFaceVelocity1;

#ifndef DEEFC
        //-----------------------------------------------------------
        //CONTACT 3: ONLY FOR DE-DE CONTACT (activeC[index] == 1 or 2 && !is_fe1 && !is_fe2):
        //           try using Hertzian contact for DE-DE contact
        if (!is_fe1 && !is_fe2)
        {
          localIndex ide1 =
              nodeToDE[deFaceToNodes[faceIndex1][0]];
          localIndex ide2 =
              nodeToDE[deFaceToNodes[faceIndex2][0]];
          if (ide1 > ide2)
          {
            localIndex tmp = ide1;
            ide1 = ide2;
            ide2 = tmp;
          }

          //set discrete element contact contribution
          {
            R1Tensor tmp = cpFacePosition1;
            tmp -= pos_de[ide1];
            const realT r0 = tmp.L2_Norm();

            tmp = cpFacePosition2;
            tmp -= pos_de[ide2];
            const realT r1 = tmp.L2_Norm();

            //add the contribution to the contact
            domain.m_discreteElementManager.m_discreteElementToDiscreteElementContactsMap[ide1][ide2].AddContribution(
                areaC[index],
                normalApproachC[index],
                r0, r1,
                applicationPointC[index],
                normalC[index]);

            activeC[index] = demcontact;
          }
          continue;
        }
#endif
      }
    } //end foreach ixfc2
  } //end foreach ixfc1

  //-----------------------------------------------------------
  //CONTACT 5: redistribute the common plane states to the faces for evolution to the next step
  {
    normalApproach = 0;
    normalApproachMax = 0;
  }

  index = 0;
  for (Array1dT<lArray1d>::size_type ixfc1 = 0; ixfc1 < this->m_neighborList.size(); ++ixfc1)
  {
    for (lArray1d::size_type it = 0; it < this->m_neighborList[ixfc1].size(); ++it, ++index)
    {
      const localIndex ixfc2 = this->m_neighborList[ixfc1][it];

      if (index >= domain.m_contactManager.DataLengths())
        throw GPException("index incremented past the end of the contact manager fields!");

      if (IsMechanicalContact(activeC[index]))
      {
        //some states are homogenized based on the weighted average of the
        //properties from all of the common planes of which I am a part
        const realT alpha0 = areaC[index] / area[ixfc1];
        const realT alpha1 = areaC[index] / area[ixfc2];

        normalApproach[ixfc1] += normalApproachC[index] * alpha0;
        normalApproachMax[ixfc1] += normalApproachMaxC[index] * alpha0;

        normalApproach[ixfc2] += normalApproachC[index] * alpha1;
        normalApproachMax[ixfc2] += normalApproachMaxC[index] * alpha1;
      }
    }
  }
}

void ExternalFaceManagerT::SetCommonPlaneGeometryAsOverlap(const localIndex index,
                                                           const localIndex ixfc1,
                                                           const localIndex ixfc2,
                                                           PhysicalDomainT& domain,
                                                           const Array1dT<Array1dT<R1Tensor> >& xs)
{
  //TODO: FOR NOW WE LEAVE ALL IN xpolyPts, BUT WE SHOULD REMOVE THE DISCARDED POINTS AT SOME POINT

  //------------------------
  // TO RESET -->
  //------------------------
  Array1dT<R1Tensor>& xpolyPts = domain.m_contactManager.m_intersectionPolygonPoints;
  Array1dT<lArray1d>& contactToXPoly =
      domain.m_contactManager.m_contactToIntersectionPolygonPointsMap;
  contactToXPoly[index].clear();
  rArray1d& areaC = domain.m_contactManager.GetFieldData<realT>("area");
  rArray1d& area = this->GetFieldData<realT>("area");
  iArray1d& activeC = domain.m_contactManager.GetFieldData<int>("active");

  //------------------------
  //------------------------
  //------------------------
  const Array1dT<R1Tensor>& applicationPointC = domain.m_contactManager.GetFieldData<R1Tensor>("applicationPoint");
  const Array1dT<R1Tensor>& normalC = domain.m_contactManager.GetFieldData<R1Tensor>("normal");

  Array1dT<R1TensorT<2> > inPoints;
  R1Tensor e1, e2;
  {
    R1TensorT<2> min, max, tmin, tmax;
    Array1dT<R1TensorT<2> > xsl1, xsl2;
    Array1dT<R1Tensor> xs1(xs[ixfc1]), xs2(xs[ixfc2]);
    for(Array1dT<R1Tensor>::size_type i = 0; i < xs1.size(); i++)
      xs1[i] -= applicationPointC[index];
    for(Array1dT<R1Tensor>::size_type i = 0; i < xs2.size(); i++)
      xs2[i] -= applicationPointC[index];

    GeometryUtilities::VectorsInPlane(normalC[index], e1, e2);
    GeometryUtilities::CartesianPointsProjectedToPlanarPoints(xs1, e1, e2, xsl1, tmin, tmax);
    GeometryUtilities::CartesianPointsProjectedToPlanarPoints(xs2, e1, e2, xsl2, min, max);

    min.SetMin(tmin);
    max.SetMax(tmax);

    const realT minDim = (max(0) - min(0)) < (max(1) - min(1)) ?
        (max(0) - min(0)) :
        (max(1) - min(1));
    const realT positionTolerance = this->m_tol.spatial * minDim;
    //    const realT penetrationTolerance = this->tol.penetration * minDim;
    //    const realT areaTolerance = this->tol.area * minDim * minDim;

    if(IsMechanicalContact(activeC[index]))
    {
      area[ixfc1] -= areaC[index];
      area[ixfc2] -= areaC[index];
      activeC[index] = MechanicalToOpen();//activeC[index]);
    }

    areaC[index] = GeometryUtilities::Intersection_2DPolygons(xsl1, xsl2, positionTolerance, inPoints);
    GeometryUtilities::RemoveCoincidentOrderedPoints(inPoints, positionTolerance);
  }

  //set polygon points for future visualization
  contactToXPoly[index].reserve(inPoints.size());
  for(Array1dT<R1Tensor>::size_type ii = 0; ii < inPoints.size(); ++ii)
  {
    R1Tensor newPoint = e1;
    newPoint *= inPoints[ii](0);
    R1Tensor tmp = e2;
    tmp *= inPoints[ii](1);
    newPoint += tmp;
    newPoint += applicationPointC[index];
    contactToXPoly[index].push_back(xpolyPts.size());
    xpolyPts.push_back(newPoint);
  }
}

bool ExternalFaceManagerT::Criterion1(  const R1Tensor& nx1,
                                        const R1Tensor& nx2,
                                        const realT tolCosMin)
{
  return Dot(nx1, nx2) < -tolCosMin;
}

bool ExternalFaceManagerT::Criterion2(const R1Tensor& e1,
                                      const R1Tensor& e2,
                                      const Array1dT<R1Tensor>& xs1,
                                      const Array1dT<R1Tensor>& xs2,
                                      Array1dT<R1TensorT<2> >& xsl1,
                                      Array1dT<R1TensorT<2> >& xsl2,
                                      realT& minDim)
{
  R1TensorT<2> min1, min2, max1, max2;
  GeometryUtilities::CartesianPointsProjectedToPlanarPoints( xs1, e1, e2, xsl1, min1, max1);
  GeometryUtilities::CartesianPointsProjectedToPlanarPoints( xs2, e1, e2, xsl2, min2, max2);
  for(localIndex i = 0; i < 2; ++i)
    if(min1(i) > max2(i) || max1(i) < min2(i))
      return false;

  //Get the minimum dimension
  {
      minDim = max1(0) - min1(0);
      realT m = max1(1) - min1(1);
      if(m < minDim)
        minDim = m;
      m = max2(0) - min2(0);
      if(m < minDim)
        minDim = m;
      m = max2(1) - min2(1);
      if(m < minDim)
        minDim = m;
      if(isZero(minDim))
        return false;
  }
  return true;
}

bool ExternalFaceManagerT::Criterion3(const Array1dT<R1TensorT<2> >& xsl1,
                                      const Array1dT<R1TensorT<2> >& xsl2,
                                      Array1dT<R1TensorT<2> >& xsl,
                                      realT& area,
                                      const realT positionTolerance,
                                      const realT areaTolerance)
{
  xsl.clear();
  area = GeometryUtilities::Intersection_2DPolygons(xsl1, xsl2, positionTolerance, xsl);
  if(area < areaTolerance)
    return false;
  return true;
}


/**
 * @brief Determine the geometry (if it exists) for the interference of two faces on their common plane
 * @author Scott Johnson
 *
 *
 * @param[in] dt Timestep
 * @param[in] aper Aperture assigned to the contact
 * @param[in] xfc1 Center of the face
 * @param[in] dxfc1 Velocity of the face at the mid-step
 * @param[in] nx1 Face normal
 * @param[in] efs1 Structure to hold miscellaneous face properties for face 1
 * @param[in] xfc2 Center of the face
 * @param[in] dxfc2 Velocity of the face at the mid-step
 * @param[in] nx2 Face normal
 * @param[in] efs2 Structure to hold miscellaneous face properties for face 2
 * @param[in] tol Structure to hold tolerance values for contact
 * @param[out] centerCommonPlane Center of the shared common plane if in contact or the center of the geometry for flow contact
 * @param[out] normalCommonPlane Normal to the common plane or the geometry for flow contact
 * @param[out] areaCommonPlane Area of the common plane or the geometry for flow contact
 * @param[out] pointsCommonPlane Points in clockwise order formed by the interference of the faces on the common plane of for the interference geometry of the flow contact
 * @return 0 if no contact, 1 if mechanical contact, 2 if flow contact
 */
int ExternalFaceManagerT::CommonPlaneInterferenceGeometry(const R1Tensor& xfc1,
                                                          const R1Tensor& nx1,
                                                          const ExternalFaceStruct& efs1,

                                                          const R1Tensor& xfc2,
                                                          const R1Tensor& nx2,
                                                          const ExternalFaceStruct& efs2,

                                                          const ToleranceStruct& tol,

                                                          R1Tensor& centerCommonPlane,
                                                          R1Tensor& applicationPoint,
                                                          R1Tensor& normalCommonPlane,
                                                          realT& areaCommonPlane,
                                                          Array1dT<R1Tensor>& pointsCommonPlane)
{

  const int verbose = 0;

  //--------------------------------------------------------
  // CHECK FACE PAIR
  //--------------------------------------------------------

  //At this point, we only assure that the loose bounding spheres of the faces are contacting
  //Let's first make sure that the normals are opposite (within a tolerance) ...

  //+++++++++++++++++++++++++++++++
  //CRITERION 1: CHECK OPPOSITE POINTING NORMAL
  //note: -1 is perfect normal contact and 0 is orthogonal
  //      (positive is non-convex contact and inaccessible)
  //      i.e., if cosMin == 1, then the faces must be perfectly parallel
  if(!Criterion1(nx1, nx2, tol.cosMin))
  {
    if( verbose == 2 )
    {
      std::cout<<"  disqualification due to criterion 1 (opposite normals) "<<std::endl;
      std::cout<<"    xfc1 = "<<xfc1<<std::endl;
      std::cout<<"    xfc2 = "<<xfc2<<std::endl;
      std::cout<<"    nx1 = "<<nx1<<std::endl;
      std::cout<<"    nx2 = "<<nx2<<std::endl;
    }
    return nocontact;
  }

  //They are now oppositely oriented and somewhat "close"
  //Let's next make sure they are close enough before
  //going through all of the effort of detailed checking

  // estimate cp normal from face normals
  normalCommonPlane = nx1;
  normalCommonPlane -= nx2;
  normalCommonPlane.Normalize();

  //Calculate e1 and e2 (unit vectors in the common plane)
  //as well as lists of the projected points
  R1Tensor e1, e2;
  if(!GeometryUtilities::VectorsInPlane(normalCommonPlane, e1, e2))
    throw GPException("Normal to common plane is 0 ... cannot do anything with this!");

  //+++++++++++++++++++++++++++++++
  //CRITERION 2: CHECK FACES VERSUS COMMON PLANE PROJECTION
  // faces projected onto common plane should have overlapping
  // Plane-Aligned Bounding Boxes (PABB's) ... check that this is the case

  Array1dT<R1TensorT<2> > xsl1, xsl2, xsl;//hold the points projected onto the common plane
  realT minDim=0;
  if(!Criterion2(e1, e2, efs1.xs, efs2.xs, xsl1, xsl2, minDim))
  {
    if( verbose == 2 )
    {
      std::cout<<"  disqualification due to criterion 2 (PABB's not aligned) "<<std::endl;
      std::cout<<"    xfc1 = "<<xfc1<<std::endl;
      std::cout<<"    xfc2 = "<<xfc2<<std::endl;
      std::cout<<"    nx1 = "<<nx1<<std::endl;
      std::cout<<"    nx2 = "<<nx2<<std::endl;
    }
    return nocontact;
  }

  const realT positionTolerance = tol.spatial * minDim;
  const realT penetrationTolerance = tol.penetration * minDim;
  const realT areaTolerance = tol.area * minDim * minDim;

  //+++++++++++++++++++++++++++++++
  //CRITERION 3: MAKE SURE PROJECTIONS HAVE AREA
  // faces projected onto common plane should have a non-zero area
  realT xslArea = 0;
  if(!Criterion3(xsl1, xsl2, xsl, xslArea, positionTolerance, areaTolerance))
  {
    if( verbose == 2 )
    {
      std::cout<<"  disqualification due to criterion 3 (overlap of projection of faces has no area) "<<std::endl;
      std::cout<<"    xfc1 = "<<xfc1<<std::endl;
      std::cout<<"    xfc2 = "<<xfc2<<std::endl;
      std::cout<<"    nx1 = "<<nx1<<std::endl;
      std::cout<<"    nx2 = "<<nx2<<std::endl;
    }
    return nocontact;
  }


//  @ annavarapusr1: HACK FOR TETS TO WORK
//   When every face between two planes was fractured, the routines were predicting erroneous contact
//   This hack checks if the common plane area is equal to the face area - if not, it returns no contact
//  if(std::fabs(efs1.area - xslArea)>1e-14)
//  {
//    return nocontact;
//  }


  //Now, things get interesting; we need both the negative and positive overlaps
  //the negative is used for the flow solver, while the positive is for the mechanics solver
  //go through each point, and if it penetrates the common plane, it is part of the positive volume
  //otherwise the negative volume; identify intersection points
  //
  //To do any of this, we need the common plane center and the penetration distances
  //We are going to do the punishing operation of finding the overlap region's centroid
  //as the starting point of determining the actual common plane's center
  //
  // (note: if the faces are relatively close in size and location, we could replace this
  //        with an average of the face centers, but we may have significant size and
  //        location distribution, which precludes such an approach)
  R1TensorT<2> xcpl;
  {
    {
      GeometryUtilities::Centroid_2DPolygon(xsl, xcpl);
      GeometryUtilities::PlanarPointProjectedToCartesianCoordinates(xcpl, e1, e2, applicationPoint);
    }
    //get first attempt at common plane center
    {
      GeometryUtilities::ProjectPointToPlaneAlongUnitVector(applicationPoint,
                                                            xfc1, nx1,
                                                            normalCommonPlane,
                                                            centerCommonPlane);
      R1Tensor ctmp;
      GeometryUtilities::ProjectPointToPlaneAlongUnitVector(applicationPoint,
                                                            xfc2, nx2,
                                                            normalCommonPlane,
                                                            ctmp);
      centerCommonPlane += ctmp;
      centerCommonPlane *= 0.5;
    }
  }

  rArray1d penetration1(xsl1.size(),0.0); //positive into the common plane
  rArray1d penetration2(xsl2.size(),0.0); //positive into the common plane
  realT penetrationSummation1 = 0.0, penetrationSummation2 = 0.0;
  bool allOutOfRange = true;
  bool allOutOfContact = true;
  bool allInContact = true;
  {
    //    R1TensorT<2> xcpl;
    //    xcpl(0) = Dot(centerCommonPlane, e1);
    //    xcpl(1) = Dot(centerCommonPlane, e2);
    for(Array1dT<R1TensorT<2> >::size_type a = 0; a < xsl1.size(); ++a) {
      xsl1[a] -= xcpl;
      penetration1[a] = GeometryUtilities::DistanceToPlane(normalCommonPlane, centerCommonPlane, efs1.xs[a]);
      if(penetration1[a] > 0) {
        penetrationSummation1 += penetration1[a];
        allOutOfContact = false;
        if(penetration1[a] <= penetrationTolerance)
          allOutOfRange = false;
      }
      else if(penetration1[a] >= -tol.maximumSeparation) {
        allOutOfRange = false;
        allInContact = false;
      }
      else
      {
        allInContact = false;
      }
    }

    for(Array1dT<R1TensorT<2> >::size_type a = 0; a < xsl2.size(); ++a) {
      xsl2[a] -= xcpl;
      penetration2[a] = -GeometryUtilities::DistanceToPlane(normalCommonPlane, centerCommonPlane, efs2.xs[a]);
      if(penetration2[a] > 0) {
        penetrationSummation2 += penetration2[a];
        allOutOfContact = false;
        if(penetration2[a] <= penetrationTolerance)
          allOutOfRange = false;
      }
      else if(penetration2[a] >= -tol.maximumSeparation) {
        allOutOfRange = false;
        allInContact = false;
      }
      else
      {
        allInContact = false;
      }
    }
  }

  //+++++++++++++++++++++++++++++++
  //CRITERION 4: MAKE SURE AT LEAST SOME OF THE NODES ARE CLOSE ENOUGH TO THE COMMON PLANE
  if(allOutOfRange)
  {
    if( verbose == 2 )
    {
      std::cout<<"  disqualification due to criterion 4 (all nodes are too far from the common plane) "<<std::endl;
      std::cout<<"    xfc1 = "<<xfc1<<std::endl;
      std::cout<<"    xfc2 = "<<xfc2<<std::endl;
      std::cout<<"    nx1 = "<<nx1<<std::endl;
      std::cout<<"    nx2 = "<<nx2<<std::endl;
    }
    return nocontact;
  }

  //Faces are now oppositely oriented and quite "close"
  //We also have lists of the nodes projected to the common plane
  //coordinates where the 2D coordinates are relative to the candidate common plane center
  //We also know that the projection onto the common plane has finite area
  //for both faces, and we have just established the normal penetration distances on each face
  //and at least one of these normal distances is within contact range of the common plane
  //so ... we know we will either have opencontact or mechanicalcontact ... resolve which
  {
    //get detailed in and out of common plane geometry

    Array1dT<R1TensorT<2> > inContact1, inContact2;

    //if all of the nodes are in contact, there is no need to go through the expense of clipping
    //...there's nothing to be clipped, so just copy the old points to the new
    if(allInContact)
    {
      inContact1.resize(xsl1.size());
      inContact2.resize(xsl2.size());
      std::copy(xsl1.begin(), xsl1.end(), inContact1.begin());
      std::reverse_copy(xsl2.begin(), xsl2.end(), inContact2.begin());
    }
    else if(!allOutOfContact)
    {
      //here we have the case where some are in contact and some are not, so we need to do the clipping
      {
        Array1dT<R1TensorT<2> > outContact1;
        rArray1d inDistance1, outDistance1;
        GeometryUtilities::ClipByPlane_3DPolygon(normalCommonPlane,
                                                 centerCommonPlane,
                                                 efs1.xs, e1, e2,
                                                 penetration1,
                                                 inContact1, outContact1,
                                                 inDistance1, outDistance1,
                                                 false);
        GeometryUtilities::RemoveCoincidentOrderedPoints(inContact1, positionTolerance);
      }
      {
        Array1dT<R1TensorT<2> > outContact2;
        rArray1d inDistance2, outDistance2;
        GeometryUtilities::ClipByPlane_3DPolygon(normalCommonPlane,
                                                 centerCommonPlane,
                                                 efs2.xs, e1, e2,
                                                 penetration2,
                                                 inContact2, outContact2,
                                                 inDistance2, outDistance2,
                                                 true);
        GeometryUtilities::RemoveCoincidentOrderedPoints(inContact2, positionTolerance);
      }
    }

    //get detailed intersection and difference
    Array1dT<R1TensorT<2> > inPoints;
    areaCommonPlane = 0.0;
    if(inContact1.size() > 2 && inContact2.size() > 2)
    {
      //only look for the intersection of the polygons if they have enough points
      areaCommonPlane = GeometryUtilities::Intersection_2DPolygons(inContact1,
                                                                   inContact2,
                                                                   positionTolerance,
                                                                   inPoints);
      GeometryUtilities::RemoveCoincidentOrderedPoints(inPoints, positionTolerance);
    }

    //resolve contact geometry in global frame
    pointsCommonPlane.clear();

    //--NOW, DO THE FINAL DETERMINATION OF WHETHER THIS IS MECHANICAL OR OPEN--
    if(areaCommonPlane < areaTolerance || inPoints.size() < 3)
    {
      //THIS _MUST_ BE AN OPEN CONTACT, SO LET'S JUST USE THE OVERLAP REGION CALCULATED
      //FOR CRITERION 3: xsl, xcpl, etc.
      R1Tensor pt;
      for(Array1dT<R1TensorT<2> >::size_type i = 0; i < xsl.size(); ++i) {
        xsl[i] -= xcpl;//note: xsl now contains the 2D coordinates relative to the center of the region
        GeometryUtilities::PlanarPointProjectedToCartesianCoordinates(xsl[i], e1, e2, pt);
        pt += centerCommonPlane;
        pointsCommonPlane.push_back(pt);
      }
      return opencontact;
    }
    else
    {
      //You've now graduated ... contact is assured, so just resolve the common plane geometry
      R1TensorT<2> xcpNew;
      GeometryUtilities::Centroid_2DPolygon(inPoints, xcpNew);

      //note: since inPoints is relative to the old common plane center estimate, this is just an offset
      R1Tensor pt;
      GeometryUtilities::PlanarPointProjectedToCartesianCoordinates(xcpNew, e1, e2, pt);
      centerCommonPlane += pt;

      for(Array1dT<R1TensorT<2> >::size_type i = 0; i < inPoints.size(); ++i) {
        //translate the points to the newly determined common plane center
        inPoints[i] -= xcpNew;
        //transform to the planar coordinate system rotationally relative to the global frame
        GeometryUtilities::PlanarPointProjectedToCartesianCoordinates(inPoints[i], e1, e2, pt);
        //translate to the global frame
        pt += centerCommonPlane;
        pointsCommonPlane.push_back(pt);
      }

      return mechanicalcontact;
    }
  }
}

/**
 * @brief Determine the geometry (if it exists) for the interference of two 2D edges on their common edge
 * @author Scott Johnson
 *
 *
 * @param[in] dt Timestep
 * @param[in] aper Aperture assigned to the contact
 * @param[in] xfc1 Center of the edge
 * @param[in] dxfc1 Velocity of the edge at the mid-step
 * @param[in] nx1 Edge normal
 * @param[in] efs1 Structure to hold miscellaneous edge properties for edge 1
 * @param[in] xfc2 Center of the face
 * @param[in] dxfc2 Velocity of the face at the mid-step
 * @param[in] nx2 Face normal
 * @param[in] efs2 Structure to hold miscellaneous face properties for edge 2
 * @param[in] tol Structure to hold tolerance values for contact
 * @param[out] centerCommonEdge Center of the shared common plane if in contact or the center of the geometry for flow contact
 * @param[out] normalCommonEdge Normal to the common plane or the geometry for flow contact
 * @param[out] lengthCommonEdge Area of the common plane or the geometry for flow contact
 * @param[out] pointCommonEdge1 Terminus of common edge
 * @param[out] pointCommonEdge2 Terminus of common edge
 * @return 0 if no contact, 1 if mechanical contact, 2 if flow contact
 */
int ExternalFaceManagerT::CommonEdgeInterferenceGeometry(  const R1Tensor& xfc1,
                                                           const R1Tensor& nx1,
                                                           const ExternalFaceStruct& efs1,

                                                           const R1Tensor& xfc2,
                                                           const R1Tensor& nx2,
                                                           const ExternalFaceStruct& efs2,

                                                           const ToleranceStruct& tol,

                                                           R1Tensor& centerCommonEdge,
                                                           R1Tensor& normalCommonEdge,
                                                           realT& lengthCommonEdge,
                                                           Array1dT<R1Tensor>& pointsCommonEdge)
{
  pointsCommonEdge.clear();
  pointsCommonEdge.resize(2);

  const int verbose = 0;
  int ret = nocontact;

  //--------------------------------------------------------
  // CHECK EDGE PAIR
  //--------------------------------------------------------

  //At this point, we only assure that the loose bounding circles of the edges are contacting
  //Let's first make sure that the normals are opposite (within a tolerance) ...

  //+++++++++++++++++++++++++++++++
  //CRITERION 1: CHECK OPPOSITE POINTING NORMAL
  //note: -1 is perfect normal contact and 0 is orthogonal
  //      (positive is non-convex contact and inaccessible)
  //      i.e., if cosMin == 1, then the edges must be perfectly parallel
  if(!Criterion1(nx1, nx2, tol.cosMin))
  {
    if( verbose == 2 )
    {
      std::cout<<"  disqualification due to criterion 1 (opposite normals) "<<std::endl;
      std::cout<<"    xfc1 = "<<xfc1<<std::endl;
      std::cout<<"    xfc2 = "<<xfc2<<std::endl;
      std::cout<<"    nx1 = "<<nx1<<std::endl;
      std::cout<<"    nx2 = "<<nx2<<std::endl;
    }
    return ret;
  }

  // get common edge normal
  normalCommonEdge = nx1;
  normalCommonEdge -= nx2;
  normalCommonEdge.Normalize();

  // project points to the plane
  realT ndist[4];
  R1Tensor xplane(xfc1), zero(0.0);
  xplane += xfc2;
  xplane *= 0.5;

  R1Tensor t[4];
  ndist[0] = -GeometryUtilities::ProjectPointToPlane(efs1.xs[0], xplane, normalCommonEdge, t[0]);
  ndist[1] = -GeometryUtilities::ProjectPointToPlane(efs1.xs[1], xplane, normalCommonEdge, t[1]);
  ndist[2] = GeometryUtilities::ProjectPointToPlane(efs2.xs[0], xplane, normalCommonEdge, t[2]);
  ndist[3] = GeometryUtilities::ProjectPointToPlane(efs2.xs[1], xplane, normalCommonEdge, t[3]);
  //ndist is now (-) if the point is in the planar halfspace

  R1Tensor tmp(0.0);
  tmp(0) = normalCommonEdge(1);
  tmp(1) = -normalCommonEdge(0);

  realT dot[] = {Dot(t[0], tmp), Dot(t[1], tmp), Dot(t[2], tmp), Dot(t[3], tmp)};
  localIndex order[] = {0, 1, 2, 3};

  //put the edge points in order - efs1 & efs2 may not be ordered consistently!
  if(dot[0] > dot[1])
  {
    order[0] = 1;//first face, first point
    order[1] = 0;//first face, second point
  }
  if(dot[3] > dot[2])
  {
    order[2] = 3;//second face, first point
    order[3] = 2;//second face, second point
  }

  //The line segments don't overlap when projected to the common plane
//  if(dot[order[2] > dot[order[1]]] || dot[order[0]] > dot[order[3]])
//    return ret;
  if(fabs(dot[order[2]] - dot[order[1]]) > 1e-15 || fabs(dot[order[0]] - dot[order[3]]) > 1e-15)
    return ret;

  const bool isSmallerThanMachinePrecision =  (fabs(ndist[0])<1e-16 && fabs(ndist[1])<1e-16) || (fabs(ndist[2])<1e-16 && fabs(ndist[3])<1e-16);
  //The line segments do not both penetrate the common plane, so find the overlap
  const bool isOpen = (ndist[0] >= 0 && ndist[1] >=0) || (ndist[2] >= 0 && ndist[3] >=0) || isSmallerThanMachinePrecision;

  if(isOpen)
  {
    ret = opencontact;
  }
  else
  {
    ret = mechanicalcontact;
    if(ndist[0] * ndist[1] < 0)
    {
      if(ndist[0] < 0)
      {
        t[1] -= t[0];
        t[1] *= ndist[0] / (ndist[0] - ndist[1]);
        t[1] += t[0];
        ndist[1] = 0;
      }
      else
      {
        t[0] -= t[1];
        t[0] *= ndist[1] / (ndist[1] - ndist[0]);
        t[0] += t[1];
        ndist[0] = 0;
      }
    }
    if(ndist[2] * ndist[3] < 0)
    {
      if(ndist[2] < 0)
      {
        t[3] -= t[2];
        t[3] *= ndist[2] / (ndist[2] - ndist[3]);
        t[3] += t[2];
        ndist[3] = 0;
      }
      else
      {
        t[2] -= t[3];
        t[2] *= ndist[3] / (ndist[3] - ndist[2]);
        t[2] += t[3];
        ndist[2] = 0;
      }
    }
    for(localIndex i = 0; i < 4; i++)
      dot[i] = Dot(t[i], tmp);

    //put the edge points in order
    if(dot[0] > dot[1])
    {
      order[0] = 1;//first face, first point
      order[1] = 0;//first face, second point
    }
    if(dot[3] > dot[2])
    {
      order[2] = 3;//second face, first point
      order[3] = 2;//second face, second point
    }

    //The line segments don't overlap when projected to the common plane
    //order is 0:1, 3:2
//    if(dot[order[2] > dot[order[1]]] || dot[order[0]] > dot[order[3]])
//      ret = nocontact;
    if(fabs(dot[order[2]] - dot[order[1]]) > 1e-15 || fabs(dot[order[0]] - dot[order[3]]) > 1e-15)
      ret = nocontact;
  }

  const localIndex i0 = dot[order[0]] < dot[order[2]] ? order[2] : order[0];
  const localIndex i1 = dot[order[1]] < dot[order[3]] ? order[1] : order[3];

  pointsCommonEdge[0](0) = t[i0](0);
  pointsCommonEdge[0](1) = t[i0](1);
  pointsCommonEdge[0](2) = 0.0;

  pointsCommonEdge[1](0) = t[i1](0);
  pointsCommonEdge[1](1) = t[i1](1);
  pointsCommonEdge[1](2) = 0.0;

  tmp = pointsCommonEdge[1];
  tmp += pointsCommonEdge[0];
  tmp *= 0.5;
  GeometryUtilities::ProjectPointToPlane(centerCommonEdge, xplane, normalCommonEdge, centerCommonEdge);

  tmp = pointsCommonEdge[1];
  tmp -= pointsCommonEdge[0];
  lengthCommonEdge = tmp.L2_Norm();

  //TODO: check signs rather than fabs'ing it
  //May be able to move these checks up to save computational expense
  //SJ 1/15/14
  const realT minDim = fabs(dot[order[1]] - dot[order[0]]) < fabs(dot[order[3]] - dot[order[2]]) ? fabs(dot[order[1]] - dot[order[0]]) : fabs(dot[order[3]] - dot[order[2]]);
  const realT penetrationTolerance = tol.penetration * minDim;
  if(0.5*fabs(ndist[i1] + ndist[i0]) > penetrationTolerance)
    return nocontact;

  const realT lengthTolerance = tol.area * minDim;
  if(fabs(dot[i1] - dot[i0]) < lengthTolerance)
    return nocontact;

  return ret;
}




/**
 * @author settgast
 * @param nodeManager
 *
 * function to set "excludeFromContact" field based on the "NodeManagerT::nodeKCBC" field.
 */
void ExternalFaceManagerT::SetExcludeFromContact( const NodeManagerT& nodeManager, bool reset )
{

  if( nodeManager.HasField<int>("KCBC") )
  {
    iArray1d& excludeFromContact = this->GetFieldData<int>("excludeFromContact");
    if( reset )
      excludeFromContact = 1;

    const iArray1d& nodeKCBC = nodeManager.GetFieldData<int>("KCBC");

    for( localIndex ixfc=0 ; ixfc<this->DataLengths() ; ++ixfc )
    {
      bool isFE = false;
      const localIndex faceIndex = this->FaceIndex( ixfc, isFE );
      if( isFE )
      {
        for( lArray1d::const_iterator a=m_faceManager->m_toNodesRelation[faceIndex].begin() ; a!=m_faceManager->m_toNodesRelation[faceIndex].end() ; ++a )
        {
          if( nodeKCBC[*a] == 0 )
          {

            if( !reset && excludeFromContact[ixfc]!=0 )
              this->m_sorted = false;

            excludeFromContact[ixfc] = 0;
            break;
          }
        }
      }
    }
  }
}

void ExternalFaceManagerT::WriteSiloExternalFaces( SiloFile& siloFile,
                                                   const std::string& siloDirName,
                                                   const std::string& meshname,
                                                   const int centering,
                                                   const int cycleNum,
                                                   const realT problemTime,
                                                   const bool isRestart,
                                                   const std::string& regionName,
                                                   const lArray1d& mask )
{
  std::string subDirectory = siloDirName;
  std::string rootDirectory = "/" + siloDirName;
  siloFile.MakeSubDirectory( subDirectory, rootDirectory );
  DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

  //-----------------------------------------------------
  //HANDLE THE VARIABLES ASSOCIATED WITH THE JOINT STATE
  sArray1d intVarNames;
  sArray1d realVarNames;
  sArray1d R1TensorVarNames;
  sArray1d R2TensorVarNames;
  sArray1d R2SymTensorVarNames;

  Array1dT<iArray1d*> intVars;
  Array1dT<rArray1d*> realVars;
  Array1dT<Array1dT<R1Tensor>*> R1Vars;
  Array1dT<Array1dT<R2Tensor>*> R2Vars;
  Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

  m_contact->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateTemporaryFields( intVarNames, intVars );
  AllocateTemporaryFields( realVarNames, realVars );
  AllocateTemporaryFields( R1TensorVarNames, R1Vars );
  AllocateTemporaryFields( R2TensorVarNames, R2Vars );
  AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

  m_contact->Serialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

  ObjectDataStructureBaseT::WriteSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DeallocateTemporaryFields<int>(intVarNames);
  DeallocateTemporaryFields<realT>(realVarNames);
  DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
  DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
  DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

  WriteNonManagedDataMembersToSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DBSetDir(siloFile.m_dbFilePtr, "..");
}

void ExternalFaceManagerT::ReadSiloExternalFaces( const SiloFile& siloFile,
                                                  const std::string& siloDirName,
                                                  const std::string& meshname,
                                                  const int centering,
                                                  const int cycleNum,
                                                  const realT problemTime,
                                                  const bool isRestart,
                                                  const std::string& ,
                                                  const lArray1d& mask)

{
  if( DBSetDir(siloFile.m_dbFilePtr, siloDirName.c_str()) != -1 )
  {
    sArray1d intVarNames;
    sArray1d realVarNames;
    sArray1d R1TensorVarNames;
    sArray1d R2TensorVarNames;
    sArray1d R2SymTensorVarNames;

    Array1dT<iArray1d*> intVars;
    Array1dT<rArray1d*> realVars;
    Array1dT<Array1dT<R1Tensor>*> R1Vars;
    Array1dT<Array1dT<R2Tensor>*> R2Vars;
    Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

    m_contact->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

    AllocateTemporaryFields( intVarNames, intVars );
    AllocateTemporaryFields( realVarNames, realVars );
    AllocateTemporaryFields( R1TensorVarNames, R1Vars );
    AllocateTemporaryFields( R2TensorVarNames, R2Vars );
    AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

    ObjectDataStructureBaseT::ReadSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, "none", mask);

    ReadNonManagedDataMembersFromSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, "none", mask );

    m_contact->Deserialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

    DeallocateTemporaryFields<int>(intVarNames);
    DeallocateTemporaryFields<realT>(realVarNames);
    DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
    DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
    DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
}

void ExternalFaceManagerT::WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                                             const std::string& ,
                                                             const std::string& ,
                                                             const int ,
                                                             const int ,
                                                             const realT ,
                                                             const bool ,
                                                             const std::string& ,
                                                             const std::string& ,
                                                             const lArray1d& )
{
  siloFile.DBWriteWrapper("nfe",nfe);
  siloFile.DBWriteWrapper("nde",nde);

  int sorted = m_sorted;
  siloFile.DBWriteWrapper("sorted",sorted);
}

void ExternalFaceManagerT::ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                                              const std::string& ,
                                                              const std::string& ,
                                                              const int ,
                                                              const int ,
                                                              const realT ,
                                                              const bool ,
                                                              const std::string& ,
                                                              const lArray1d& )
{
  siloFile.DBReadWrapper("nfe",nfe);
  siloFile.DBReadWrapper("nde",nde);

  int sorted;
  siloFile.DBReadWrapper("sorted",sorted);
  m_sorted = sorted;
}

#ifdef SRC_INTERNAL
void ExternalFaceManagerT::GeodynCoupling( NodeManagerT& nodeManager )
{
  if (!doBackgroundAMR) return;
  if (!p_stepbamr) throw GPException("stepbamr pointer is not set...");
  A2Gtype A2G;

  const Array1dT<R1Tensor>& nodalReferencePosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& nodalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();
  const Array1dT<R1Tensor>& nodalVelocity = nodeManager.GetFieldData<FieldInfo::velocity>();
  Array1dT<R1Tensor>& nodalForce = nodeManager.GetFieldData<FieldInfo::force>();

  std::map<localIndex,localIndex> face,mirror;
  typedef std::map<localIndex,localIndex>::iterator faceit;

  localIndex numPackedNodes = 0;
  for( localIndex ixfc=0 ; ixfc<this->DataLengths() ; ++ixfc )
  {
    bool isFE = false;
    const localIndex faceIndex = FaceIndex(ixfc,isFE);
    if( !isFE ) continue;

    lArray1d& nodelist = m_faceManager->m_toNodesRelation[faceIndex];

    R1Tensor xmin= static_cast<R1Tensor>( std::numeric_limits<realT>::max() );
    R1Tensor xmax= xmin;
    xmax *= -1;

    for( localIndex a=0 ; a<nodelist.size() ; ++a )
    {
      xmin.SetMin(nodalReferencePosition[nodelist[a]]);
      xmax.SetMax(nodalReferencePosition[nodelist[a]]);
    }
    R1Tensor xdif=xmax-xmin;
    if(xdif[2]<std::numeric_limits<realT>::epsilon()) continue;
    localIndex n2n[2];
    int n=0;
    for( localIndex a=0 ; a<nodelist.size() ; ++a,++numPackedNodes )
    {
      if (nodalReferencePosition[nodelist[a]][2]<0) continue;
      n2n[n++]=a;
    }
    int ndif=n2n[1]-n2n[0];
    for( localIndex a=0 ; a<nodelist.size() ; ++a,++numPackedNodes )
    {
      if (nodalReferencePosition[nodelist[a]][2]<0) continue;
      n2n[n++]=a;
      localIndex ap=(a+1)%nodelist.size(), am=(a+nodelist.size()-1)%nodelist.size();
      if (nodalReferencePosition[nodelist[ap]][2]<0)
        mirror[nodelist[a]]=nodelist[ap];
      else
        mirror[nodelist[a]]=nodelist[am];
    }
    if (ndif==1) face[nodelist[n2n[1]]]=nodelist[n2n[0]];
    else         face[nodelist[n2n[0]]]=nodelist[n2n[1]];
  }

  A2G.reserve(face.size());
  int f(face.begin()->first),f0(f);
  static int dump_shape(1);
  std::ofstream ofshape;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(dump_shape && rank==0)
  {
    ofshape.open("bamr.shape.last");
    ofshape << face.size() << std::endl;
  }
  do
  {
    R1Tensor x(nodalReferencePosition[f]);
    x += nodalDisplacement[f];
    A2G.push_back(A2Gitem(f,x[0],x[1],nodalVelocity[f][0],nodalVelocity[f][1]));
    if(dump_shape && rank==0)
    {
      ofshape<<f<<" "<<x[0]<<" "<<x[1]<<" "<<nodalVelocity[f][0]<<" "<<nodalVelocity[f][1]<<std::endl;
    }
    faceit fit=face.find(f);
    f=fit->second;
    face.erase(fit);
  } while (f!=f0);
  dump_shape=0;
  G2Atype G2A;
  //  int res=p_stepbamr(A2G,G2A);
  assert ( G2A.size()==A2G.size() );
  for (G2Atype::iterator a(G2A.begin()),ae(G2A.end());a!=ae;++a)
  {
    localIndex li=a->first;
    R1Tensor nf(.5*a->second.fx,.5*a->second.fy,realT(0));
    nodalForce[li] += nf;
    nodalForce[mirror[li]] += nf;
    //    std::cout << li << nf <<std::endl;
  }
}

void ExternalFaceManagerT::GeodynCouplingParallel( NodeManagerT& nodeManager )
{
  if (!doBackgroundAMR)
    return;
  if (!p_stepbamr)
    throw GPException("stepbamr pointer is not set...");

  //Define temporary field arrays
  const Array1dT<R1Tensor>& nodalReferencePosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& nodalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();
  const Array1dT<R1Tensor>& nodalVelocity = nodeManager.GetFieldData<FieldInfo::velocity>();
  const Array1dT<R1Tensor>& faceNormals = this->GetFieldData<R1Tensor>("faceNormal");
  //  Array1dT<R1Tensor>& nodalForces = nodeManager.GetFieldData<FieldInfo::force>();

  //Define Bifroest
  GeodynBifroest share;
  share.Initialize();

  //Start filling A2G type structures
  TempBifroestNodeSendData sendNode;
  Array1dT<TempBifroestNodeSendData> keepNodes;
  TempBifroestFaceSendData sendFace;
  Array1dT<TempBifroestFaceSendData> keepFaces;
  {
    std::map<int, lSet> partitionNodes;
    std::map<int, lSet> partitionFaces;
    for( localIndex ixfc=this->IndexFEFirst() ; ixfc<this->IndexFEAfterLast() ; ++ixfc )
    {
      //go through all nodes associated with the face
      const lArray1d& nodelist = m_faceManager->m_toNodesRelation[ixfc];
      iSet partitions;
      const localIndex sz(nodelist.size());

      //set face properties
      bool is_fe;
      sendFace.faceIndex = this->m_faceManager->m_localToGlobalMap[this->FaceIndex(ixfc, is_fe)];
      sendFace.numberOfNodes = sz;
      sendFace.normal = faceNormals[ixfc];
      //set partitions
      {
        for( localIndex a=0 ; a<nodelist.size() ; ++a)
        {
          if(a > 3)
            throw GPException("Cannot yet handle faces with more than 4 nodes!");
          sendFace.nodeIndices[a] = nodeManager.m_localToGlobalMap[nodelist[a]];
          R1Tensor x(nodalReferencePosition[nodelist[a]]);
          //FIXME: add in GEODYN function to determine which partition the node belongs to
          partitions.insert(GeodynPartition(x));
        }
      }

      for(iSet::const_iterator iter = partitions.begin(); iter != partitions.end(); ++iter)
      {
        if(*iter == share.Rank())
          keepFaces.push_back(sendFace);
        else
          share.AddShare(*iter, sendFace);
        for( localIndex a=0 ; a<nodelist.size() ; ++a )
        {
          partitionNodes[*iter].insert(nodelist[a]);
        }
      }
    }

    //add partition nodes
    for(std::map<int, lSet>::const_iterator iter = partitionNodes.begin(); iter != partitionNodes.end(); ++iter)
    {
      const int irank = iter->first;
      if(irank == share.Rank())
      {
        for(lSet::const_iterator iter1 = iter->second.begin(); iter1 != iter->second.end(); ++iter1)
        {
          sendNode.nodeIndex = nodeManager.m_localToGlobalMap[*iter1];
          sendNode.x = nodalReferencePosition[*iter1];
          sendNode.x += nodalDisplacement[*iter1];
          sendNode.dx = nodalVelocity[*iter1];
          keepNodes.push_back(sendNode);
        }
      }
      else
      {
        for(lSet::const_iterator iter1 = iter->second.begin(); iter1 != iter->second.end(); ++iter1)
        {
          sendNode.nodeIndex = nodeManager.m_localToGlobalMap[*iter1];
          sendNode.x = nodalReferencePosition[*iter1];
          sendNode.x += nodalDisplacement[*iter1];
          sendNode.dx = nodalVelocity[*iter1];
          share.AddShare(irank, sendNode);
        }
      }
    }
  }//finished filling A2G type structures

  //Create communication topology for A2G and vice versa
  share.CoordinateShares();
  //FIXME: add logic here to deal with those nodes and faces not sent to other processes
  //SEE GEODYNBIFROEST.CPP AND CHANGE PROCESSPACKAGE AS APPROPRIATE TO DEAL WITH SLAVED OBJECTS

  //Communicate, process nodes, and deal with received information (i.e., G2A ...)
  //SEE GEODYNBIFROEST.CPP AND CHANGE SYNCHRONIZE AS APPROPRIATE TO DEAL WITH SLAVED OBJECTS
  share.Synchronize();
  for (Array1dT<TempBifroestFaceSendData>::const_iterator iter = keepFaces.begin(); iter != keepFaces.end(); ++iter)
  {
    //FIXME: add logic to update forces for the given keepFaces and keepNodes
    //globalIndex gi = iter->nodeIndices[0];
    //nodalForces[li] += nf;
    //nodalForces[mirror[li]] += nf;
  }
}
#endif

/*@annavarapusr1:
 * Get Weighting and stabilization parameters for Nitsche's method. Currently uses only analytical estimates
 * Not rigorous for quads and hexes - do an eigenvalue estimate if performance is poor. */
void ExternalFaceManagerT::GetProjectionTensorAndWeightingAndStabilizationParameters(const int dim,
                                                                                     const bool planeStress,
                                                                                     PhysicalDomainT& domain)
{
  const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

  rArray1d& nitscheStab_n  = this->GetFieldData<realT>("nitscheStab_n");
  rArray1d& nitscheStab_t1 = this->GetFieldData<realT>("nitscheStab_t1");
  rArray1d& nitscheStab_t2 = this->GetFieldData<realT>("nitscheStab_t2");
  rArray1d& nitscheGamma = this->GetFieldData<realT>("nitscheGamma");

  unsigned int totContFaces;
//  if(dim==2 or dim==3)
//  {
//    totContFaces = domain.m_feFaceManager.m_numFaces;
//  }
  if(dim==2)
  {
    totContFaces = domain.m_feFaceManager.m_numFaces;
  }
  else
  {
    totContFaces = domain.m_contactManager.DataLengths();
  }

//  unsigned int totContFaces = domain.m_contactManager.DataLengths();

  for (localIndex iContFace =  0; iContFace < totContFaces; iContFace++)
  {
    bool contActiv = false;
    localIndex kf1, kf2;

//    if(dim==2 or dim==3)
//    {
//      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//      if (!(childFaceIndex[iContFace].empty()))
//      {
//        contActiv = true;
////        kf1 = childFaceIndex[iContFace][0];
////        kf2 = childFaceIndex[iContFace][1];
//        if(dim==2)
//        {
//          kf1 = childFaceIndex[iContFace][0];
//          kf2 = childFaceIndex[iContFace][1];
//        }
//        if(dim==3)
//        {
//          kf1 = iContFace;
//          kf2 = childFaceIndex[iContFace][0];
//        }
//      }
//    }

    if(dim==2)
    {
      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

      if (!(childFaceIndex[iContFace].empty()))
      {
        contActiv = true;
        kf1 = childFaceIndex[iContFace][0];
        kf2 = childFaceIndex[iContFace][1];
      }
    }
    else
    {
      const iArray1d& active = domain.m_contactManager.GetFieldData<int>("activeInit");
      const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
      const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");

      if(active[iContFace]!=0)
      {
        contActiv = true;
        bool fe;
        kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
        kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
      }
    }

//    const iArray1d& active = domain.m_contactManager.GetFieldData<int>("active");
//    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//    const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//
//    if(active[iContFace]!=0)
//    {
//      contActiv = true;
//      bool fe;
//      kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
//      kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
//    }

    if(contActiv)
    {
      const realT measGam = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf1, true);
      realT measOmg1, measOmg2, NormD1, NormD2, stabAnalytical;
      GetMeasOmgAndNormD(dim, kf1, domain, planeStress, measOmg1, NormD1);
      GetMeasOmgAndNormD(dim, kf2, domain, planeStress, measOmg2, NormD2);

      if(domain.m_contactManager.m_weighted_nitsche_active)
      {
        nitscheGamma(faceToExternalFaceMap(kf1)) = NormD2*measOmg1/(NormD2*measOmg1+NormD1*measOmg2);
        nitscheGamma(faceToExternalFaceMap(kf2)) = NormD1*measOmg2/(NormD2*measOmg1+NormD1*measOmg2);
        stabAnalytical = (measGam*NormD1*NormD2)/(NormD2*measOmg1 + NormD1*measOmg2);
      }
      else
      {
        nitscheGamma(faceToExternalFaceMap(kf1)) = 0.5;
        nitscheGamma(faceToExternalFaceMap(kf2)) = 0.5;
        stabAnalytical = (measGam*0.25)*(NormD1*measOmg2 + NormD2*measOmg1)/(measOmg1*measOmg2);
      }

      // Analytical estimate might underestimate the parameter for quads/hexes.
      // Consequently multiply this estimate by a magnifying factor specified in the xml file.
      // Default value for this factor is unity.
      stabAnalytical = domain.m_contactManager.m_stab_magnifying_factor*stabAnalytical;

      nitscheStab_n(faceToExternalFaceMap(kf1)) = stabAnalytical;
      nitscheStab_n(faceToExternalFaceMap(kf2)) = stabAnalytical;

      nitscheStab_t1(faceToExternalFaceMap(kf1)) = stabAnalytical;
      nitscheStab_t1(faceToExternalFaceMap(kf2)) = stabAnalytical;

      if(dim==3)
      {
        nitscheStab_t2(faceToExternalFaceMap(kf1)) = stabAnalytical;
        nitscheStab_t2(faceToExternalFaceMap(kf2)) = stabAnalytical;
      }
    }
  }
}
//@author: annavarapusr1
//@brief: Function to calculate quantities required for the analytical estimate for the stabilization parameter
//        || D ||_2 is evaluated analytically to avoid calculating the maximum eigenvalue
void ExternalFaceManagerT::GetMeasOmgAndNormD(const int dim,
                                              const localIndex FaceID,
                                              const PhysicalDomainT& domain,
                                              const bool planeStress,
                                              realT& measOmg,
                                              realT& NormD)
{
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[FaceID][0].first;

  const rArray1d& elemVol = elemRegion->GetFieldData<FieldInfo::volume>();
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[FaceID][0].second;
  measOmg = elemVol(EleID);

  const localIndex paramIndex = elemRegion->m_mat->NumParameterIndex0() > 1 ? EleID : 0 ;
  const MaterialBaseParameterData& matParams = *(elemRegion->m_mat->ParameterData(paramIndex));

  const realT G = matParams.init_shearModulus;
  realT lambda = matParams.Lame;

  if(planeStress)
  {
    lambda = 2*lambda*G / ( lambda + 2*G );
  }

  rArray2d D(0.5*dim*(dim+1),0.5*dim*(dim+1));

  if(dim==2)
  {
    D(0,0) = lambda+2*G; D(0,1) = lambda; D(1,0) = lambda; D(1,1) = lambda+2*G; D(2,2) = G;

    if(planeStress)
    {
      realT L = 2*G*lambda/(2*G-lambda);
      NormD = (2*sqrt(4*pow(G,4) + 12*pow(G,3)*L +9*pow(G,2)*pow(L,2)))/(2*G + L);
    }
    else
    {
      NormD = 2*(G+lambda);
    }
  }
  if(dim==3)
  {
    D(0,0) = lambda+2*G; D(0,1) = lambda;     D(0,2) = lambda;
    D(1,0) = lambda;     D(1,1) = lambda+2*G; D(1,2) = lambda;
    D(2,0) = lambda;     D(2,1) = lambda;     D(2,2) = lambda+2*G;
    D(3,3) = G;          D(4,4) = G;          D(5,5) = G;
    NormD = 2*G+3*lambda;
  }
}
