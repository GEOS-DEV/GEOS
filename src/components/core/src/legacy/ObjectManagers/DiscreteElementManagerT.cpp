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

/*
 * DiscreteElementManagerT.cpp
 *
 *  Created on: Sep 20, 2010
 *      Author: settgast1
 */

#include "DiscreteElementManagerT.h"
#include "Utilities/GeometryUtilities.h"

/**
 * @brief Constructor to set internal pointers to the external node and face
 * managers
 * @author Scott Johnson
 *
 * @param nm External node manager pointer
 * @param fm External face manager pointer
 */
DiscreteElementManagerT::DiscreteElementManagerT(NodeManager* nm = 0, FaceManagerT* fm = 0):
  DiscreteElementManagerBaseT(ObjectDataStructureBaseT::DiscreteElementManager),
  m_discreteElementToExternalFacesMap(m_VariableOneToManyMaps["discreteElementToExternalFaces"]),
  m_discreteElementToExternalNodesMap(m_VariableOneToManyMaps["discreteElementToExternalNodes"]),
  m_discreteElementToDiscreteElementContactsMap()
{
  this->m_faceManager = fm;
  this->m_nodeManager = nm;

  //add fields for rotation, translation, and tractions
  this->AddBaseFields();

  //add fields to nodemanager - used during reading
  this->m_nodeManager->AddKeyedDataField<FieldInfo::relativePosition>();
  this->m_nodeManager->AddKeyedDataField<FieldInfo::currentPosition>();

  //needed to resolve face geometry and dynamics for contact
  this->m_nodeManager->AddKeyedDataField<FieldInfo::velocity>();
  this->m_nodeManager->AddKeyedDataField<FieldInfo::acceleration>();
  this->m_nodeManager->AddKeyedDataField<FieldInfo::force>();
  this->m_nodeManager->AddKeyedDataField<FieldInfo::contactForce>();
  this->m_nodeManager->AddKeyedDataField<FieldInfo::mass>();

  //needed for mapping between discrete element manager and dependent face and
  // node managers
  this->m_nodeManager->AddKeyedDataField<FieldInfo::demIndex>();
  this->m_faceManager->AddKeyedDataField<FieldInfo::demIndex>();
}

DiscreteElementManagerT::~DiscreteElementManagerT()
{
  // TODO Auto-generated destructor stub
}

void DiscreteElementManagerT::WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                                                const std::string& siloDirName,
                                                                const std::string& meshname,
                                                                const int centering,
                                                                const int cycleNum,
                                                                const realT problemTime,
                                                                const bool isRestart,
                                                                const std::string& multiRoot,
                                                                const std::string& regionName,
                                                                const lArray1d& mask)
{
  localIndex ncontacts = 0;
  std::string s;

  if(m_discreteElementToDiscreteElementContactsMap.size() == 0)
    return;

  std::string subDirectory =   "DiscreteElementContacts";
  if ( DBMkDir(siloFile.m_dbFilePtr, subDirectory.c_str()) != 0)
    throw GPException("Cannot make directory DiscreteElementContacts");
  DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());
  for(std::map<localIndex, std::map<localIndex, DiscreteElementContact> >::const_iterator it0 =
        m_discreteElementToDiscreteElementContactsMap.begin() ;
      it0 != m_discreteElementToDiscreteElementContactsMap.end() ; ++it0)
  {
    for(std::map<localIndex, DiscreteElementContact>::const_iterator it1 = it0->second.begin() ;
        it1 != it0->second.end() ; ++it1, ++ncontacts)
    {
      std::stringstream ss;
      ss << ncontacts;
      s = ss.str();
      if ( DBMkDir(siloFile.m_dbFilePtr, s.c_str()) != 0)
        throw GPException("Cannot make directory for individual contact under DiscreteElementContacts");
      DBSetDir(siloFile.m_dbFilePtr, s.c_str());
      siloFile.DBWriteWrapper("i0",it0->first);
      siloFile.DBWriteWrapper("i1",it1->first);
      it1->second.WriteSilo(siloFile);
      DBSetDir(siloFile.m_dbFilePtr, "..");
    }
  }
  siloFile.DBWriteWrapper("n_discreteElementToDiscreteElementContactsMap",ncontacts);
  DBSetDir(siloFile.m_dbFilePtr, "..");
}

void DiscreteElementManagerT::ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                                                 const std::string& siloDirName,
                                                                 const std::string& meshname,
                                                                 const int centering,
                                                                 const int cycleNum,
                                                                 const realT problemTime,
                                                                 const bool isRestart,
                                                                 const std::string& regionName,
                                                                 const lArray1d& mask)
{
  localIndex ncontacts, i0, i1;
  std::string s;

  std::string subDirectory =   "DiscreteElementContacts";
  if ( DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str()) == 0)
  {
    siloFile.DBReadWrapper("n_discreteElementToDiscreteElementContactsMap",ncontacts);
    for(localIndex i = 0 ; i < ncontacts ; i++)
    {
      std::stringstream ss;
      ss << i;
      s = ss.str();
      if ( DBSetDir(siloFile.m_dbFilePtr, s.c_str()) == 0)
      {
        siloFile.DBReadWrapper("i0",i0);
        siloFile.DBReadWrapper("i1",i1);
        m_discreteElementToDiscreteElementContactsMap[i0][i1].ReadSilo(siloFile);
        DBSetDir(siloFile.m_dbFilePtr, "..");
      }
      else
      {
        throw GPException("Cannot find requested contact index");
      }
    }
    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
}


unsigned int
DiscreteElementManagerT::Unpack( const char*& buffer, lArray1d& elementReceiveLocalIndices )
{
  unsigned int sizeOfUnpacked = DiscreteElementManagerBaseT::Unpack(buffer, elementReceiveLocalIndices);

  //now also pack the number of nodes and faces associated with each DE
  int itmp;
  globalIndex gtmp;
  lArray1d& nodeToDE = this->m_nodeManager->GetFieldData<FieldInfo::demIndex>();
  lArray1d& faceToDE = this->m_faceManager->GetFieldData<FieldInfo::demIndex>();
  std::map<globalIndex,localIndex>::const_iterator gtol;
  for(lArray1d::const_iterator it = elementReceiveLocalIndices.begin() ; it != elementReceiveLocalIndices.end() ; ++it)
  {
    sizeOfUnpacked += bufvector::Unpack( buffer, itmp );
    for(int i = 0 ; i < itmp ; ++i)
    {
      sizeOfUnpacked += bufvector::Unpack( buffer, gtmp );
      gtol = this->m_nodeManager->m_globalToLocalMap.find(gtmp);
      if(gtol == this->m_nodeManager->m_globalToLocalMap.end())
        throw GPException("DiscreteElementManagerT::Unpack : cannot find the requested global node id to reference");
      this->m_discreteElementToExternalNodesMap[*it].push_back(gtol->second);
      nodeToDE(gtol->second) = *it;
    }

    sizeOfUnpacked += bufvector::Unpack( buffer, itmp );
    for(int i = 0 ; i < itmp ; ++i)
    {
      sizeOfUnpacked += bufvector::Unpack( buffer, gtmp );
      gtol = this->m_faceManager->m_globalToLocalMap.find(gtmp);
      if(gtol == this->m_faceManager->m_globalToLocalMap.end())
        throw GPException("DiscreteElementManagerT::Unpack : cannot find the requested global face id to reference");
      this->m_discreteElementToExternalFacesMap[*it].push_back(gtol->second);
      faceToDE(gtol->second) = *it;
    }
  }
  return sizeOfUnpacked;
}

/**
 * @author Scott Johnson
 * @param sendElements local indices of elements to pack and send
 * @param buffer the buffer to pack the elements into
 * @return size of characters packed into the buffer.
 *
 * This function packs complete elements into a buffer. this should include all
 * information needed to reconstruct the element
 * on a remote process domain. This does not include maps to other objects, as
 * those are locally indexed relations and
 * must be constructed on the receiving domain.
 */
unsigned int
DiscreteElementManagerT::Pack( const lArray1d& sendElements, bufvector& buffer ) const
{
  unsigned int sizeOfPacked = DiscreteElementManagerBaseT::Pack(sendElements, buffer);

  //now also pack the number of nodes and faces associated with each DE
  int itmp;
  globalIndex gtmp;
  for(lArray1d::const_iterator it = sendElements.begin() ; it != sendElements.end() ; ++it)
  {
    lArray1d& nodeIndices = this->m_discreteElementToExternalNodesMap(*it);
    itmp = nodeIndices.size();
    sizeOfPacked += buffer.Pack( itmp );
    for(lArray1d::const_iterator iter = nodeIndices.begin() ; iter != nodeIndices.end() ; ++iter)
    {
      gtmp = this->m_nodeManager->m_localToGlobalMap(*iter);
      sizeOfPacked += buffer.Pack( gtmp );
    }

    lArray1d& faceIndices = this->m_discreteElementToExternalFacesMap(*it);
    itmp = faceIndices.size();
    sizeOfPacked += buffer.Pack( itmp );
    for(lArray1d::const_iterator iter = faceIndices.begin() ; iter != faceIndices.end() ; ++iter)
    {
      gtmp = this->m_faceManager->m_localToGlobalMap(*iter);
      sizeOfPacked += buffer.Pack( gtmp );
    }
  }
  return sizeOfPacked;
}

/**
 * @brief Apply a given force in the global frame at the given nodal point to
 * the discrete element
 * @author Scott Johnson
 * @param[in] nodeIndex nodal index in the node manager
 * @param[in] f global frame force vector to apply to the discrete element
 * associated with the given node
 */
void DiscreteElementManagerT::ApplyNodalForce(const localIndex nodeIndex, const R1Tensor& f)
{
  const localIndex ide = this->m_nodeManager->GetFieldData<FieldInfo::demIndex>()[nodeIndex];

  //(1) increment force
  this->GetFieldData<FieldInfo::force> ()[ide] += f;

  //(2) increment moment
  R1Tensor dm;
  {
    //get the local frame force
    R1Tensor fl;
    GlobalToLocalDirection(ide, f, fl);

    //get the local frame moment arm
    const R1Tensor& relativePosition =
      this->m_nodeManager->GetFieldData<FieldInfo::relativePosition> ()[nodeIndex];

    //get the incremental moment
    dm.Cross(relativePosition, fl);
  }
  this->GetFieldData<FieldInfo::moment> ()[ide] += dm;
}

void DiscreteElementManagerT::CalculateEnergy()
{
  const array<real64>& mass = GetFieldData<FieldInfo::mass>();
  const array<R1Tensor>& vels = GetFieldData<FieldInfo::velocity> ();
  const array<R1Tensor>& rvels = GetFieldData<FieldInfo::rotationalVelocity> ();
  const array<R1Tensor>& inertias = GetFieldData<FieldInfo::rotationalInertia> ();

  array<real64>& energy = GetFieldData<realT> ("kinetic_energy");
  for(localIndex a = 0 ; a < this->DataLengths() ; a++)
  {
    R1Tensor omega2(rvels[a]);
    for(localIndex j = 0 ; j < nsdof ; j++)
      omega2(j) *= omega2(j);
    energy[a] = mass[a] * Dot(vels[a],vels[a]) + Dot(inertias[a],omega2);
    energy[a] *= 0.5;
  }
}

void DiscreteElementManagerT::ApplyDiscreteElementContactForces(const realT youngs, const realT poissons)
{
  std::map<localIndex, std::map<localIndex, DiscreteElementContact> >::const_iterator iter0;
  std::map<localIndex, DiscreteElementContact>::const_iterator iter1;

  array<R1Tensor>& deMoment    = GetFieldData<FieldInfo::moment> ();
  deMoment = 0.0;
  array<R1Tensor>& deForce     = GetFieldData<FieldInfo::force> ();
  deForce = 0.0;

  const array<R1Tensor>& dePosition  = GetFieldData<FieldInfo::currentPosition> ();

  R2Tensor r1, r2;
  R1Tensor m, ml, frc, tmp, tmp1;

  const realT hsc = (2.0/3.0)*(youngs/(1-poissons*poissons));

  for(iter0 = m_discreteElementToDiscreteElementContactsMap.begin() ; iter0 != m_discreteElementToDiscreteElementContactsMap.end() ; ++iter0)
  {
    this->RotationTensor(iter0->first, r1);
    for(iter1 = iter0->second.begin() ; iter1 != iter0->second.end() ; ++iter1)
    {
      this->RotationTensor(iter1->first, r2);

      const DiscreteElementContact& dec = iter1->second;

      //get contact values: force (frc) and contact point (cp)
      R1Tensor cp;
      {
        realT dn = 0.0, reff = 0.0;
        if(!dec.Values(dn, reff, cp, frc))
          continue;

        const realT fmag = -hsc * dn * sqrt(reff * dn);
        frc *= fmag;
      }

      //FIXME: this only works in serial ... think about parallel!!!

      //first
      {
        const localIndex a = iter0->first;
        deForce[a] += frc;
        tmp1 = cp;
        tmp1 -= dePosition[a];
        m.Cross(tmp1, frc);
        ml.AijBj(r1, m);
        deMoment[a] += ml;
      }

      //second
      {
        const localIndex a = iter1->first;
        deForce[a] -= frc;
        tmp1 = cp;
        tmp1 -= dePosition[a];
        m.Cross(tmp1, frc);
        ml.AijBj(r1, m);
        deMoment[a] -= ml;
      }
    }
  }
  CalculateEnergy();
}

/**
 * @brief Apply a given force in the global frame for the given face of the
 * discrete element
 * @author Scott Johnson
 * @param[in] faceIndex Local face index in the face manager
 * @param[in] x global frame application point of the force
 * @param[in] f global frame force vector to apply to the discrete element
 * associated with the given face
 */
void DiscreteElementManagerT::ApplyFaceForce(const localIndex faceIndex,
                                             const R1Tensor& x,
                                             const R1Tensor& f)
{
  lArray1d& faceToDE = this->m_faceManager->GetFieldData<FieldInfo::demIndex>();
  const localIndex ide = faceToDE(faceIndex);

  //(1) increment force
  this->GetFieldData<FieldInfo::force> ()[ide] += f;

  //(2) increment moment
  R1Tensor dm;
  R1Tensor r(x);
  r -= this->GetFieldData<FieldInfo::currentPosition>()[ide];
  //get global frame moment
  dm.Cross(r, f);
  {
    //get the local frame moment
    R1Tensor ml;
    GlobalToLocalDirection(ide, dm, ml);
    this->GetFieldData<FieldInfo::moment> ()[ide] += ml;
  }
}

/**
 * @brief Apply all nodal forces to the parent discrete element
 * @author Scott Johnson
 */
void DiscreteElementManagerT::ApplyNodalForces()
{
  if(this->DataLengths()==0)
    return;

  //for contact resolution of faces
  const array<R1Tensor>& nodalForce = m_nodeManager->GetFieldData<FieldInfo::force> ();
  const array<R1Tensor>& relativePosition = m_nodeManager->GetFieldData<FieldInfo::relativePosition> ();

  array<R1Tensor>& deMoment               = GetFieldData<FieldInfo::moment> ();
  deMoment = 0.0;
  array<R1Tensor>& deForce                = GetFieldData<FieldInfo::force> ();
  deForce = 0.0;

  R2Tensor r;
  R1Tensor m, fl;

  //--FOR ENERGY CALCULATION--
  const array<real64>& mass = GetFieldData<FieldInfo::mass>();
  const array<R1Tensor>& vels = GetFieldData<FieldInfo::velocity> ();
  const array<R1Tensor>& rvels = GetFieldData<FieldInfo::rotationalVelocity> ();
  const array<R1Tensor>& inertias = GetFieldData<FieldInfo::rotationalInertia> ();
  array<real64>& energy = GetFieldData<realT> ("kinetic_energy");

  for (localIndex a = 0 ; a < this->DataLengths() ; ++a)
  {
    //get the rotation from global to local frame
    this->RotationTensor(a, r);

    for(localIndex b=0 ; b<m_discreteElementToExternalNodesMap[a].size() ; b++)
    {
      //nodal index
      const localIndex ni = m_discreteElementToExternalNodesMap[a][b];

      //get contribution to body force
      deForce[a] += nodalForce[ni];

      //get contribution to moment
      fl.AijBj(r, nodalForce[ni]);
      m.Cross(relativePosition[ni], fl);
      deMoment[a] += m;
    }

    //--FOR ENERGY CALCULATION--
    R1Tensor omega2(rvels[a]);
    for(localIndex j = 0 ; j < nsdof ; j++)
      omega2(j) *= omega2(j);
    energy[a] = mass[a] * Dot(vels[a],vels[a]) + Dot(inertias[a],omega2);
    energy[a] *= 0.5;
  }
}

/**
 * @brief Recalculates several physical properties given the geometry
 * @author Scott Johnson
 * Recalculates the nodal relative positions
 * Recalculates the de reference and current positions, the rotation, the
 * rotational inertia, and the mass
 */
void DiscreteElementManagerT::RecalculatePhysicalProperties()
{
  //mass, volume, position
  array<real64>& mass = GetFieldData<FieldInfo::mass> ();
  array<R1Tensor>& pos = GetFieldData<FieldInfo::currentPosition> ();
  array<R1Tensor>& rpos = GetFieldData<FieldInfo::referencePosition> ();

  //rotational inertia
  array<R1Tensor>& I = GetFieldData<FieldInfo::rotationalInertia> ();

  //nodal position and mass
  array<R1Tensor>& nodeRelativePosition = m_nodeManager->GetFieldData<
    FieldInfo::relativePosition> ();
  array<R1Tensor>& nodeCurrentPosition =
    m_nodeManager->GetFieldData<FieldInfo::currentPosition> ();
  array<real64>& nodeMass =
    m_nodeManager->GetFieldData<FieldInfo::mass> ();

  //go through each DE
  for (localIndex idem = 0 ; idem != this->DataLengths() ; ++idem)
  {
    //----------------------------------------------------------------------
    // 1) Find the volume, centroid, and moment of inertia in space frame
    //    a) Set the mass + the current and reference positions

    //find the number of triangular facets and fill the facet indices
    localIndex tmax = 0;
    lArray1d index;
    realT totalSurfaceArea = 0.;
    {
      const lArray1d& externalFaces = this->m_discreteElementToExternalFacesMap[idem];
      for (lArray1d::const_iterator ite = externalFaces.begin() ; ite != externalFaces.end() ; ++ite)
      {
        const lArray1d& tarr = this->m_faceManager->m_toNodesRelation[*ite];

        R1Tensor center, normal;
        const realT faceArea = GeometryUtilities::Centroid_3DPolygon(tarr, nodeCurrentPosition, center, normal );
        totalSurfaceArea += faceArea;

        //get triangular facets in the index array
        {
          const localIndex ttmp = tarr.size() - 2;
          tmax += ttmp;

          //get the triangular facet index list as well as the facial area
          for(localIndex inode = 0 ; inode < ttmp ; inode++)
          {
            //get the triangular facet indices
            index.push_back(tarr[0]);
            index.push_back(tarr[inode + 1]);
            index.push_back(tarr[inode + 2]);
          }
        }

        //get the facial area and allocate DE mass to nodes based on
        // sum(faceArea_i/nnodes_in_face_i)_node_j
        {
          //go through every node in the face and add in the (faceArea) /
          // (number of faces in the node)
          const realT ttmp = faceArea / tarr.size();
          for(lArray1d::const_iterator it = tarr.begin() ; it != tarr.end() ; ++it)
            nodeMass[*it] += ttmp;
        }
      }
    }

    //get the physical properties
    R2SymTensor inertia;
    CalculatePhysicalProperties(nodeCurrentPosition, tmax, index, mass[idem], pos[idem], inertia);

    //right now, the nodeMass holds the total surface area attributed to a node
    //to convert to equivalent mass we need to multiply by mass / total area
    // (see below)
    const realT nodalMassFactor = mass[idem] / totalSurfaceArea;

    //fill the de reference position with the current position
    rpos[idem] = pos[idem];
    {
      R2Tensor rotation;
      realT lambda[nsdof];

      // 2) Find the Eigenvectors of the moment of inertia
      inertia.EigenVals(lambda, 0.);
      R1Tensor v[nsdof];
      inertia.EigenVecs(lambda, v);

      // 3) Resolve (and set) the moment of inertia in local frame
      for (localIndex i = 0 ; i < nsdof ; i++)
      {
        I[idem](i) = lambda[i];
        for (localIndex j = 0 ; j < nsdof ; ++j)
        {
          // 4) Get the rotation based on the Eigenvectors
          rotation(i, j) = v[i](j);
        }
      }
      this->SetRotation(idem, rotation);
      this->RotationTensor(idem, rotation);

      // 5) Find the nodal relative positions
      const lArray1d& tmp2 = this->m_discreteElementToExternalNodesMap[idem];
      for (lArray1d::const_iterator it = tmp2.begin() ; it != tmp2.end() ; ++it)
      {
        //GlobalToLocal(idem, nodeCurrentPosition[ii],
        // nodeRelativePosition[ii]);
        //since you're doing the same operation repeatedly, just do the
        // operation
        R1Tensor tmppos(nodeCurrentPosition[*it]);
        tmppos -= pos[idem];
        nodeRelativePosition[*it].AijBj(rotation, tmppos);
        nodeMass[*it] *= nodalMassFactor;
      }
    }
  }
}

/**
 * @brief Update nodal states and zero out the discrete element forces and
 * accelerations
 * @author Scott Johnson
 */
void DiscreteElementManagerT::UpdateNodalStatesZeroForcesAndAccelerations()
{
  //zero the DE accelerations and forces
  array<R1Tensor>& acceleration = this->GetFieldData<FieldInfo::acceleration> ();
  array<R1Tensor>& force = this->GetFieldData<FieldInfo::force> ();

  //zero the DE rotation accelerations and moments
  array<R1Tensor>& racc = this->GetFieldData<FieldInfo::rotationalAcceleration> ();
  array<R1Tensor>& moment = this->GetFieldData<FieldInfo::moment> ();

  //for contact resolution of faces
  array<R1Tensor>& nodeVelocity       = m_nodeManager->GetFieldData<FieldInfo::velocity> ();
  array<R1Tensor>& nodeAcceleration   = m_nodeManager->GetFieldData<FieldInfo::acceleration> ();
  array<R1Tensor>& nodeForce          = m_nodeManager->GetFieldData<FieldInfo::force> ();
  array<R1Tensor>& nodeDisplacement   = m_nodeManager->GetFieldData<FieldInfo::displacement> ();

  const array<R1Tensor>& nodeCurrentPosition   = m_nodeManager->GetFieldData<FieldInfo::currentPosition> ();
  const array<R1Tensor>& nodeReferencePosition = m_nodeManager->GetFieldData<FieldInfo::referencePosition> ();

  R2Tensor r;
  acceleration = 0.0;
  force = 0.0;
  racc = 0.0;
  moment = 0.0;
  nodeAcceleration = 0.0;
  nodeForce = 0.0;

  for (localIndex a = 0 ; a < DataLengths() ; ++a)
  {
    this->RotationTensor(a, r);
    for(localIndex b=0 ; b<m_discreteElementToExternalNodesMap[a].size() ; b++)
    {
      int nn = m_discreteElementToExternalNodesMap[a][b];
      nodeDisplacement[nn] = nodeCurrentPosition[nn];
      nodeDisplacement[nn] -= nodeReferencePosition[nn];
      NodalVelocity(a,b,r,nodeVelocity[nn]);
    }
  }
}
