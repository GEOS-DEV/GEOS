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
 * @file ApplyBoundaryConditions.cpp
 * @author settgast1
 * @date Apr 21, 2011
 */

#include "ApplyBoundaryConditions.h"
#include "BoundaryConditions.h"

namespace BoundaryConditionFunctions
{
void ApplyRigidWallBoundaryCondition( ObjectDataStructureBaseT& object, const realT dt )
{

  if( dt>0.0 )
  {

    // iterate over all boundary conditions.
    for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
    {
      // check if the requested field has a wall boundary condition applied to
      // it.
      WallBoundaryCondition* bc =  (*bcItr)->UpcastActiveBCPointer<WallBoundaryCondition>(0);   // dynamic_cast<WallBoundaryCondition*>
                                                                                                // (*bcItr);
      if( bc )
      {
        array<R1Tensor>& v = object.GetFieldData<FieldInfo::velocity> ();
        const array<R1Tensor>& u = object.GetFieldData<FieldInfo::displacement> ();
        const array<R1Tensor>& X = object.GetFieldData<FieldInfo::referencePosition> ();

        R1Tensor x;
        R1Tensor vn;

        const R1Tensor& b = bc->m_position;
        const R1Tensor& n = bc->GetDirection(0.0);   // fixme need time

        for(array<string>::size_type i =0 ; i < bc->m_setNames.size() ; ++i)
        {
          std::string setName = bc->m_setNames[i];

          std::map< std::string, set<localIndex> >::const_iterator setMap = object.m_Sets.find( setName );
          if( setMap != object.m_Sets.end() )
          {

            const set<localIndex>& set = setMap->second;

            for ( set<localIndex>::const_iterator a=set.begin() ; a!=set.end() ; ++a )
            {
              // get the projected position of the node at the end of the
              // timestep
              x.cA( dt, v[*a] );
              x += X[*a];
              x += u[*a];

              // make "x" the vector from the known point on the plane to the
              // end of step position.
              x -= b;

              realT residual = Dot( x, n );

              // if the residual is negative, then that means that the position
              // is penetrated into the wall
              // by the distance defined by the residual.
              if( residual < 0.0 )
              {
                // set the velocity so that the point ends up on the wall
                // surface.
                vn.cA( -residual / dt, n );
                v[*a] += vn;
              }
            }
          }
        }
      }
    }
  }
}



void ApplyTractionBoundaryCondition( PhysicalDomainT& domain, realT time)
{

  FaceManagerT& faceManager = domain.m_feFaceManager;
  NodeManager& nodeManager = domain.m_feNodeManager;

  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=faceManager.m_bcData.begin() ; bcItr!=faceManager.m_bcData.end() ; ++bcItr )
  {

    // check if the requested field has a boundary condition applied to it.
    TractionBoundaryCondition* bc = (*bcItr)->UpcastActiveBCPointer<TractionBoundaryCondition>(time);  // dynamic_cast<TractionBoundaryCondition*>
                                                                                                       // (*bcItr);
    if( bc )
    {

      array<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force> ();

      //  const R1Tensor& n = bc->GetDirection(time); //oldway


      for(array<string>::size_type i =0 ; i < bc->m_setNames.size() ; ++i)
      {
        int findSet = 1;

        std::string setName = bc->m_setNames[i];
        std::map< std::string, set<localIndex> >::const_iterator setMap = faceManager.m_Sets.find( setName );
        if( setMap != faceManager.m_Sets.end() )
        {
          const set<localIndex>& set = setMap->second;

          for ( set<localIndex>::const_iterator k=set.begin() ; k!=set.end() ; ++k )
          {
            if (domain.m_feFaceManager.m_toElementsRelation[*k].size()>0)
            {
              const int numNodesOnFace = faceManager.m_toNodesRelation[*k].size();
              const realT area = faceManager.SurfaceArea( nodeManager, *k );

              /* Old way
                 realT value = bc->GetValue(faceManager,k,time);

                 R1Tensor traction = n;
                 if( bc->IsNormalTraction() )
                 {
                 traction = faceManager.FaceNormal( nodeManager, *k );
                 }
                 traction *= value * area / numNodesOnFace;
               */
              R1Tensor traction = bc->GetTractionOnFace(domain,k,time);
              traction *= area / numNodesOnFace;

              for( lArray1d::const_iterator a=faceManager.m_toNodesRelation[*k].begin() ;
                   a!=faceManager.m_toNodesRelation[*k].end() ; ++a )
              {
                force[*a] += traction;
              }

            }


          }
        }
        else
        {
          findSet = 0;
        }

        int findSetAll;
        MPI_Allreduce(&findSet, &findSetAll, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (findSetAll == 0)
        {
          std::stringstream msg;
          msg << "Traction boundary condition is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
          throw GPException(msg.str());
        }


      }

    }
  }
}



/**
 * @author settgast
 * @param nodeManager
 * @param KinematicConstraintNodes
 * @param tiedNodeTolerance
 * function to fill the KinematicConstrainNodes array with groups of nodes that
 * are to be tied together by
 * a kinematic boundary condition.
 */
void BuildKinematicConstraintBoundaryCondition( NodeManager& nodeManager,
                                                FaceManagerT& faceManager,
                                                array<set<localIndex>>& KinematicConstraintNodes,
                                                const realT tiedNodeTolerance )
{

  array<R1Tensor>& X = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const array<integer>& isExternal = nodeManager.m_isExternal;

  nodeManager.AddKeylessDataField<int>("KCBC",true,true);
  array<integer>& nodeKCBC = nodeManager.GetFieldData<int>("KCBC");
  nodeKCBC = 0;

  //***** first we will set up spatial bins to store nodes
  R1Tensor xmax, xmin, range;
  int numBinX = 2;
  int numBinY = 2;
  int numBinZ = 2;

  // set spatial min/max and range of the problem.
  for( array<R1Tensor>::const_iterator iterX=X.begin() ; iterX!=X.end() ; ++iterX )
  {
    xmin.SetMin(*iterX);
    xmax.SetMax(*iterX);
  }
  range  = xmax;
  range -= xmin;


  array<set<localIndex>> allBins(numBinX*numBinY*numBinZ);

  // fill the effing bins with node indexes that belong in that bin
  for( int i=0 ; i<numBinX ; ++i )
  {
    for( int j=0 ; j<numBinY ; ++j )
    {
      for( int k=0 ; k<numBinZ ; ++k )
      {
        R1Tensor min, max, size;
        min[0] = xmin[0] + range[0] * i/numBinX;
        max[0] = xmin[0] + range[0] * (i+1)/numBinX;
        min[1] = xmin[1] + range[1] * j/numBinY;
        max[1] = xmin[1] + range[1] * (j+1)/numBinY;
        min[2] = xmin[2] + range[2] * k/numBinZ;
        max[2] = xmin[2] + range[2] * (k+1)/numBinZ;
        size  = max;
        size -= min;
        size *= 0.01;
        min -= size;
        max += size;

        for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
        {
          if( isExternal[a] == 1 )
          {
            if( X[a][0]>=min[0] && X[a][0]<=max[0] &&
                X[a][1]>=min[1] && X[a][1]<=max[1] &&
                X[a][2]>=min[2] && X[a][2]<=max[2] )
            {
              allBins[i*numBinY*numBinZ + j*numBinZ + k].insert(a);
            }
          }
        }
      }
    }
  }


  //***** now we need to loop over all nodes in a bin, and check if two or more
  // within the tolerance to be tied
  //      together.

  // set a tolerance that scales with the size of the problem... TODO this
  // should be something a little more logical
  const realT tol = Dot(range,range) * tiedNodeTolerance * tiedNodeTolerance;
  R1Tensor diff;
  set<localIndex> alreadyTied;

  // loop over each bin.
  for( array<set<localIndex>>::const_iterator ibin=allBins.begin() ; ibin!=allBins.end() ; ++ibin )
  {
    // loop over each node in the bin, from beginning to end.
    for( set<localIndex>::const_iterator a=ibin->begin() ; a!=ibin->end() ; ++a )
    {
      // check to see if the node has already been tied. you will see how this
      // is possible in a second.
      if( alreadyTied.count(*a) == 0 )
      {
        // hold a set of nodes that stores all indices that have been tied
        // together.
        set<localIndex> tiedNodes;

        // now we set a new iterator "b" to the next bin node iterator "a"
        set<localIndex>::const_iterator b = a; ++b;

        // now we loop over all the remaining nodes in the ibin.
        for( ; b!=ibin->end() ; ++b )
        {
          // now check to see if this node has already been tied in a group.
          if( alreadyTied.count(*b) == 0 )
          {
            // calculate the distance between the nodes
            diff  = X[*a];
            diff -= X[*b];


//              if( Dot(diff,diff) < 0.01 )
//                std::cout<<"Dot(diff("<<*a<<","<<*b<<")) =
// "<<Dot(diff,diff)<<" "<<tol<<std::endl;

            // check to see if the distance is lower than the tolerance.
            if( Dot(diff,diff) < tol )
            {
              // Tie the nodes together!!

              // set the position to the position of the first node.
              X[*b] = X[*a];
//                std::cout<<"MATCH FOUND"<<std::endl;

              // add both nodes to the "tiedNodes" set.
              tiedNodes.insert(*a);
              tiedNodes.insert(*b);

              nodeKCBC[*a] = 1;
              nodeKCBC[*b] = 1;

              // add both nodes to the "alreadyTied" set. This is used to skip
              // nodes that have alreay been placed into
              // a tied node set.
              alreadyTied.insert(*a);
              alreadyTied.insert(*b);
            }
          }
        }
        // add the tiedNodes set to the KinematicConstraintNodes array.
        if( !(tiedNodes.empty()) )
          KinematicConstraintNodes.push_back(tiedNodes);
      }
    }
  }



  // set the number of nodes with kinematic constrain boundary conditions for
  // each face.
  array<integer>& numBCNodesInFace = faceManager.GetFieldData<int>("numKCBCnodesInFace");
  const array<integer>& isExternalFace = faceManager.m_isExternal;
  numBCNodesInFace = 0;

  for( array<set<localIndex>>::iterator iterGroup=KinematicConstraintNodes.begin() ;
       iterGroup!=KinematicConstraintNodes.end() ; ++iterGroup )
  {
    for( set<localIndex>::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
    {
      for( set<localIndex>::const_iterator kf=nodeManager.m_nodeToFaceMap[*a].begin() ;
           kf!=nodeManager.m_nodeToFaceMap[*a].end() ; ++kf )
      {
        if( isExternalFace[*kf] )
        {
          ++(numBCNodesInFace[*kf]);
        }
      }
    }
  }



}

/**
 *
 * @param faceManager
 * @param nodeManager
 * @param KinematicConstraintNodes
 * @param tiedNodeRuptureStress
 */
void ApplyKinematicConstraintBoundaryCondition( FaceManagerT& faceManager,
                                                NodeManager& nodeManager,
                                                array<set<localIndex>>& KinematicConstraintNodes,
                                                const realT tiedNodeNormalRuptureStress,
                                                const realT tiedNodeShearRuptureStress )
{

  const array<real64>& mass = nodeManager.GetFieldData<FieldInfo::mass>();
  array<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
  array<R1Tensor>& velocity = nodeManager.GetFieldData<FieldInfo::velocity>();
  const array<integer>& isExternalFace = faceManager.m_isExternal;
  array<integer>& numKCBCnodesInFace = faceManager.GetFieldData<int>("numKCBCnodesInFace");

  array<integer>& nodeKCBC = nodeManager.GetFieldData<int>("KCBC");

  // loop over all groups of tied nodes
  for( array<set<localIndex>>::iterator iterGroup=KinematicConstraintNodes.begin() ;
       iterGroup!=KinematicConstraintNodes.end() ; ++iterGroup )
  {

    // variables to hold the sum of conserved quantities
    realT sumMass = 0.0;
    R1Tensor sumForce(0.0);
    R1Tensor sumMomentum(0.0);

    // set to hold list of broken nodes
    set<localIndex> brokenNodes;

    // loop over all nodes in the group of tied nodes
    for( set<localIndex>::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
    {
      // flag to see if the node should be broken based on topology. it is set
      // to true by default, and changed to
      // false by the following face loop.
      bool breakNode = true;

      realT areaish = 0.0;
      R1Tensor normalish;

      // loop over all faces connected to this node
      for( set<localIndex>::const_iterator faceIndex=nodeManager.m_nodeToFaceMap[*a].begin() ;
           faceIndex!=nodeManager.m_nodeToFaceMap[*a].end() ; ++faceIndex )
      {
        // of course, we only care about external faces
        if( isExternalFace[*faceIndex] )
        {

          // effective "area" of the node.
          const int numNodesInFace = faceManager.m_toNodesRelation[*faceIndex].size();
          R1Tensor normal;
          R1Tensor center;
          areaish += faceManager.FaceCenterAndNormal( nodeManager, *faceIndex, center, normal ) / numNodesInFace;

          normalish += normal;

          // if the face has more than 3 nodes that are tied with KCBC, then
          // don't split the node.
          if( numKCBCnodesInFace[*faceIndex] == numNodesInFace )
          {
            breakNode = false;
          }
        }
      }

      normalish.Normalize();

      R1Tensor stressTypeThingie = force[*a];
      stressTypeThingie /= -areaish;

      const realT normalStressish = Dot( normalish, stressTypeThingie );
      R1Tensor junk;
      junk = normalish;
      junk *= -normalStressish;
      junk += stressTypeThingie;
      const realT shearStressish = junk.L2_Norm();



      // evaluate whether to release the nodes from their tied group.
      if( breakNode || ( normalStressish > tiedNodeNormalRuptureStress ) || (shearStressish>tiedNodeShearRuptureStress) )
      {
        // insert the node into the broken nodes list.
        brokenNodes.insert(*a);

      }
      else
      {
        // the nodes are not released, and the summation of the conserved
        // quantities is calculated
        sumMass += mass[*a];
        sumForce += force[*a];

        R1Tensor momentum  = velocity[*a];
        momentum *= mass[*a];
        sumMomentum += momentum;
      }

    }



    // erase the nodes that are broken from the group
    for( set<localIndex>::const_iterator b=brokenNodes.begin() ; b!=brokenNodes.end() ; ++b )
    {

      for( set<localIndex>::const_iterator faceIndex=nodeManager.m_nodeToFaceMap[*b].begin() ;
           faceIndex!=nodeManager.m_nodeToFaceMap[*b].end() ; ++faceIndex )
      {
        if( isExternalFace[*faceIndex] )
        {
          --(numKCBCnodesInFace[*faceIndex]);
        }
      }
      nodeKCBC[*b] = 0;

      iterGroup->erase(*b);

    }

    // if there are still more than 1 node in the group, then we have to
    // distribute the conserved quantities among
    // the members of the tied group.
    if( iterGroup->size() > 1 )
    {
      // momentum / mass = velocity
      sumMomentum /= sumMass;

      // now set each node in the group to value consistent with kinematic
      // constrain... All nodes will have same velocity,
      // and acceleration.
      for( set<localIndex>::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
      {
        velocity[*a] = sumMomentum;
        force[*a]  = sumForce;
        force[*a] *= mass[*a]/sumMass;
      }
    }
  }


  // now delete tied groups with less than 2 nodes in the group.
  for( array<set<localIndex>>::size_type k=0 ; k<KinematicConstraintNodes.size() ; ++k )
  {
    set<localIndex>& group = KinematicConstraintNodes[k];
    array<set<localIndex>>::iterator iterGroup = KinematicConstraintNodes.begin() + k;

    if( group.size() < 2 )
    {

      for( set<localIndex>::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
      {
        for( set<localIndex>::const_iterator faceIndex=nodeManager.m_nodeToFaceMap[*a].begin() ;
             faceIndex!=nodeManager.m_nodeToFaceMap[*a].end() ; ++faceIndex )
        {
          if( isExternalFace[*faceIndex] )
          {
            --(numKCBCnodesInFace[*faceIndex]);
          }
        }
        nodeKCBC[*a] = 0;

      }


      KinematicConstraintNodes.erase(iterGroup);
    }
  }

}

}
