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
      for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++ bcItr )
      {
        // check if the requested field has a wall boundary condition applied to it.
        WallBoundaryCondition* bc =  (*bcItr)->UpcastActiveBCPointer<WallBoundaryCondition>(0); // dynamic_cast<WallBoundaryCondition*> (*bcItr);
        if( bc )
        {
          Array1dT<R1Tensor>& v = object.GetFieldData<FieldInfo::velocity> ();
          const Array1dT<R1Tensor>& u = object.GetFieldData<FieldInfo::displacement> ();
          const Array1dT<R1Tensor>& X = object.GetFieldData<FieldInfo::referencePosition> ();

          R1Tensor x;
          R1Tensor vn;

          const R1Tensor& b = bc->m_position;
          const R1Tensor& n = bc->GetDirection(0.0); // fixme need time
          
          for(sArray1d::size_type i =0; i < bc->m_setNames.size(); ++i){
          	std::string setName = bc->m_setNames[i];
            
            std::map< std::string, lSet >::const_iterator setMap = object.m_Sets.find( setName );
            if( setMap != object.m_Sets.end() )
            {

              const lSet& set = setMap->second;

              for ( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a )
              {
                // get the projected position of the node at the end of the timestep
                x.cA( dt, v[*a] );
                x += X[*a];
                x += u[*a];

                // make "x" the vector from the known point on the plane to the end of step position.
               x -= b;
 
                realT residual = Dot( x, n );
 
                // if the residual is negative, then that means that the position is penetrated into the wall
                // by the distance defined by the residual.
                if( residual < 0.0 )
                {
                  // set the velocity so that the point ends up on the wall surface.
                  vn.cA( -residual / dt , n );
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
    for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=faceManager.m_bcData.begin() ; bcItr!=faceManager.m_bcData.end() ; ++ bcItr )
    {

      // check if the requested field has a boundary condition applied to it.
      TractionBoundaryCondition* bc = (*bcItr)->UpcastActiveBCPointer<TractionBoundaryCondition>(time);// dynamic_cast<TractionBoundaryCondition*> (*bcItr);
      if( bc )
      {

        Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force> ();

      //  const R1Tensor& n = bc->GetDirection(time); //oldway


        for(sArray1d::size_type i =0; i < bc->m_setNames.size(); ++i)
        {
          int findSet = 1;

          std::string setName = bc->m_setNames[i];
          std::map< std::string, lSet >::const_iterator setMap = faceManager.m_Sets.find( setName );
          if( setMap != faceManager.m_Sets.end() )
          {
          const lSet& set = setMap->second;

            for ( lSet::const_iterator k=set.begin() ; k!=set.end() ; ++k )
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
   * function to fill the KinematicConstrainNodes array with groups of nodes that are to be tied together by
   * a kinematic boundary condition.
   */
  void BuildKinematicConstraintBoundaryCondition( NodeManager& nodeManager,
                                                  FaceManagerT& faceManager,
                                                  Array1dT<lSet>& KinematicConstraintNodes,
                                                  const realT tiedNodeTolerance )
  {

    Array1dT<R1Tensor>& X = nodeManager.GetFieldData<FieldInfo::referencePosition>();
    const iArray1d& isExternal = nodeManager.m_isExternal;

    nodeManager.AddKeylessDataField<int>("KCBC",true,true);
    iArray1d& nodeKCBC = nodeManager.GetFieldData<int>("KCBC");
    nodeKCBC = 0;

    //***** first we will set up spatial bins to store nodes
    R1Tensor xmax, xmin, range;
    int numBinX = 2;
    int numBinY = 2;
    int numBinZ = 2;

    // set spatial min/max and range of the problem.
    for( Array1dT<R1Tensor>::const_iterator iterX=X.begin() ; iterX!=X.end() ; ++iterX )
    {
      xmin.SetMin(*iterX);
      xmax.SetMax(*iterX);
    }
    range  = xmax;
    range -= xmin;


    Array1dT<lSet> allBins(numBinX*numBinY*numBinZ);

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

          for( localIndex a=0; a<nodeManager.m_numNodes ; ++a )
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


    //***** now we need to loop over all nodes in a bin, and check if two or more within the tolerance to be tied
    //      together.

    // set a tolerance that scales with the size of the problem... TODO this should be something a little more logical
    const realT tol = Dot(range,range) * tiedNodeTolerance * tiedNodeTolerance ;
    R1Tensor diff;
    lSet alreadyTied;

    // loop over each bin.
    for( Array1dT<lSet>::const_iterator ibin=allBins.begin() ; ibin!=allBins.end() ; ++ibin )
    {
      // loop over each node in the bin, from beginning to end.
      for( lSet::const_iterator a=ibin->begin() ; a!=ibin->end() ; ++a )
      {
        // check to see if the node has already been tied. you will see how this is possible in a second.
        if( alreadyTied.count(*a) == 0 )
        {
          // hold a set of nodes that stores all indices that have been tied together.
          lSet tiedNodes;

          // now we set a new iterator "b" to the next bin node iterator "a"
          lSet::const_iterator b = a; ++b;

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
//                std::cout<<"Dot(diff("<<*a<<","<<*b<<")) = "<<Dot(diff,diff)<<" "<<tol<<std::endl;

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

                // add both nodes to the "alreadyTied" set. This is used to skip nodes that have alreay been placed into
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



    // set the number of nodes with kinematic constrain boundary conditions for each face.
    iArray1d& numBCNodesInFace = faceManager.GetFieldData<int>("numKCBCnodesInFace");
    const iArray1d& isExternalFace = faceManager.m_isExternal;
    numBCNodesInFace = 0;

    for( Array1dT<lSet>::iterator iterGroup=KinematicConstraintNodes.begin() ;
         iterGroup!=KinematicConstraintNodes.end() ; ++iterGroup )
    {
      for( lSet::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
      {
        for( lSet::const_iterator kf=nodeManager.m_nodeToFaceMap[*a].begin() ;
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
                                                  Array1dT<lSet>& KinematicConstraintNodes,
                                                  const realT tiedNodeNormalRuptureStress,
                                                  const realT tiedNodeShearRuptureStress )
  {

    const rArray1d& mass = nodeManager.GetFieldData<FieldInfo::mass>();
    Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
    Array1dT<R1Tensor>& velocity = nodeManager.GetFieldData<FieldInfo::velocity>();
    const iArray1d& isExternalFace = faceManager.m_isExternal;
    iArray1d& numKCBCnodesInFace = faceManager.GetFieldData<int>("numKCBCnodesInFace");

    iArray1d& nodeKCBC = nodeManager.GetFieldData<int>("KCBC");

    // loop over all groups of tied nodes
    for( Array1dT<lSet>::iterator iterGroup=KinematicConstraintNodes.begin() ;
        iterGroup!=KinematicConstraintNodes.end() ; ++iterGroup )
    {

      // variables to hold the sum of conserved quantities
      realT sumMass = 0.0;
      R1Tensor sumForce(0.0);
      R1Tensor sumMomentum(0.0);

      // set to hold list of broken nodes
      lSet brokenNodes;

      // loop over all nodes in the group of tied nodes
      for( lSet::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
      {
        // flag to see if the node should be broken based on topology. it is set to true by default, and changed to
        // false by the following face loop.
        bool breakNode = true;

        realT areaish = 0.0;
        R1Tensor normalish;

        // loop over all faces connected to this node
        for( lSet::const_iterator faceIndex=nodeManager.m_nodeToFaceMap[*a].begin() ;
             faceIndex!=nodeManager.m_nodeToFaceMap[*a].end() ; ++faceIndex )
        {
          // of course, we only care about external faces
          if( isExternalFace[*faceIndex] )
          {

            // effective "area" of the node.
            const int numNodesInFace = faceManager.m_toNodesRelation[*faceIndex].size();
            R1Tensor normal;
            R1Tensor center;
            areaish += faceManager.FaceCenterAndNormal( nodeManager, *faceIndex, center , normal ) / numNodesInFace;

            normalish += normal;

            // if the face has more than 3 nodes that are tied with KCBC, then don't split the node.
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
          // the nodes are not released, and the summation of the conserved quantities is calculated
          sumMass += mass[*a];
          sumForce += force[*a];

          R1Tensor momentum  = velocity[*a];
          momentum *= mass[*a];
          sumMomentum += momentum;
        }

      }







      // erase the nodes that are broken from the group
      for( lSet::const_iterator b=brokenNodes.begin() ; b!=brokenNodes.end() ; ++b )
      {

        for( lSet::const_iterator faceIndex=nodeManager.m_nodeToFaceMap[*b].begin() ;
             faceIndex!=nodeManager.m_nodeToFaceMap[*b].end() ; ++faceIndex )
        {
          if( isExternalFace[*faceIndex] )
          {
            --(numKCBCnodesInFace[*faceIndex]) ;
          }
        }
        nodeKCBC[*b] = 0;

        iterGroup->erase(*b);

      }

      // if there are still more than 1 node in the group, then we have to distribute the conserved quantities among
      // the members of the tied group.
      if( iterGroup->size() > 1 )
      {
        // momentum / mass = velocity
        sumMomentum /= sumMass;

        // now set each node in the group to value consistent with kinematic constrain... All nodes will have same velocity,
        // and acceleration.
        for( lSet::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
        {
          velocity[*a] = sumMomentum;
          force[*a]  = sumForce;
          force[*a] *= mass[*a]/sumMass;
        }
      }
    }


    // now delete tied groups with less than 2 nodes in the group.
    for( Array1dT<lSet>::size_type k=0 ; k<KinematicConstraintNodes.size() ; ++k )
    {
      lSet& group = KinematicConstraintNodes[k];
      Array1dT<lSet>::iterator iterGroup = KinematicConstraintNodes.begin() + k;

      if( group.size() < 2 )
      {

        for( lSet::const_iterator a=iterGroup->begin() ; a!=iterGroup->end() ; ++a )
        {
          for( lSet::const_iterator faceIndex=nodeManager.m_nodeToFaceMap[*a].begin() ;
               faceIndex!=nodeManager.m_nodeToFaceMap[*a].end() ; ++faceIndex )
          {
            if( isExternalFace[*faceIndex] )
            {
              --(numKCBCnodesInFace[*faceIndex]) ;
            }
          }
          nodeKCBC[*a] = 0;

        }


        KinematicConstraintNodes.erase(iterGroup);
      }
    }

  }

}
