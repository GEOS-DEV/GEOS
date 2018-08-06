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
 * @file ApplyBoundaryConditions.h
 * @author settgast1
 * @date Apr 21, 2011
 */

#ifndef APPLYBOUNDARYCONDITIONS_H_
#define APPLYBOUNDARYCONDITIONS_H_

#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ObjectManagers/TableManager.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "BoundaryConditions.h"
#include <mpi.h>

namespace BoundaryConditionFunctions
{

/////////////////////////////////////////////////
// Forward declaration of templated functions
template< typename T >
inline void SetValue( const lSet& set,
                      ObjectDataStructureBaseT& object,
                      BoundaryConditionBase* bc,
                      array<T>& fieldData,
                      realT time );

template< typename T >
void ApplyDirichletBoundaryCondition( ObjectDataStructureBaseT& object,
                                      const std::string fieldName,
                                      realT time );


template< typename T, typename Solver, typename BCFunctionPtr>
void ApplyBoundaryCondition(Solver* solverPtr,
                            BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const std::string& fieldName, realT time );


template< typename T, typename Solver, typename BCFunctionPtr>
void ApplyMultiSetBoundaryCondition(Solver* solverPtr,
                                    BCFunctionPtr boundaryConditionFunctionPtr,
                                    PhysicalDomainT& domain,
                                    ObjectDataStructureBaseT& object,
                                    const std::string& fieldName, realT time );

template<typename BCFunctionPtr>
void ApplyBoundaryCondition(BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const std::string& fieldName,
                            realT time );

/////////////////////////////////////////////////


template< typename T >
inline void SetValue( const lSet& set,
                      ObjectDataStructureBaseT& object,
                      BoundaryConditionBase* bc,
                      array<T>& fieldData,
                      realT time )
{
  T replace;
  for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a )
  {
    T& fieldValue = fieldData[*a];

    replace = bc->GetDirection(time);
    replace *= Dot( bc->GetDirection(time), fieldValue );
    fieldValue -= replace;

    replace = bc->GetDirection(time);
    replace *= bc->GetValue(object,a,time);
    fieldValue += replace;
  }
}

template<>
inline void SetValue<realT>( const lSet& set,
                             ObjectDataStructureBaseT& object,
                             BoundaryConditionBase* bc,
                             array<realT>& fieldData,
                             realT time )
{
  for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a )
  {
    realT& fieldValue = fieldData[*a];
    fieldValue = bc->GetValue(object,a,time);
  }

}

template< typename T >
void ApplyDirichletBoundaryCondition( ObjectDataStructureBaseT& object,
                                      const std::string fieldName,
                                      realT time )
{

//    const TableManager& tableManager = TableManager::Instance();

  array<T>& fieldData = object.GetFieldData<T>(fieldName);

  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {
    // check to see if the requested field has a boundary condition applied to
    // it.
    BoundaryConditionBase* bc = *bcItr;
//      if( streq( bc->GetFieldName(time), fieldName) )
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && streq( bc->GetFieldName(time), fieldName) )
    {

      for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
      {
        int findSet = 1;

        std::map< std::string, lSet >::iterator setMap = object.m_Sets.find( bc->m_setNames[i] );
        if( setMap != object.m_Sets.end() )
        {
          lSet& set = setMap->second;

          if( bc->GetComponent(time) != -1 )
          {
            lSet::iterator aend = set.end();
            for( lSet::iterator a=set.begin() ; a!=aend ; ++a )
            {
              realT& fieldComponent = FieldComponent(fieldData[*a], bc->GetComponent(time) );
              fieldComponent = bc->GetValue(object,a,time);
            }
          }
          else
          {
            SetValue( set, object, bc, fieldData,time );
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
          msg << "Boundary condition for field \'" << fieldName  << "\' is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
          throw GPException(msg.str());
        }
      }
    }
  }
}



void ApplyRigidWallBoundaryCondition( ObjectDataStructureBaseT& object, const realT dt );

/*
   void ApplyTractionBoundaryCondition( FaceManagerT& faceManager,
                                       NodeManagerT& nodeManager,
                                       realT time);
 */

void ApplyTractionBoundaryCondition( PhysicalDomainT& domain,
                                     realT time);


//void ApplyTractionBoundaryCondition(PhysicalDomainT& domain,
// ObjectDataStructureBaseT& object, BoundaryConditionBase* bc_base, const lSet&
// set);

/***************
 *
 * Old BC methods
 *
 * NB these function calls use the fieldName to check if the boundary condition
 * should be applied
 * This is a legacy from the implementation in our early explicit solvers
 *
 * In general it is better (more extensible and maintainable) to check if the BC
 * is appropriate for the solver by checking if the
 * BC can be upcast to the desired BC
 *
 ****************
 */

/**
 * @author walsh24
 *
 * Apply a solver-defined boundary condition to a set of objects in a field
 *
 * @param solverPtr Pointer to the solver (usually "this")
 * @param boundaryConditionFunctionPtr  Pointer to the boundary-condition
 * function to be applied
 * @param domain The physical domain
 * @param object The object to apply the boundary condition to.
 * @param fieldName The name of the field to apply the boundary condition to.
 *
 * Use:
 *
 * Call used by the implicit ADRSolver to apply a dirichlet boundary condition
 * to the node field named in the concentrationFieldName_ string:
 * ApplyBoundaryCondition<realT>(this,&ADRSolver::DirichletBoundaryCondition,
 * domain, domain.m_nodeManager, concentrationFieldName_);
 *
 */
template< typename T, typename Solver, typename BCFunctionPtr>
void ApplyBoundaryCondition(Solver* solverPtr,
                            BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const std::string& fieldName,
                            const realT time )
{


  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {

    // check if field has the boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && streq( bc->GetFieldName(time), fieldName) )
    {

      for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
      {
        std::map< std::string, lSet >::iterator setMap = object.m_Sets.find( bc->m_setNames[i] );

        int findSet = 1;
        if( setMap != object.m_Sets.end() )
        {

          lSet& set = setMap->second;

          // std::cout << bc->m_fieldName << " " << bc->m_setName  << " "<<
          // bc->m_component << " " << value << std::endl;
          (solverPtr->*boundaryConditionFunctionPtr)(domain,object,bc,set,time);
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
          msg << "Boundary condition for field \'" << fieldName  << "\' is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
          throw GPException(msg.str());
        }
      }
    }
  }
}

// as above but with dt passed to the boundary condition function
template< typename T, typename Solver, typename BCFunctionPtr>
void ApplyBoundaryCondition(Solver* solverPtr,
                            BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const std::string& fieldName,
                            const realT time,
                            const realT dt )
{


  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {

    // check if field has the boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
//      if( streq( bc->GetFieldName(time), fieldName) ){
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && streq( bc->GetFieldName(time), fieldName) )
    {

      for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
      {
        int findSet = 1;
        std::map< std::string, lSet >::iterator setMap = object.m_Sets.find( bc->m_setNames[i] );
        if( setMap != object.m_Sets.end() )
        {

          lSet& set = setMap->second;

          // std::cout << bc->m_fieldName << " " << bc->m_setName  << " "<<
          // bc->m_component << " " << value << std::endl;
          (solverPtr->*boundaryConditionFunctionPtr)(domain,object,bc,set,time,dt);
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
          msg << "Boundary condition for field \'" << fieldName  << "\' is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
          throw GPException(msg.str());
        }

      }
    }
  }
}

// as above but with dt and dofOffset passed to the boundary condition function
template< typename T, typename Solver, typename BCFunctionPtr>
void ApplyBoundaryCondition(Solver* solverPtr,
                            BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const std::string& fieldName,
                            const realT time,
                            const realT dt,
                            const int dofOffset )
{


  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {

    // check if field has the boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
//      if( streq( bc->GetFieldName(time), fieldName) ){
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && streq( bc->GetFieldName(time), fieldName) )
    {

      for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
      {
        int findSet = 1;

        std::map< std::string, lSet >::iterator setMap = object.m_Sets.find( bc->m_setNames[i] );
        if( setMap != object.m_Sets.end() )
        {

          lSet& set = setMap->second;

          // std::cout << bc->m_fieldName << " " << bc->m_setName  << " "<<
          // bc->m_component << " " << value << std::endl;
          (solverPtr->*boundaryConditionFunctionPtr)(domain,object,bc,set,time,dt, dofOffset );
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
          msg << "Boundary condition for field \'" << fieldName  << "\' is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
          throw GPException(msg.str());
        }

      }
    }
  }
}


/**
 * @author walsh24
 *
 * Apply a solver boundary condition defined for multiple sets to a field.
 *
 * This function call should be employed if the boundary condition requires
 * interaction between two or more sets
 * simultaneously (eg. contact between two surfaces). Otherwise the
 * ApplyBoundaryCondition function is more appropriate.
 *
 * @param solverPtr Pointer to the solver (usually "this")
 * @param boundaryConditionFunctionPtr  Pointer to the boundary-condition
 * function to be applied
 * @param domain The physical domain
 * @param object The object to apply the boundary condition to.
 * @param fieldName The name of the field to apply the boundary condition to.
 *
 *
 */
template< typename T, typename Solver, typename BCFunctionPtr>
void ApplyMultiSetBoundaryCondition(Solver* solverPtr,
                                    BCFunctionPtr boundaryConditionFunctionPtr,
                                    PhysicalDomainT& domain,
                                    ObjectDataStructureBaseT& object,
                                    const std::string& fieldName, realT time ){

  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {

    // check if field has the boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && streq( bc->GetFieldName(time), fieldName) )
    {
//      if( streq( bc->GetFieldName(time), fieldName) ){
      (solverPtr->*boundaryConditionFunctionPtr)(domain,object, bc,time);
    }
  }
}


/**
 * @author walsh24
 *
 * Apply a boundary condition function to a field
 *
 * @param boundaryConditionFunctionPtr  Pointer to the boundary-condition
 * function to be applied
 * @param domain The physical domain
 * @param object The object to apply the boundary condition to.
 * @param fieldName The name of the field to apply the boundary condition to.
 *
 * Use this call when the boundary condition function is not part of the solver
 */
template<typename BCFunctionPtr>
void ApplyBoundaryCondition(BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const std::string& fieldName,
                            realT time ){


  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {

    // check if field has the boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && streq( bc->GetFieldName(time), fieldName) )
    {
//      if( streq( bc->GetFieldName(time), fieldName) ){

      for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
      {
        int findSet = 1;

        std::map< std::string, lSet >::iterator setMap = object.m_Sets.find( bc->m_setNames[i] );
        if( setMap != object.m_Sets.end() )
        {

          lSet& set = setMap->second;

          // std::cout << bc->m_fieldName << " " << bc->m_setName  << " "<<
          // bc->m_component << " " << value << std::endl;
          (*boundaryConditionFunctionPtr)(domain,object,&(*bc),set,time);

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
          msg << "Boundary condition for field \'" << fieldName  << "\' is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
          throw GPException(msg.str());
        }

      }
    }
  }
}

/**
 *  New BC methods
 *
 *  These BC's are distinguished only by the BCClass type
 *  The name of the field is not checked.
 *
 *  For example this allows Traction boundary conditions to be applied without
 * knowing if
 *  it is a "Traction", "UniformTraction", "RadialPressure" or some other
 * derived type
 *
 *
 */

/**
 * @author walsh24
 *
 * Apply a solver-defined boundary condition to a set of objects in a field
 *
 * @param solverPtr Pointer to the solver (usually "this")
 * @param boundaryConditionFunctionPtr  Pointer to the boundary-condition
 * function to be applied
 * @param domain The physical domain
 * @param object The object to apply the boundary condition to.
 * @param fieldName The name of the field to apply the boundary condition to.
 *
 * Use:
 *
 */
template< typename Solver, typename BCClass, typename BCFunctionPtr>
void ApplyBoundaryCondition(Solver* solverPtr,
                            BCFunctionPtr boundaryConditionFunctionPtr,
                            PhysicalDomainT& domain,
                            ObjectDataStructureBaseT& object,
                            const realT time )
{


  // iterate over all boundary conditions.
  for( array<BoundaryConditionBase*>::const_iterator bcItr=object.m_bcData.begin() ; bcItr!=object.m_bcData.end() ; ++bcItr )
  {

    // check if field has the boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;

    if( time >= bc->GetStartTime() && time < bc->GetEndTime()  )
    {

      BCClass* valid_bc = bc->UpcastActiveBCPointer<BCClass>(time); // use
                                                                    // Upcast
                                                                    // function
                                                                    // (rather
                                                                    // than
                                                                    // dynamic_cast<BCClass*>(
                                                                    // bc) ) to
                                                                    // enable
                                                                    // conditional
                                                                    // boundary
                                                                    // conditions

      if( valid_bc )
      {

        for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
        {
          int findSet = 1;

          std::map< std::string, lSet >::iterator setMap = object.m_Sets.find( bc->m_setNames[i] );
          if( setMap != object.m_Sets.end() )
          {

            lSet& set = setMap->second;

            (solverPtr->*boundaryConditionFunctionPtr)(domain,object,valid_bc,set,time);
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
            msg << "Some boundary condition is applied to set \'" << bc->m_setNames[i] << "\', which is undefined.";
            throw GPException(msg.str());
          }

        }
      }
    }
  }
}



void BuildKinematicConstraintBoundaryCondition( NodeManager& nodeManager,
                                                FaceManagerT& faceManager,
                                                array<lSet>& KinematicConstraintNodes,
                                                const realT tiedNodeTolerance );


void ApplyKinematicConstraintBoundaryCondition( FaceManagerT& faceManager,
                                                NodeManager& nodeManager,
                                                array<lSet>& KinematicConstraintNodes,
                                                const realT tiedNodeNormalRuptureStress,
                                                const realT tiedNodeShearRuptureStress  );

}
#endif /* APPLYBOUNDARYCONDITIONS_H_ */
