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

/*
 * BoundaryConditionManager.hpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#include "common/DataTypes.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "BoundaryConditionBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const boundaryConditionMananger("BoundaryConditionMananger");
}
}


class BoundaryConditionManager : public dataRepository::ManagedGroup
{
public:
  BoundaryConditionManager( string const & name, dataRepository::ManagedGroup * const parent );
  virtual ~BoundaryConditionManager() override;

  static BoundaryConditionManager * get();

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void ApplyInitialConditions( dataRepository::ManagedGroup * domain ) const;

  void ApplyBoundaryCondition( dataRepository::ManagedGroup * object,
                               std::string const & fieldName,
                               real64 const time );

//  template< typename ... ARGS>
//  void ApplyBoundaryCondition( dataRepository::ManagedGroup & object,
//                               std::string const & fieldName,
//                               real64 const time,
//                               ARGS & ... args );

  template< typename BCFunctionPtr, typename... ARGS>
  void ApplyBoundaryCondition( BCFunctionPtr boundaryConditionFunctionPtr,
                               dataRepository::ManagedGroup * object,
                               std::string const & fieldName,
                               real64 const time,
                               ARGS & ... args );

  template< typename Solver, typename BCFunctionPtr, typename... ARGS>
  void ApplyBoundaryCondition( Solver* solverPtr,
                               BCFunctionPtr boundaryConditionFunctionPtr,
                               dataRepository::ManagedGroup * object,
                               std::string const & fieldName,
                               real64 const time,
                               ARGS & ... args );


  template< typename LAMBDA >
  void ApplyBoundaryCondition( real64 const time,
                               dataRepository::ManagedGroup const * const object,
                               string const & fieldName,
                               LAMBDA && lambda )
  {
    dataRepository::ManagedGroup const * sets = object->GetGroup(dataRepository::keys::sets);

    // iterate over all boundary conditions.
    forSubGroups<BoundaryConditionBase>([&](BoundaryConditionBase * bc) -> void
    {
      if( time >= bc->GetStartTime() && time < bc->GetEndTime() && ( bc->GetFieldName()==fieldName) )
      {
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = sets->getWrapper<set<localIndex>>(setName);
          if( setWrapper != nullptr )
          {
            set<localIndex> const & set = setWrapper->reference();
            lambda( bc, set );
          }
        }
      }
    });
  }

  template< typename LAMBDA >
  void ApplyBoundaryCondition( real64 const time,
                               string const & fieldName,
                               LAMBDA && lambda )
  {
    // iterate over all boundary conditions.
    forSubGroups<BoundaryConditionBase>([&](BoundaryConditionBase * bc) -> void
    {
      if( time >= bc->GetStartTime() && time < bc->GetEndTime() && ( bc->GetFieldName()==fieldName) )
      {
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          lambda( bc, setName );
        }
      }
    });
  }
};



//
//template< typename ... ARGS>
//void BoundaryConditionManager::ApplyBoundaryCondition(
// dataRepository::ManagedGroup * object,
//                                                       std::string const &
// fieldName,
//                                                       real64 const time,
//                                                       ARGS & ... args )
//{
//  dataRepository::ManagedGroup const * sets =
// object->GetGroup(dataRepository::keys::sets);
//
//  // iterate over all boundary conditions.
//  forSubGroups<BoundaryConditionBase>([&](BoundaryConditionBase * bc) -> void
//  {
//    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && (
// bc->GetFieldName()==fieldName) )
//    {
//      string_array setNames = bc->GetSetNames();
//      for( auto & setName : setNames )
//      {
//        dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper =
// sets->getWrapperPtr<set<localIndex>>(setName);
//        if( setWrapper != nullptr )
//        {
//          set<localIndex> const & set = setWrapper->reference();
//          bc->ApplyBounaryConditionDefaultMethod(set,time,args...);
//        }
//      }
//    }
//  });
//}
//

template< typename BCFunctionPtr, typename... ARGS>
void BoundaryConditionManager::ApplyBoundaryCondition( BCFunctionPtr boundaryConditionFunctionPtr,
                                                       dataRepository::ManagedGroup * object,
                                                       std::string const & fieldName,
                                                       real64 const time,
                                                       ARGS & ... args )
{
  dataRepository::ManagedGroup const * sets = object->GetGroup(dataRepository::keys::sets);

  // iterate over all boundary conditions.
  forSubGroups<BoundaryConditionBase>([&](BoundaryConditionBase * bc) -> void
    {
      if( time >= bc->GetStartTime() && time < bc->GetEndTime() && ( bc->GetFieldName()==fieldName) )
      {
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = sets->getWrapper<set<localIndex>>(setName);
          if( setWrapper != nullptr )
          {
            set<localIndex> const & set = setWrapper->reference();
            (*boundaryConditionFunctionPtr)( bc, set, time, args... );
          }
        }
      }
    });
}



template< typename Solver, typename BCFunctionPtr, typename... ARGS>
void BoundaryConditionManager::ApplyBoundaryCondition( Solver* solverPtr,
                                                       BCFunctionPtr boundaryConditionFunctionPtr,
                                                       dataRepository::ManagedGroup * object,
                                                       std::string const & fieldName,
                                                       real64 const time,
                                                       ARGS & ... args )
{
  dataRepository::ManagedGroup const * sets = object->GetGroup(dataRepository::keys::sets);

  // iterate over all boundary conditions.
  forSubGroups<BoundaryConditionBase>([&](BoundaryConditionBase const * bc) -> void
    {
      if( time >= bc->GetStartTime() && time < bc->GetEndTime() && ( bc->GetFieldName()==fieldName) )
      {
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = sets->getWrapper<set<localIndex>>(setName);
          if( setWrapper != nullptr )
          {
            set<localIndex> const & set = setWrapper->reference();
            (solverPtr->*boundaryConditionFunctionPtr)(object, bc, set, time, args...);
          }
        }
      }
    });
}



} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
        */
