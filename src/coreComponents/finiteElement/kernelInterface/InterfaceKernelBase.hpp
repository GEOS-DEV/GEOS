/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ImplicitKernelBase.hpp"

/**
 * @file InterfaceKernelBase.hpp
 */

#ifndef GEOS_FINITEELEMENT_INTERFACEKERNELBASE_HPP_
#define GEOS_FINITEELEMENT_INTERFACEKERNELBASE_HPP_

#ifndef SELECTED_FE_TYPES_2D
#define SELECTED_FE_TYPES_2D BASE_FE_TYPES_2D 
#endif

namespace geos
{

namespace finiteElement
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class InterfaceKernelBase : public ImplicitKernelBase< SUBREGION_TYPE,
                                                       CONSTITUTIVE_TYPE,
                                                       FE_TYPE,
                                                       NUM_DOF_PER_TEST_SP,
                                                       NUM_DOF_PER_TRIAL_SP >
{
public:

  using Base = ImplicitKernelBase< SUBREGION_TYPE,
                                   CONSTITUTIVE_TYPE,
                                   FE_TYPE,
                                   NUM_DOF_PER_TEST_SP,
                                   NUM_DOF_PER_TRIAL_SP >;

  InterfaceKernelBase( NodeManager const & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       arrayView1d< globalIndex const > const inputDofNumber,
                       globalIndex const rankOffset,
                       CRSMatrixView< real64, globalIndex const > const inputMatrix,
                       arrayView1d< real64 > const inputRhs,
                       real64 const inputDt ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt )
{}

  struct StackVariables  
  {
  //public:
  //
  //  GEOS_HOST_DEVICE
  //  StackVariables():
  //  {}
  };

  protected:

};


//using NestedMapType = std::map< string, array1d< localIndex > > const;
//using NestedMapFE = std::map< string, std::unique_ptr< geos::finiteElement::FiniteElementBase > >; 

template< typename POLICY,
          typename CONSTITUTIVE_BASE,
          typename KERNEL_FACTORY >
static
real64 interfaceBasedKernelApplication( MeshLevel & mesh,
                                        string const & targetRegionName,
                                        arrayView1d< localIndex const > const & faceElementsList, 
                                        FiniteElementBase const & subRegionFE,
                                        //NestedMapType & faceTypesToFaceElements,
                                        //NestedMapFE  & faceTypeToFiniteElements,
                                        string const & constitutiveStringName,
                                        KERNEL_FACTORY & kernelFactory )
{
  GEOS_MARK_FUNCTION;


  // save the maximum residual contribution for scaling residuals for convergence criteria.
  real64 maxResidualContribution = 0;
  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementManager = mesh.getElemManager();

  SurfaceElementRegion & region = elementManager.getRegion< SurfaceElementRegion >( targetRegionName );
  FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
  localIndex const targetRegionIndex = 0;

  constitutive::ConstitutiveBase* constitutiveRelation = nullptr;
  constitutive::NullModel* nullConstitutiveModel = nullptr;
  if( subRegion.template hasWrapper< string >( constitutiveStringName ) )
  {
    string const & constitutiveName = subRegion.template getReference< string >( constitutiveStringName );
    constitutiveRelation = &subRegion.template getConstitutiveModel( constitutiveName );
    std::cout << constitutiveName << std::endl;
  }
  else
  {
    nullConstitutiveModel = &subRegion.template registerGroup< constitutive::NullModel >( "nullModelGroup" );
    constitutiveRelation = nullConstitutiveModel;
  }
  std::cout << "constitutiveSringName: " << constitutiveStringName << std::endl;

  //for (const auto& [finiteElementName, faceElementsList] : faceTypesToFaceElements)
  //{
 
    localIndex const numElems = faceElementsList.size();
 
    std::cout << "numElems:" << numElems << std::endl;
    //std::cout << "FE: " << finiteElementName << std::endl;

    // Call the constitutive dispatch which converts the type of constitutive model into a compile time constant.
    constitutive::ConstitutivePassThru< CONSTITUTIVE_BASE >::execute( *constitutiveRelation,
                                                                    [&maxResidualContribution,
                                                                     &nodeManager,
                                                                     &edgeManager,
                                                                     &faceManager,
                                                                     targetRegionIndex,
                                                                     &kernelFactory,
                                                                     &subRegion,
                                                                     &subRegionFE,
                                                                     numElems]
                                                                      ( auto & castedConstitutiveRelation )
    {
 
      //FiniteElementBase & subRegionFE = *(faceTypeToFiniteElements[finiteElementName]);
       
      finiteElement::FiniteElementDispatchHandler< SELECTED_FE_TYPES_2D >::dispatch2D( subRegionFE,
                                                                                      [&maxResidualContribution,
                                                                                       &nodeManager,
                                                                                       &edgeManager,
                                                                                       &faceManager,
                                                                                       targetRegionIndex,
                                                                                       &kernelFactory,
                                                                                       &subRegion,
                                                                                       numElems,
                                                                                       &castedConstitutiveRelation] ( auto const finiteElement )
      //GEOS_UNUSED_VAR(subRegionFE);

      //finiteElement::FiniteElementDispatchHandler< SELECTED_FE_TYPES_2D >::dispatch2D( subRegionFE,
      //                                                                                [&maxResidualContribution,
      //                                                                                 numElems] ( auto const finiteElement )
      {
        //GEOS_UNUSED_VAR(finiteElement);
        auto kernel = kernelFactory.createKernel( nodeManager,
                                                  edgeManager,
                                                  faceManager,
                                                  targetRegionIndex,
                                                  subRegion,
                                                  finiteElement,
                                                  castedConstitutiveRelation);

        //GEOS_UNUSED_VAR( kernel);
        using KERNEL_TYPE = decltype( kernel );

        maxResidualContribution =
          std::max( maxResidualContribution,
                    KERNEL_TYPE::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernel ) );

      } );
    } );
  //}

  // Remove the null constitutive model (not required, but cleaner)
  if( nullConstitutiveModel )
  {
    subRegion.deregisterGroup( "nullModelGroup" );
  }

  //GEOS_UNUSED_VAR(kernelFactory, nodeManager, edgeManager, faceManager, constitutiveRelation, targetRegionIndex);

  return maxResidualContribution;
}
//END_interfaceBasedKernelApplication

} // namespace finiteElement
} // namespace geos

#endif /* GEOS_FINITEELEMENT_INTERFACEKERNELBASE_HPP_ */