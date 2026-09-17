/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StressStrainAverageKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRESSSTRAINAVERAGEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRESSSTRAINAVERAGEKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "finiteElement/elementFormulations/FiniteElementOperators.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{


/**
 * @class AverageStressStrainOverQuadraturePoints
 * @tparam SUBREGION_TYPE the subRegion type
 * @tparam FE_TYPE the finite element type
 * @tparam SOLID_TYPE the solid mechanics constitutuve type
 */
template< typename FE_TYPE,
          typename SOLID_TYPE >
class AverageStressStrainOverQuadraturePoints :
  public AverageOverQuadraturePointsBase< CellElementSubRegion,
                                          FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = AverageOverQuadraturePointsBase< CellElementSubRegion,
                                                FE_TYPE >;

  using Base::m_elementVolume;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor for the class
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   * @param displacement the displacement solution field
   * @param avgStrain the strain averaged over quadrature points
   * @param stress the stress solution field
   * @param avgStress the stress averaged over quadrature points
   * @param temperature the temperature field (may be empty for non-thermal simulations)
   * @param temperature_n the temperature at the previous time step (may be empty for non-thermal simulations)
   */
  AverageStressStrainOverQuadraturePoints( NodeManager & nodeManager,
                                           EdgeManager const & edgeManager,
                                           FaceManager const & faceManager,
                                           CellElementSubRegion const & elementSubRegion,
                                           FE_TYPE const & finiteElementSpace,
                                           SOLID_TYPE const & solidModel,
                                           fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const displacement,
                                           fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const displacementInc,
                                           fields::solidMechanics::arrayView2dLayoutStrain const avgStrain,
                                           fields::solidMechanics::arrayView2dLayoutStrain const avgPlasticStrain,
                                           arrayView3d< real64 const, solid::STRESS_USD > const stress,
                                           fields::solidMechanics::arrayView2dLayoutAvgStress const avgStress,
                                           arrayView1d< real64 const > const temperature,
                                           arrayView1d< real64 const > const temperature_n ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace ),
    m_solidUpdate( solidModel.createKernelUpdates()),
    m_displacement( displacement ),
    m_displacementInc( displacementInc ),
    m_avgStrain( avgStrain ),
    m_avgPlasticStrain( avgPlasticStrain ),
    m_stress( stress ),
    m_avgStress( avgStress ),
    m_temperature( temperature ),
    m_temperature_n( temperature_n )
  {}

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : Base::StackVariables
  {real64 uLocal[FE_TYPE::maxSupportPoints][3];
   real64 uHatLocal[FE_TYPE::maxSupportPoints][3]; };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    Base::setup( k, stack );

    for( localIndex a = 0; a < FE_TYPE::maxSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );
      for( int i = 0; i < 3; ++i )
      {
        stack.uLocal[a][i] = m_displacement[localNodeIndex][i];
        stack.uHatLocal[a][i] = m_displacementInc[localNodeIndex][i];
      }
    }

    for( int icomp = 0; icomp < 6; ++icomp )
    {
      m_avgStrain[k][icomp] = 0.0;
      m_avgStress[k][icomp] = 0.0;
    }
  }

  /**
   * @brief Increment the average property with the contribution of the property at this quadrature point
   * @param k The element index
   * @param q The quadrature point index
   * @param stack The StackVariables object that hold the stack variables.
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    //real64 const weight = FE_TYPE::transformedQuadratureWeight( q, stack.xLocal, stack.feStack ) / m_elementVolume[k];

    real64 dNdX[ FE_TYPE::maxSupportPoints ][3];
    real64 const detJxW = FE_TYPE::calcGradN( q, stack.xLocal, stack.feStack, dNdX );
    real64 strain[6] = {0.0};
    real64 strainInc[6] = {0.0};
    finiteElement::feOps::symmetricGradient( dNdX, stack.uLocal, strain );
    finiteElement::feOps::symmetricGradient( dNdX, stack.uHatLocal, strainInc );

    real64 elasticStrainInc[6] = {0.0};
    m_solidUpdate.getElasticStrainInc( k, q, elasticStrainInc );

    real64 const thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );
    real64 const deltaTemperature = ( m_temperature.size() > 0 )
                                    ? ( m_temperature[k] - m_temperature_n[k] )
                                    : 0.0;

    real64 conversionFactor[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5}; // used for converting from engineering shear to tensor shear

    for( int icomp = 0; icomp < 6; ++icomp )
    {
      m_avgStrain[k][icomp] += conversionFactor[icomp]*detJxW*strain[icomp]/m_elementVolume[k];
      m_avgStress[k][icomp] += detJxW*m_stress[k][q][icomp]/m_elementVolume[k];

      // Thermal strain is purely volumetric: subtract thermal strain from normal components so that
      // only the mechanical (plastic) part is accumulated.
      real64 const thermalStrainInc = ( icomp < 3 ) ? thermalExpansionCoefficient * deltaTemperature : 0.0;
      real64 const mechanicalStrainInc = strainInc[icomp] - thermalStrainInc;

      // This is a hack to handle boundary conditions such as those seen in plane-strain wellbore problems
      // Essentially, if bcs are constraining the strain (and thus total displacement), we do not accumulate any plastic strain (regardless
      // of stresses in material law)
      if( std::abs( mechanicalStrainInc ) > 1.0e-8 )
      {
        m_avgPlasticStrain[k][icomp] += conversionFactor[icomp]*detJxW*(mechanicalStrainInc - elasticStrainInc[icomp])/m_elementVolume[k];
      }
    }
  }

  /**
   * @brief Launch the kernel over the elements in the subRegion
   * @tparam POLICY the kernel policy
   * @tparam KERNEL_TYPE the type of kernel
   * @param numElems the number of elements in the subRegion
   * @param kernelComponent the kernel component
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static void
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
    } );
  }

protected:

  /// The material
  typename SOLID_TYPE::KernelWrapper const m_solidUpdate;

  /// The displacement solution
  fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const m_displacement;

  /// The displacement increment
  fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const m_displacementInc;

  /// The average strain
  fields::solidMechanics::arrayView2dLayoutStrain const m_avgStrain;

  /// The average plastic strain
  fields::solidMechanics::arrayView2dLayoutStrain const m_avgPlasticStrain;

  /// The stress solution
  arrayView3d< real64 const, solid::STRESS_USD > const m_stress;

  /// The average stress
  fields::solidMechanics::arrayView2dLayoutAvgStress const m_avgStress;

  /// The temperature at the current time step (empty for non-thermal simulations)
  arrayView1d< real64 const > const m_temperature;

  /// The temperature at the previous time step (empty for non-thermal simulations)
  arrayView1d< real64 const > const m_temperature_n;

};



/**
 * @class AverageStressStrainOverQuadraturePointsKernelFactory
 * @brief Class to create and launch the kernel
 */
class AverageStressStrainOverQuadraturePointsKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam SUBREGION_TYPE the subRegion type
   * @tparam FE_TYPE the finite element type
   * @tparam SOLID_TYPE the constitutive type
   * @tparam POLICY the kernel policy
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   * @param property the property at quadrature points
   * @param averageProperty the property averaged over quadrature points
   */
  template< typename FE_TYPE,
            typename SOLID_TYPE,
            typename POLICY >
  static void
  createAndLaunch( NodeManager & nodeManager,
                   EdgeManager const & edgeManager,
                   FaceManager const & faceManager,
                   CellElementSubRegion const & elementSubRegion,
                   FE_TYPE const & finiteElementSpace,
                   SOLID_TYPE const & solidModel,
                   fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const displacement,
                   fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const displacementInc,
                   fields::solidMechanics::arrayView2dLayoutStrain const avgStrain,
                   fields::solidMechanics::arrayView2dLayoutStrain const avgPlasticStrain,
                   arrayView3d< real64 const, solid::STRESS_USD > const stress,
                   fields::solidMechanics::arrayView2dLayoutAvgStress const avgStress,
                   arrayView1d< real64 const > const temperature = {},
                   arrayView1d< real64 const > const temperature_n = {} )
  {
    AverageStressStrainOverQuadraturePoints< FE_TYPE, SOLID_TYPE >
    kernel( nodeManager, edgeManager, faceManager, elementSubRegion, finiteElementSpace,
            solidModel, displacement, displacementInc, avgStrain, avgPlasticStrain, stress, avgStress,
            temperature, temperature_n );

    AverageStressStrainOverQuadraturePoints< FE_TYPE, SOLID_TYPE >::template
    kernelLaunch< POLICY >( elementSubRegion.size(), kernel );
  }
};



}



#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRESSSTRAINAVERAGEKERNELS_HPP_ */
