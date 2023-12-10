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


/**
 * @file LayeredModel.hpp
 * @brief This class implements two-layer homogenization model
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_LAYEREDMODEL_HPP_
#define GEOS_CONSTITUTIVE_SOLID_LAYEREDMODEL_HPP_

#include "constitutive/solid/SolidBase.hpp"
#include "ElasticIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE1 type of solid model for layer 1
 * @tparam SOLID_TYPE2 type of solid model for layer 2
 */

template< typename SOLID_TYPE1 ,  
          typename SOLID_TYPE2 >
class LayeredModelUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   */
  LayeredModelUpdates( SOLID_TYPE1 const & solidModelLayer1,
                       SOLID_TYPE2 const & solidModelLayer2,
                       arrayView3d< real64, solid::STRESS_USD > const & newStress,
                       arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                       arrayView1d< real64 const > const & thermalExpansionCoefficient,
                       const bool & disableInelasticity ):
    SolidBaseUpdates( newStress, oldStress, thermalExpansionCoefficient, disableInelasticity ),                       
    m_solidUpdateLayer1( solidModelLayer1.createKernelUpdates() ),
    m_solidUpdateLayer2( solidModelLayer2.createKernelUpdates() )
  {}

  /// Deleted default constructor
  LayeredModelUpdates() = delete;

  /// Default copy constructor
  LayeredModelUpdates( LayeredModelUpdates const & ) = default;

  /// Default move constructor
  LayeredModelUpdates( LayeredModelUpdates && ) = default;

  /// Deleted copy assignment operator
  LayeredModelUpdates & operator=( LayeredModelUpdates const & ) = delete;

  /// Deleted move assignment operator
  LayeredModelUpdates & operator=( LayeredModelUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  /// Use base version of saveConvergedState
  using SolidBaseUpdates::saveConvergedState;
  


  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( & stress )[6],
                          real64 ( & stiffness )[6][6] ) const
  {
//Get trial stress and elastic stiffness by disabling inelasticity
 // m_solidUpdate.m_disableInelasticity = true;
  //m_solidUpdateLayer1.smallStrainUpdate( k, q, timeIncrement, strainIncrement, trialStress, elasticStiffness );

  //Enable inelasticity to get the rate-independent update
  //m_solidUpdate.m_disableInelasticity = false;
  std::cout<<"m_newStress size in LayeredModel = "<<m_newStress.size(0)<<std::endl;
  m_solidUpdateLayer2.smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );


  saveStress( k, q, stress );

  return;
  }

//TODO: modify implementation of smallStrainUpdate to use optimized stiffness -
// this implementation uses full stiffness tensor
//  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const
  {
    this->smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
  }

  protected:
  typename SOLID_TYPE1::KernelWrapper const m_solidUpdateLayer1;
  typename SOLID_TYPE2::KernelWrapper const m_solidUpdateLayer2;

};

/**
 * @brief LayeredModelBase class used for dispatch of all layered solids.
 */
class LayeredModelBase
{};

/**
 * @brief Class to represent a homogenization model for two-layer unit cell each made of an elastic or plastic model.
 * It is used as an interface to access all constitutive models relative to the material properties.
 *
 * @tparam SOLID_TYPE1 type of solid model for layer 1
 * @tparam SOLID_TYPE2 type of solid model for layer 2
 */

template< typename SOLID_TYPE1 ,
          typename SOLID_TYPE2 >
class LayeredModel : public SolidBase
{
public:

/// Alias for class LayeredModel : public SolidBase

using KernelWrapper = LayeredModelUpdates< SOLID_TYPE1, SOLID_TYPE2 >;
  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  LayeredModel( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~LayeredModel() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Layered" ) + SOLID_TYPE1::m_catalogNameString + SOLID_TYPE2::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

 // real64 relaxationTime() const { return m_relaxationTime; }

  string solidModelNameLayer1() { return m_solidModelNameLayer1; }

  string solidModelNameLayer2() { return m_solidModelNameLayer2; }

  /**
   * @brief Create a instantiation of the LayeredModelUpdates class
   *        that refers to the data in this.
   * @return An instantiation of LayeredModelSolidUpdates.
   */
   
  KernelWrapper createKernelUpdates() const
  {
    return LayeredModelUpdates< SOLID_TYPE1 , SOLID_TYPE2 >( getSolidModelLayer1(),
                                                  getSolidModelLayer2(),
                                                  m_newStress,
                                                  m_oldStress,
                                                  m_thermalExpansionCoefficient,
                                                  m_disableInelasticity );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams ) const
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          getSolidModelLayer1(),
                          getSolidModelLayer2(),
                          m_newStress,
                          m_oldStress,
                          m_thermalExpansionCoefficient,
                          m_disableInelasticity  );
  }

  //virtual void saveConvergedState() const override;

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for solid model name
    static constexpr char const * solidModelNameLayer1String() { return "solidModelNameLayer1"; }
    static constexpr char const * solidModelNameLayer2String() { return "solidModelNameLayer2"; }

  };

  std::vector< string > getSubRelationNames() const override final
  {
   std::vector< string > subRelationNames = { m_solidModelNameLayer1 , 
                                              m_solidModelNameLayer2 };
   return subRelationNames;
  }
  // KernelWrapper createKernelUpdates() const
  // {
  //   return DuvautLionsSolidUpdates< SOLID_TYPE >( getSolidModel(),
  //                         m_relaxationTime,
  //                         m_newStress,
  //                         m_oldStress,
  //                         m_disableInelasticity );
  //  }



protected:
  string m_solidModelNameLayer1;

  string m_solidModelNameLayer2;

  SOLID_TYPE1 & getSolidModelLayer1()
  { return this->getParent().template getGroup< SOLID_TYPE1 >( m_solidModelNameLayer1 ); }

  SOLID_TYPE2 & getSolidModelLayer2() 
  { return this->getParent().template getGroup< SOLID_TYPE2 >( m_solidModelNameLayer2 ); }

  SOLID_TYPE1 const & getSolidModelLayer1() const
  { return this->getParent().template getGroup< SOLID_TYPE1 >( m_solidModelNameLayer1 ); }

  SOLID_TYPE2 const & getSolidModelLayer2() const
  { return this->getParent().template getGroup< SOLID_TYPE2 >( m_solidModelNameLayer2 ); }
};

}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_LAYEREDMODEL_HPP_ */
