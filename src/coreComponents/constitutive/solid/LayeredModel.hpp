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
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_LAYEREDMODEL_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LAYEREDMODEL_HPP_

#include "constitutive/solid/SolidBase.hpp"
#include "ElasticIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"


namespace geosx
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
                       const bool & disableInelasticity ):
    SolidBaseUpdates( newStress, oldStress, disableInelasticity ),                       
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

//TODO: modify implementation of smallStrainUpdate to use optimized stiffness -
// this implementation uses full stiffness tensor
GEOSX_HOST_DEVICE
void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}

GEOSX_HOST_DEVICE
void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override
{


  //Get trial stress and elastic stiffness by disabling inelasticity
 // m_solidUpdate.m_disableInelasticity = true;
  //m_solidUpdateLayer1.smallStrainUpdate( k, q, timeIncrement, strainIncrement, trialStress, elasticStiffness );

  //Enable inelasticity to get the rate-independent update
  //m_solidUpdate.m_disableInelasticity = false;
  m_solidUpdateLayer2.smallStrainUpdate( k, q, strainIncrement, stress, stiffness );


  // for( localIndex i=0; i<6; ++i )
  // {
  //   stress[i] = 0.5 *  trialStress[i] + (1-timeRatio) * stress[i];
  //   for( localIndex j=0; j<6; ++j )
  //   {
  //     stiffness[i][j] = timeRatio * elasticStiffness[i][j]  + (1 - timeRatio) * stiffness[i][j];
  //   }
  // }

 // saveStress( k, q, stress );

  return;
}     
                            

protected:
  typename SOLID_TYPE1::KernelWrapper m_solidUpdateLayer1;
  typename SOLID_TYPE2::KernelWrapper m_solidUpdateLayer2;
 // real64 const m_relaxationTime;
};




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
  LayeredModelUpdates< SOLID_TYPE1 , SOLID_TYPE2 > createKernelUpdates() const
  {
    return LayeredModelUpdates< SOLID_TYPE1 , SOLID_TYPE2 >( getSolidModelLayer1(),
                                                  getSolidModelLayer2(),
                                                  m_newStress,
                                                  m_oldStress,
                                                  m_disableInelasticity );
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

  /// the name of the solid model
  string m_solidModelNameLayer1;

  string m_solidModelNameLayer2;

  //real64 m_relaxationTime;
  
  SOLID_TYPE1 const & getSolidModelLayer1() const
  { return this->getParent().template getGroup< SOLID_TYPE1 >( m_solidModelNameLayer1 ); }

  SOLID_TYPE2 const & getSolidModelLayer2() const
  { return this->getParent().template getGroup< SOLID_TYPE2 >( m_solidModelNameLayer2 ); }

};

// template< typename SOLID_TYPE >
// DuvautLionsSolid< SOLID_TYPE >::DuvautLionsSolid( string const & name, Group * const parent ):
//   SolidBase( name, parent )
// {}

// template< typename SOLID_TYPE >
// DuvautLionsSolid< SOLID_TYPE >::~DuvautLionsSolid() = default;



// template< typename SOLID_TYPE >
// void DuvautLionsSolid< SOLID_TYPE >::initializePreSubGroups()
// {
//   if( SOLID_TYPE::catalogName() != getSolidModel().getCatalogName() )
//   {
//     GEOSX_ERROR( " The coupled solid "<<this->getName()<<
//                  " expects a solid model of type "<<SOLID_TYPE::catalogName()<<
//                  " but the specified solid model \""<<this->m_solidModelName<<
//                  "\" is of type" << getSolidModel().getCatalogName() );
//   }
// }

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LAYEREDMODEL_HPP_ */
