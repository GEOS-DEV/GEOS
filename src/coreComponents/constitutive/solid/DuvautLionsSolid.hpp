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
 * @file DuvautLionsSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_

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
 * @tparam SOLID_TYPE type of solid model
 */
template< typename SOLID_TYPE >
class DuvautLionsSolidUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   */
  DuvautLionsSolidUpdates( SOLID_TYPE const & solidModel,
                           real64 const & relaxationTime,      
                           arrayView3d< real64, solid::STRESS_USD > const & newStress,
                           arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                           const bool & disableInelasticity ):
    SolidBaseUpdates( newStress, oldStress, disableInelasticity ),                       
    m_solidUpdate( solidModel.createKernelUpdates() ),
    m_relaxationTime( relaxationTime )
  {}

  /// Deleted default constructor
  DuvautLionsSolidUpdates() = delete;

  /// Default copy constructor
  DuvautLionsSolidUpdates( DuvautLionsSolidUpdates const & ) = default;

  /// Default move constructor
  DuvautLionsSolidUpdates( DuvautLionsSolidUpdates && ) = default;

  /// Deleted copy assignment operator
  DuvautLionsSolidUpdates & operator=( DuvautLionsSolidUpdates const & ) = delete;

  /// Deleted move assignment operator
  DuvautLionsSolidUpdates & operator=( DuvautLionsSolidUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

//TODO: modify implementation of smallStrainUpdate to use optimized stiffness -
// this implementation uses full stiffness tensor
GEOSX_HOST_DEVICE
void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
}

GEOSX_HOST_DEVICE
void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override
{
  real64 trialStress[6];   // Trial stress (elastic predictor)
  real64 elasticStiffness[6][6];  //Elastic stiffness
  real64 timeRatio = 1 / (1 + timeIncrement / m_relaxationTime);

  for( localIndex i=0; i<6; ++i )
  {
    trialStress[i] = stress[i];
  }

  //Get trial stress and elastic stiffness by disabling inelasticity
 // m_solidUpdate.m_disableInelasticity = true;
  m_solidUpdate.smallStrainUpdate( k, q, timeIncrement, strainIncrement, trialStress, elasticStiffness );

  //Enable inelasticity to get the rate-independent update
  //m_solidUpdate.m_disableInelasticity = false;
  m_solidUpdate.smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );


  for( localIndex i=0; i<6; ++i )
  {
    stress[i] = timeRatio *  trialStress[i] + (1-timeRatio) * stress[i];
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] = timeRatio * elasticStiffness[i][j]  + (1 - timeRatio) * stiffness[i][j];
    }
  }

  saveStress( k, q, stress );
  viscousStateUpdate( k, q, timeRatio );
  return;
}     
                            

protected:
  typename SOLID_TYPE::KernelWrapper m_solidUpdate;
  real64 const m_relaxationTime;
};




/**
 * @brief Class to represent a rate-dependent Duvaut-Lions material coupled with a plasticity model.
 * It is used as an interface to access all constitutive models relative to the material properties.
 *
 * @tparam SOLID_TYPE type of solid model
 */
//START_SPHINX_INCLUDE_00
template< typename SOLID_TYPE >
class DuvautLionsSolid : public SolidBase
//END_SPHINX_INCLUDE_00
{
public:

/// Alias for DuvautLionsSolidUpdates
using KernelWrapper = DuvautLionsSolidUpdates< SOLID_TYPE >;
  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  DuvautLionsSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~DuvautLionsSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Visco" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  //virtual void initializePreSubGroups() override;

  real64 relaxationTime() const { return m_relaxationTime; }

  string solidModelName() { return m_solidModelName; }

  /**
   * @brief Create a instantiation of the DuvautLionsSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of DuvautLionsSolidUpdates.
   */
  DuvautLionsSolidUpdates< SOLID_TYPE > createKernelUpdates() const
  {
    return DuvautLionsSolidUpdates< SOLID_TYPE >( getSolidModel(),
                                                  m_relaxationTime,
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
    static constexpr char const * solidModelNameString() { return "solidModelName"; }
    /// string/key for relaxation time
    static constexpr char const * relaxationTimeString() { return "relaxationTime"; }
  };

  std::vector< string > getSubRelationNames() const override final
  {
   std::vector< string > subRelationNames = { m_solidModelName };
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
  

  //START_SPHINX_INCLUDE_01
protected:

  /// the name of the solid model
  string m_solidModelName;

  real64 m_relaxationTime;
  
  SOLID_TYPE const & getSolidModel() const
  { return this->getParent().template getGroup< SOLID_TYPE >( m_solidModelName ); }
  //END_SPHINX_INCLUDE_01

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

#endif /* GEOSX_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_ */
