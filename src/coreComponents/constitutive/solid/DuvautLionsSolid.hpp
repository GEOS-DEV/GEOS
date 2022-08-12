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

#include "SolidBase.hpp"

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
class DuvautLionsSolidUpdates: public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   */
  DuvautLionsSolidUpdates( SOLID_TYPE const & solidModel):
    m_solidUpdate( solidModel.createKernelUpdates() )
  {}

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override;

protected:
  typename SOLID_TYPE::KernelWrapper const m_solidUpdate;
  real64 m_relaxationTime; 
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DuvautLionsSolidUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & timeIncrement,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 real64 ( & stiffness )[6][6] ) const
{
  real64 trialStress[6];   // Trial stress
  real64 elasticStiffness[6][6];  //Elastic stiffness
  real64 internalVariable;        //internal state variable
  real64 trialInternalVariable
  real64 timeRatio = timeIncrement / m_relaxationTime;

  for( localIndex i=0; i<6; ++i )
  {
    trialStress[i] = stress[i];
  }

  //Get trial stress and elastic stiffness by disabling inelasticity
  m_solidUpdate.m_disableInelasticity = true; 
  m_solidUpdate.smallStrainUpdate( k, q, timeIncrement, strainIncrement, trialStress, elasticStiffness );

  //enable inelasticity to get the rate-independent update
  m_solidUpdate.m_disableInelasticity = false;
  m_solidUpdate.smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );

  stress = (trialStress + timeRatio * stress) / (1 + timeRatio);

  for ( localIndex i=0; i<6; ++i )
  {
    for ( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] = (elasticStiffness[i][j] + timeRatio * stiffness[i][j]) / (1+timeRatio);
    }
  }


  //TO DO: add if statement to update the internal state variables depending on what SOLID_TYPE is.
  if( SOLID_TYPE::catalogName() == 'DruckerPrager' )
  {
   internalVariable =  m_solidUpdate.m_oldCohesion;
   trialInternalVariable = m_solidUpdate.m_newCohesion;
   internalVariable = (internalVariable + timeRatio * trialInternalVariable) / (1 + timeRatio);
   m_solidUpdate.m_newCohesion = internalVariable;

  }else
  {
    GEOSX_ERROR( " The Duvaut-Lions solid for "<<SOLID_TYPE::catalogName()<<" is not implemented." );
  }
}


/**
 * @brief Class to represent a rate-dependent Duvaut-Lions material coupled with a plasticity model.
 * It is used as an interface to access all constitutive models relative to the material properties.
 *
 * @tparam SOLID_TYPE type of solid model
 */
//START_SPHINX_INCLUDE_00
template< typename SOLID_TYPE>
class DuvautLionsSolid : public SolidBase
//END_SPHINX_INCLUDE_00
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
 DuvautLionsSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~DuvautLionsSolid() override;

  virtual void initializePreSubGroups() override;

  virtual void saveConvergedState() const override;

  /**
   * @brief Create a instantiation of the DuvautLionsSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of DuvautLionsSolidUpdates.
   */
  DuvautLionsSolidUpdates< SOLID_TYPE > createKernelUpdates() const
  {

    return DuvautLionsSolidUpdates< SOLID_TYPE >( getSolidModel());
  }

  //START_SPHINX_INCLUDE_01
protected:
  SOLID_TYPE const & getSolidModel() const
  { return this->getParent().template getGroup< SOLID_TYPE >( m_solidModelName ); }
  //END_SPHINX_INCLUDE_01

};


template< typename SOLID_TYPE>
DuvautLionsSolid< SOLID_TYPE >::DuvautLionsSolid( string const & name, Group * const parent ):
  SolidBase( name, parent )
{}

template< typename SOLID_TYPE >
DuvautLionsSolid< SOLID_TYPE >::~DuvautLionsSolid() = default;


template< typename SOLID_TYPE >
void DuvautLionsSolid< SOLID_TYPE >::initializePreSubGroups()
{
  if( SOLID_TYPE::catalogName() != getSolidModel().getCatalogName() )
  {
    GEOSX_ERROR( " The coupled solid "<<this->getName()<<
                 " expects a solid model of type "<<SOLID_TYPE::catalogName()<<
                 " but the specified solid model \""<<this->m_solidModelName<<
                 "\" is of type" << getSolidModel().getCatalogName() );
  }
}

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_ */
