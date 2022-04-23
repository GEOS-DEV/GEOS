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
 * @file CoupledSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_

#include "constitutive/solid/CoupledSolidBase.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE type of solid model
 * @tparam PORO_TYPE type of porosity model
 * @tparam PERM_TYPE type of permeability model
 */
template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
class CoupledSolidUpdates
{
public:
  /**
   * @brief Constructor
   */
  CoupledSolidUpdates( SOLID_TYPE const & solidModel,
                       PORO_TYPE const & porosityModel,
                       PERM_TYPE const & permModel ):
    m_solidUpdate( solidModel.createKernelUpdates() ),
    m_porosityUpdate( porosityModel.createKernelUpdates() ),
    m_permUpdate( permModel.createKernelWrapper() )
  {}

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_porosityUpdate.numGauss(); }

  GEOSX_HOST_DEVICE
  real64 getOldPorosity( localIndex const k,
                         localIndex const q ) const
  {
    return m_porosityUpdate.getOldPorosity( k, q );
  }

  GEOSX_HOST_DEVICE
  real64 getPorosity( localIndex const k,
                      localIndex const q ) const
  {
    return m_porosityUpdate.getPorosity( k, q );
  }

  GEOSX_HOST_DEVICE
  real64 getInitialPorosity( localIndex const k,
                             localIndex const q ) const
  {
    return m_porosityUpdate.getInitialPorosity( k, q );
  }

  GEOSX_HOST_DEVICE
  virtual void updateStateFromPressure( localIndex const k,
                                        localIndex const q,
                                        real64 const & pressure ) const
  {
    GEOSX_UNUSED_VAR( k, q, pressure );
  }

protected:
  typename SOLID_TYPE::KernelWrapper const m_solidUpdate;

  typename PORO_TYPE::KernelWrapper const m_porosityUpdate;

  typename PERM_TYPE::KernelWrapper const m_permUpdate;
};



/**
 * @brief Class to represent a material which is formed by coupling several constitutive models.
 * It is used as an interface to access all constitutive models relative to the material properties.
 *
 * @tparam SOLID_TYPE type of solid model
 * @tparam PORO_TYPE type of porosity model
 * @tparam PERM_TYPE type of permeability model
 */
//START_SPHINX_INCLUDE_00
template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
class CoupledSolid : public CoupledSolidBase
//END_SPHINX_INCLUDE_00
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CoupledSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CoupledSolid() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
   */
  CoupledSolidUpdates< SOLID_TYPE, PORO_TYPE, PERM_TYPE > createKernelUpdates() const
  {

    return CoupledSolidUpdates< SOLID_TYPE, PORO_TYPE, PERM_TYPE >( getSolidModel(),
                                                                    getPorosityModel(),
                                                                    getPermModel() );
  }

  //START_SPHINX_INCLUDE_01
protected:
  SOLID_TYPE const & getSolidModel() const
  { return this->getParent().template getGroup< SOLID_TYPE >( m_solidModelName ); }

  PORO_TYPE const & getPorosityModel() const
  { return this->getParent().template getGroup< PORO_TYPE >( m_porosityModelName ); }

  PERM_TYPE const & getPermModel() const
  { return this->getParent().template getGroup< PERM_TYPE >( m_permeabilityModelName ); }
  //END_SPHINX_INCLUDE_01

};


template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE, PERM_TYPE >::CoupledSolid( string const & name, Group * const parent ):
  CoupledSolidBase( name, parent )
{}

template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE, PERM_TYPE >::~CoupledSolid() = default;


template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
void CoupledSolid< SOLID_TYPE, PORO_TYPE, PERM_TYPE >::initializePreSubGroups()
{
  if( PORO_TYPE::catalogName() != getPorosityModel().getCatalogName() )
  {
    GEOSX_ERROR( " The coupled solid "<< this->getName()<<
                 " expects a porosity model of type "<<PORO_TYPE::catalogName()<<
                 " but the specified porosity model \""<<m_porosityModelName<<
                 "\" is of type " << getPorosityModel().getCatalogName() );
  }
  if( PERM_TYPE::catalogName() != getPermModel().getCatalogName() )
  {
    GEOSX_ERROR( " The coupled solid "<<this->getName()<<
                 " expects a permeability model of type "<<PERM_TYPE::catalogName()<<
                 " but the specified permeability model \""<<m_permeabilityModelName<<
                 "\" is of type " << getPermModel().getCatalogName() );
  }
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

#endif /* GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_ */
