/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

template< typename SOLID_TYPE,
          typename PORO_TYPE >
class CoupledSolidUpdates
{
public:
  /**
   * @brief Constructor
   */
  CoupledSolidUpdates( SOLID_TYPE * solidModel,
                       PORO_TYPE * porosityModel ):
    m_solidUpdate( solidModel->createKernelUpdates() ),
    m_porosityUpdate( porosityModel->createKernelUpdates() )
  {}

  /// Deleted default constructor
  CoupledSolidUpdates() = delete;

  /// Default copy constructor
  CoupledSolidUpdates( CoupledSolidUpdates const & ) = default;

  /// Default move constructor
  CoupledSolidUpdates( CoupledSolidUpdates && ) = default;

  /// Deleted copy assignment operator
  CoupledSolidUpdates & operator=( CoupledSolidUpdates const & ) = delete;

  /// Deleted move assignment operator
  CoupledSolidUpdates & operator=( CoupledSolidUpdates && ) =  delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  typename SOLID_TYPE::KernelWrapper const getSolidModel() const
  {
    return m_solidUpdate;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  typename PORO_TYPE::KernelWrapper const getPorosityModel() const
  {
    return m_porosityUpdate;
  }

protected:
  typename SOLID_TYPE::KernelWrapper const m_solidUpdate;

  typename PORO_TYPE::KernelWrapper const m_porosityUpdate;
};



/**
 * @brief Class to represent a coupled solid model
 */
template< typename SOLID_TYPE,
          typename PORO_TYPE >
class CoupledSolid : public ConstitutiveBase
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

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Poro" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// Post-process XML input
  virtual void postProcessInput() override;

  struct viewKeyStruct
  {
    static constexpr char const * solidModelNameString() { return "solidModelName"; }
    static constexpr char const * porosityModelNameString() { return "porosityModelName"; }
  };

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
   */
  CoupledSolidUpdates< SOLID_TYPE, PORO_TYPE > createKernelUpdates() const
  {

    return CoupledSolidUpdates< SOLID_TYPE, PORO_TYPE >( m_solidModel,
                                                         m_porosityModel );
  }

protected:

  // the solid model
  SOLID_TYPE * m_solidModel;

  // the porosity model
  PORO_TYPE * m_porosityModel;

  // the name of the solid model
  string m_solidModelName;

  // the name of the porosity model
  string m_porosityModelName;


  // PERMEABILITY_TYPE * m_permModel;
};


template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::CoupledSolid( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_solidModel( nullptr ),
  m_porosityModel( nullptr ),
  m_solidModelName(),
  m_porosityModelName()
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model." );

  registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the porosity model." );
}

template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::~CoupledSolid()
{}

template< typename SOLID_TYPE,
          typename PORO_TYPE >
void CoupledSolid< SOLID_TYPE, PORO_TYPE >::postProcessInput()
{
  m_solidModel = &this->getParent().template getGroup< SOLID_TYPE >( m_solidModelName );
  m_porosityModel = &this->getParent().template getGroup< PORO_TYPE >( m_porosityModelName );
}


}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_ */
