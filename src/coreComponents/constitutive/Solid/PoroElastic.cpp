/*
 * PoroElastic.cpp
 *
 *  Created on: Mar 28, 2019
 *      Author: rrsettgast
 */

#include "PoroElastic.hpp"

#include "LinearElasticAnisotropic.hpp"
#include "LinearElasticIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
PoroElastic<BASE>::PoroElastic( string const & name, ManagedGroup * const parent ):
  BASE( name, parent ),
  m_compressibility(),
  m_referencePressure(),
  m_biotCoefficient(),
  m_poreVolumeMultiplier(),
  m_dPVMult_dPressure(),
  m_poreVolumeRelation()
{
  this->RegisterViewWrapper( viewKeyStruct::biotCoefficientString, &m_biotCoefficient, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Biot's coefficient");

  this->RegisterViewWrapper( viewKeyStruct::compressibilityString, &m_compressibility, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid Compressibilty");

  this->RegisterViewWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("ReferencePressure");


  this->RegisterViewWrapper( viewKeyStruct::poreVolumeMultiplierString, &m_poreVolumeMultiplier, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("");

  this->RegisterViewWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("");
}

template< typename BASE >
PoroElastic<BASE>::~PoroElastic()
{
}

template< typename BASE >
void PoroElastic<BASE>::PostProcessInput()
{
  //    m_compressibility = 1 / K;

  if (m_compressibility <= 0)
  {
//    string const message = std::to_string( numConstantsSpecified ) + " Elastic Constants Specified. Must specify 2 constants!";
//    GEOS_ERROR( message );
  }
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );

}

template< typename BASE >
void PoroElastic<BASE>::DeliverClone( string const & name,
                                      ManagedGroup * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<PoroElastic<BASE> >( name, parent );
  }
  BASE::DeliverClone( name, parent, clone );
  PoroElastic<BASE> * const newConstitutiveRelation = dynamic_cast<PoroElastic<BASE> *>(clone.get());

  newConstitutiveRelation->m_compressibility      = m_compressibility;
  newConstitutiveRelation->m_referencePressure    = m_referencePressure;
  newConstitutiveRelation->m_biotCoefficient      = m_biotCoefficient;
  newConstitutiveRelation->m_poreVolumeMultiplier = m_poreVolumeMultiplier;
  newConstitutiveRelation->m_dPVMult_dPressure    = m_dPVMult_dPressure;
  newConstitutiveRelation->m_poreVolumeRelation   = m_poreVolumeRelation;
}

template< typename BASE >
void PoroElastic<BASE>::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  BASE::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_poreVolumeMultiplier.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier = 1.0;

}

typedef PoroElastic<LinearElasticIsotropic> PoroLinearElasticIsotropic;
typedef PoroElastic<LinearElasticAnisotropic> PoroLinearElasticAnisotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroLinearElasticIsotropic, string const &, ManagedGroup * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroLinearElasticAnisotropic, string const &, ManagedGroup * const )

}
} /* namespace geosx */
