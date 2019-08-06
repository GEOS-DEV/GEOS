/*
 * ContactRelationBase.hpp
 *
 *  Created on: Aug 3, 2019
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_

#include "../ConstitutiveBase.hpp"
#include "managers/Functions/FunctionBase.hpp"

namespace geosx
{


namespace constitutive
{


class ContactRelationBase : public constitutive::ConstitutiveBase
{
public:
  ContactRelationBase( string const & name,
                       ManagedGroup * const parent );
  virtual ~ContactRelationBase() override;

  static string CatalogName() { return "Contact"; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override {}

  virtual ManagedGroup * CreateChild( string const & catalogKey,
                                      string const & name ) override;

  virtual void InitializePreSubGroups( ManagedGroup * const ) override;


  inline real64 stiffness() const { return m_penaltyStiffness; }

  inline real64 effectiveAperture( real64 const aperture ) const { return m_apertureFunction->Evaluate( & aperture ); }
  inline real64 dEffectiveAperture_dAperture( real64 const aperture ) const
  {
    real64 aperPlus = aperture ;
    real64 aperMinus = aperture - 1.0e-6;
    real64 slope = (m_apertureFunction->Evaluate( &aperPlus ) - m_apertureFunction->Evaluate( &aperMinus ) ) / 1.0e-6;
    return slope;
  }

  inline void apertureForPermeablityCalculation( real64 const aper0,
                                                 real64 const aper,
                                                 integer const integrationOption,
                                                 real64 & aperTerm,
                                                 real64 & dAperTerm_dAper )
  {
    // forward euler
    if( integrationOption == 0 )
    {
      aperTerm = aper0*aper0*aper0 ;
      dAperTerm_dAper = 0.0;
    }
    // simpsons rule / exact
    else if ( integrationOption == 1 )
    {
      aperTerm = 0.25 * ( aper0*aper0*aper0 +
                          aper0*aper0*aper +
                          aper0*aper*aper +
                          aper*aper*aper );
      dAperTerm_dAper = 0.25 * ( aper0*aper0 +
                                 2*aper0*aper +
                                 3*aper*aper );
    }
    else
    {
      GEOS_ERROR("invalid integrationOption");
    }

    return;
  }

  struct viewKeyStruct: public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto penaltyStiffnessString  = "penaltyStiffness";
  };

private:
  real64 m_penaltyStiffness;
  FunctionBase * m_apertureFunction;

};

}
} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_ */
