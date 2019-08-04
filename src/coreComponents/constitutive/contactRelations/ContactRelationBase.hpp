/*
 * ContactRelationBase.hpp
 *
 *  Created on: Aug 3, 2019
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_

#include "../ConstitutiveBase.hpp"

namespace geosx
{

class FunctionBase;

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
