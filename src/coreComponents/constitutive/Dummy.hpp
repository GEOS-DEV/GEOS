/*
 * Dummy.hpp
 *
 *  Created on: Apr 14, 2020
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_DUMMY_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_DUMMY_HPP_

#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{

class Dummy : public constitutive::ConstitutiveBase
{
public:

  Dummy( string const & name,
         Group * const parent );

  virtual ~Dummy();


  virtual void DeliverClone( string const & GEOSX_UNUSED_PARAM(name),
                             Group * const GEOSX_UNUSED_PARAM(parent),
                             std::unique_ptr< ConstitutiveBase > & GEOSX_UNUSED_PARAM(clone) ) const override final;

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Dummy";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string CatalogName() { return m_catalogNameString; }

  virtual string GetCatalogName() override { return CatalogName(); }


  /// @typedef Alias for LinearElasticIsotropicUpdates
  using KernelWrapper = double;

  double createKernelWrapper( bool const GEOSX_UNUSED_PARAM(includeState) = false )
  {
    return 0.0;
  }
};

} // constitutive
} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_DUMMY_HPP_ */
