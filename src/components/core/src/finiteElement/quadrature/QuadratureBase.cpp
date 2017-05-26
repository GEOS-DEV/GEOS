#include "QuadratureBase.hpp"

namespace geosx
{

QuadratureBase::CatalogInterface::CatalogType& QuadratureBase::GetCatalog()
{
  static QuadratureBase::CatalogInterface::CatalogType catalog;
  return catalog;
}
}
