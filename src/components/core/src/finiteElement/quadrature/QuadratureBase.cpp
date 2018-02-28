#include "QuadratureBase.hpp"

namespace geosx
{
QuadratureBase::~QuadratureBase()
{}


QuadratureBase::CatalogInterface::CatalogType& QuadratureBase::GetCatalog()
{
  static QuadratureBase::CatalogInterface::CatalogType catalog;
  return catalog;
}
}
