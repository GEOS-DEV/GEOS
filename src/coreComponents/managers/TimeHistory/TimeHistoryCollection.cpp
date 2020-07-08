#include "TimeHistoryCollection.hpp"
#include "managers/Tasks/TaskBase.hpp"

namespace geosx
{
REGISTER_CATALOG_ENTRY( TaskBase, HistoryCollection, std::string const &, Group * const );
}
