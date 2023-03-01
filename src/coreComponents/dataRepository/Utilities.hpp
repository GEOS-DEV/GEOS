
#include "common/DataTypes.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;

  /**
   */
/**
 * @brief Prints the memory allocations for a group in the data repository recursively
 * 
 * @param group The group to print
 * @param[in] indent The level of indentation to add to this level of output.
 * @param[in] threshold The allocation size required to output a allocation size.
 */
void printMemoryAllocation( Group const & group, integer const indent, real64 const threshold );

}
}