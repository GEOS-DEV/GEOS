#ifndef RESTARTFLAGS_H_
#define RESTARTFLAGS_H_

namespace geosx
{
namespace dataRepository
{

enum class RestartFlags : unsigned char
{
  NO_WRITE,
  WRITE,
  WRITE_AND_READ
};

}   /* namespace dataRepository */
}   /* namespace geosx */

#endif  /* RESTARTFLAGS_H_ */
