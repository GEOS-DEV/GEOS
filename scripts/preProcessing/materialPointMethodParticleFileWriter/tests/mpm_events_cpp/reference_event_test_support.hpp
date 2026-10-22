// CPU repository adapters for the production ReferenceCohesiveZones event.
// Scheduling, activation and the event class itself are compiled from GEOS.
#pragma once
#include <cassert>
#include <cfloat>
#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#define GEOS_ERROR_IF(condition, message) do { if(condition) { \
  std::ostringstream stream; stream << message; throw std::runtime_error(stream.str()); } } while(false)
#define GEOS_LOG_RANK_0_IF(condition, message) ((void)(condition))
#define GEOS_LOG_RANK_0(message) ((void)0)
#define REGISTER_CATALOG_ENTRY(...)

namespace geos
{
using real64 = double;
using string = std::string;
using string_array = std::vector<string>;
template<class T, class U> T dynamicCast(U & value) { return dynamic_cast<T>(value); }
namespace dataRepository
{
enum class InputFlags { FALSE, OPTIONAL, REQUIRED };
enum class RestartFlags { NO_WRITE, WRITE_AND_READ };
struct Wrapper
{
  InputFlags input = InputFlags::OPTIONAL;
  RestartFlags restart = RestartFlags::NO_WRITE;
  std::function<int()> read;
  std::function<void(int)> write;
  Wrapper & setInputFlag(InputFlags value) { input=value; return *this; }
  Wrapper & setRestartFlags(RestartFlags value) { restart=value; return *this; }
  template<class T> Wrapper & setApplyDefaultValue(T const &) { return *this; }
  Wrapper & setDescription(string const &) { return *this; }
};
struct Group
{
  virtual ~Group() = default;
  std::map<string, Wrapper> wrappers;
  template<class T> Wrapper & registerWrapper(string const & key, T *value)
  {
    auto &wrapper=wrappers[key];
    if constexpr(std::is_same_v<T,int>)
    {
      wrapper.read=[value]() { return *value; };
      wrapper.write=[value](int v) { *value=v; };
    }
    return wrapper;
  }
};
}
using dataRepository::Group;
namespace mpm
{
enum class NormalsAndPositionsMethodOption { LogisticRegression };
enum class CohesiveSurfaceDisplacementUpdateOption { Nodal, TypeA, TypeB };
}
template<class T> struct EnumStrings { static string concat(string const &) { return "test enums"; } };
class MPMEventBase : public Group
{
public:
  MPMEventBase(string const & name, Group const *) : name_(name) {}
  virtual string getCatalogName() const { return "Dependency"; }
  virtual void postInputInitialization() {}
  string const & getName() const { return name_; }
  bool hasDependencies() const { return !dependencies.empty(); }
  string_array getDependencies() const { return dependencies; }
  real64 getDelay() const { return delay; }
  real64 getDuration() const { return duration; }
  void setStartTime(real64 t) { start=t; }
  void setEndTime(real64 t) { end=t; }
  real64 getStartTime() const { return start; }
  real64 getEndTime() const { return end; }
  real64 getTimeInterval() const { return end-start; }
  int hasStarted() const { return started; }
  void setHasStarted(int flag) { started=flag; }
  int isComplete() const { return complete; }
  void setIsComplete(int flag) { complete=flag; }
  string_array dependencies;
  real64 start=0, end=DBL_MAX, delay=0, duration=1;
  int started=0, complete=0;
private:
  string name_;
};
struct EventManager
{
  MPMEventBase dependency{"ready",nullptr};
  template<class T> T & getEvent(string const &) { return dependency; }
};
struct CohesiveZoneRegion
{
  int enabled=0, initialized=0, referenceRequests=0, displacementMode=0;
  void setEnabled(int flag) { enabled=flag; }
  void setInitialized(int flag) { initialized=flag; if(!flag) ++referenceRequests; }
  void setCZSurfaceDisplacementUpdate(int mode) { displacementMode=mode; }
};
struct CohesiveZoneManager
{
  std::map<string,CohesiveZoneRegion> regions{{"bond",{}},{"second",{}}};
  template<class T> T & getRegion(string const & name) { return regions.at(name); }
  template<class T> bool hasRegion(string const & name) const { return regions.count(name)>0; }
};
struct EventState
{
  real64 start, end;
  int started, complete;
  std::map<string,int> registered;
};
// Exercise the production wrapper's restart registration. This replaces the
// repository serializer only; it is not a full GEOS/HDF5 restart test.
inline EventState save(MPMEventBase const & event)
{
  EventState state{event.start,event.end,event.started,event.complete,{}};
  for(auto const & entry:event.wrappers)
    if(entry.second.restart==dataRepository::RestartFlags::WRITE_AND_READ && entry.second.read)
      state.registered[entry.first]=entry.second.read();
  return state;
}
inline void restore(MPMEventBase & event, EventState const & state)
{
  event.start=state.start; event.end=state.end; event.started=state.started; event.complete=state.complete;
  for(auto const & entry:state.registered) event.wrappers.at(entry.first).write(entry.second);
}
}
