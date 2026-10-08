/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MPM_MORTAR_STATE_VALIDATION_HPP
#define GEOS_MPM_MORTAR_STATE_VALIDATION_HPP

#include <algorithm>
#include <cstddef>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace geos { namespace mpm { namespace mortarState {
using Key=std::array<long long,3>;

/** The global, uniform interpolation lattice, not a rank's floating-point copy.
 * Morton/node IDs from MPI are NOT lattice coordinates. A valid physical buffer
 * is in [0,N]; a periodic key is in [0,N). An absent record is never zero mass.
 */
struct Lattice
{
  std::array<double,3> low{},high{},spacing{};
  Key cells{};
  std::array<bool,3> periodic{{false,false,false}};

  static Lattice make(std::array<double,3> const & lo, std::array<double,3> const & hi,
                      Key const & n, std::array<bool,3> const & p)
  {
    Lattice out;out.low=lo;out.high=hi;out.cells=n;out.periodic=p;
    for(int d=0;d<3;++d)
    {
      if(!std::isfinite(lo[d])||!std::isfinite(hi[d])||!(hi[d]>lo[d])||n[d]<1||(!p[d]&&n[d]<3))
        throw std::invalid_argument("Invalid global mortar lattice; a nonperiodic B-spline direction needs >=3 cells.");
      out.spacing[d]=(hi[d]-lo[d])/static_cast<double>(n[d]);
    }
    return out;
  }
  bool contains(Key const & key) const
  {
    for(int d=0;d<3;++d)
      if(key[d]<0||key[d]>cells[d]||(periodic[d]&&key[d]==cells[d]))return false;
    return true;
  }
  Key key(std::array<double,3>const & x)const
  {
    Key out{};
    for(int d=0;d<3;++d)
    {
      double u=(x[d]-low[d])/spacing[d];
      if(!std::isfinite(u)||std::abs(u)>9.e15)throw std::invalid_argument("Unrepresentable mortar node coordinate.");
      out[d]=std::llround(u);
      if(periodic[d])out[d]=(out[d]%cells[d]+cells[d])%cells[d];
    }
    return out;
  }
  std::array<double,3> point(Key const & key)const
  {std::array<double,3> out{};for(int d=0;d<3;++d)out[d]=low[d]+key[d]*spacing[d];return out;}
  // Snap only the mortar's geometric working copy. Preserve unwrapped periodic
  // images here so tetrahedra at a periodic partition face do not span the box.
  std::array<double,3> snap(std::array<double,3>const & x)const
  {
    std::array<double,3> out{};
    for(int d=0;d<3;++d)out[d]=low[d]+std::round((x[d]-low[d])/spacing[d])*spacing[d];
    return out;
  }
  double mismatch(std::array<double,3>const & x)const
  {auto y=snap(x);double e=0;for(int d=0;d<3;++d)e=std::max(e,std::abs(x[d]-y[d])/spacing[d]);return e;}
};

inline std::set<Key> requestUnion(std::vector<Key>const & allRequests)
{return std::set<Key>(allRequests.begin(),allRequests.end());}

/** Stable fingerprint over individually encoded values; never hash padded structs.
 * This is an audit of replicated inputs before component workers run, not a
 * replacement for the physical residual checks after solving.
 */
struct Fingerprint
{
  unsigned long long value=1469598103934665603ULL;
  void integer(unsigned long long x)
  {for(int i=0;i<8;++i){value^=(x>>(i*8))&255ULL;value*=1099511628211ULL;}}
  void scalar(double x)
  {if( std::fpclassify( x ) == FP_ZERO ) { x = 0.0; }std::uint64_t bits=0;std::memcpy(&bits,&x,sizeof bits);integer(bits);}
  void key(Key const & k){for(auto v:k)integer(static_cast<unsigned long long>(v));}
};

inline std::string keyText(Key const & k)
{std::ostringstream s;s<<"("<<k[0]<<","<<k[1]<<","<<k[2]<<")";return s.str();}

struct Missing
{Key key{};int field=0;bool addedByClosure=false;};

// priority=-1 means absent, NOT inactive. A genuine inactive/zero-mass record has
// an owner/ghost priority >=0 and must pass this check unchanged.
inline std::vector<Missing> missingRecords(std::vector<Key> const & keys,
    std::size_t fields, std::vector<int>const & participating,
    std::vector<int>const & priority, std::set<Key>const & originalRequests)
{
  if(fields==0||priority.size()!=keys.size()*fields)throw std::invalid_argument("Invalid mortar record dimensions.");
  std::vector<Missing> out;
  for(std::size_t i=0;i<keys.size();++i)for(int field:participating)
  {
    if(field<0||static_cast<std::size_t>(field)>=fields)throw std::invalid_argument("Invalid mortar field.");
    if(priority[i*fields+field]<0)out.push_back(Missing{keys[i],field,originalRequests.count(keys[i])==0});
  }
  return out;
}

}}} // geos::mpm::mortarState
#endif
