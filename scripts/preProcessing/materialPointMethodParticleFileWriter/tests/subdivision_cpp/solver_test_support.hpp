// Test-only host views, wrappers and collectives. This is not GEOS/LvArray/MPI.
#include "mesh/ParticleSubdivision.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
using namespace geos;
using globalIndex=long long;
using string=std::string;
constexpr int hostMemorySpace=0, MPI_COMM_GEOS=0;
#define GEOS_MARK_FUNCTION
#define GEOS_ERROR_IF(c,m) do { if(c){std::ostringstream msg;msg<<m;throw std::runtime_error(msg.str());} } while(false)
#define GEOS_LOG_RANK_0(m) do {std::ostringstream msg;msg<<m;} while(false)
#define GEOS_FORCE_INLINE inline
struct serialPolicy {};
struct parallelDeviceReduce {};
template<class Policy,class F> void forAll(int n,F f){for(int i=0;i<n;++i)f(i);}

template<class T,int N> struct Slice
{
  T * data;std::array<int,N> dimensions;
  int size(int a=0)const{return dimensions[a];}
  decltype(auto) operator[](int index)const
  {
    assert(index>=0 && index<dimensions[0]);
    if constexpr(N==1)return (data[index]);
    else
    {
      int stride=1;std::array<int,N-1> next{};
      for(int a=1;a<N;++a){stride*=dimensions[a];next[a-1]=dimensions[a];}
      return Slice<T,N-1>{data+index*stride,next};
    }
  }
  template<class...I> decltype(auto) operator()(int i,I...indices)const
  {if constexpr(sizeof...(I)==0)return (*this)[i];else return (*this)[i](indices...);}
};
template<class T,int N> struct Array
{
  using ValueType=T;static constexpr int NDIM=N;using Raw=std::remove_const_t<T>;
  std::shared_ptr<std::vector<Raw>> storage=std::make_shared<std::vector<Raw>>();
  std::array<int,N> dimensions{};
  Array()=default;
  explicit Array(int n){dimensions.fill(1);resize(n);}
  template<class U> Array(Array<U,N> const & other):storage(other.storage),dimensions(other.dimensions){}
  template<class...I> void resize(int n,I...rest)
  {
    if constexpr(sizeof...(I)>0){static_assert(sizeof...(I)==N-1);dimensions={n,rest...};}
    else dimensions[0]=n;
    int count=1;for(int x:dimensions)count*=x;storage->resize(count);
  }
  int size(int a=0)const{return dimensions[a];}
  decltype(auto) operator[](int i)const{return Slice<T,N>{storage->data(),dimensions}[i];}
  Array toView()const{return *this;}
};
template<class T> using array1d=Array<T,1>;
template<class T> using arrayView1d=Array<T,1>;
template<class T> using arrayView2d=Array<T,2>;
template<class T> using arrayView3d=Array<T,3>;
template<class T> using SortedArrayView=Array<T,1>;
namespace RAJA
{
 template<class P,class T>struct ReduceSum{std::shared_ptr<T> v;explicit ReduceSum(T x):v(std::make_shared<T>(x)){}void operator+=(T x)const{*v+=x;}T get()const{return *v;}};
 template<class P,class T>struct ReduceMax:ReduceSum<P,T>{using ReduceSum<P,T>::ReduceSum;void max(T x)const{*this->v=std::max(*this->v,x);}};
}
namespace LvArray
{
 namespace math{using std::sqrt;using std::min;using std::max;using std::pow;using std::abs;}
 namespace tensorOps
 {
 template<int N,class A,class B>void copy(A&& a,B const& b){for(int i=0;i<N;++i)a[i]=b[i];}
 template<int N,int M,class A,class B>void copy(A&& a,B const& b){for(int i=0;i<N;++i)copy<M>(a[i],b[i]);}
 template<int N,class A>double l2NormSquared(A const& a){double s=0;for(int i=0;i<N;++i)s+=a[i]*a[i];return s;}
 template<int N,class A,class B>void add(A&& a,B const& b){for(int i=0;i<N;++i)a[i]+=b[i];}
 template<int N,class A,class B>void subtract(A&& a,B const& b){for(int i=0;i<N;++i)a[i]-=b[i];}
 template<int N,class A>void scale(A&& a,double x){for(int i=0;i<N;++i)a[i]*=x;}
 template<int N,class A>void fill(A&& a,double x){for(int i=0;i<N;++i)a[i]=x;}
 template<int N,int M,class A,class B,class C>void Ri_eq_AijBj(A&& a,B const& b,C const& c){for(int i=0;i<N;++i){a[i]=0;for(int j=0;j<M;++j)a[i]+=b[i][j]*c[j];}}
 template<int N,class A>double determinant(A const& a){static_assert(N==3);return particleGeometry::affine::determinant(a);}
 template<int N,class A,class B>void invert(A&& a,B const& b){static_assert(N==3);double d=determinant<3>(b);for(int i=0;i<3;++i)for(int j=0;j<3;++j)a[i][j]=(b[(j+1)%3][(i+1)%3]*b[(j+2)%3][(i+2)%3]-b[(j+1)%3][(i+2)%3]*b[(j+2)%3][(i+1)%3])/d;}
 }
 template<class T,int N,class F,class...I>void visit(Slice<T,N> s,F f,I...indices){for(int i=0;i<s.size();++i){if constexpr(N==1)f(s[i],indices...,i);else visit(s[i],f,indices...,i);}}
 template<class T,int N,class F>void forValuesInSliceWithIndices(Slice<T,N> s,F f){visit(s,f);}
}
struct WrapperBase
{
 string name;bool sized=true;explicit WrapperBase(string n):name(std::move(n)){}virtual ~WrapperBase()=default;
 string getName()const{return name;}bool sizedFromParent()const{return sized;}void move(int,bool){}
 virtual int numArrayDims()const=0;virtual void resize(int)=0;
};
template<class A>struct Wrapper:WrapperBase
{
 A value;explicit Wrapper(string n):WrapperBase(std::move(n)){}
 static Wrapper & cast(WrapperBase & w){return dynamic_cast<Wrapper&>(w);}A & reference(){return value;}
 int numArrayDims()const override{return A::NDIM;}void resize(int n)override{if(sized)value.resize(n);}
};
namespace camp{template<class T>using first=std::tuple_element_t<0,T>;}
namespace types
{
 struct StandardArrays{};template<class T>struct ListofTypeList{};
 template<class A,class F>bool attempt(WrapperBase & w,F f){if(dynamic_cast<Wrapper<A>*>(&w)){f(std::tuple<A>{});return true;}return false;}
 template<class F>void dispatch(ListofTypeList<StandardArrays>,F f,WrapperBase & w)
 {
   bool ok=attempt<Array<int,1>>(w,f)||attempt<Array<int,2>>(w,f)||attempt<Array<globalIndex,1>>(w,f)||attempt<Array<double,1>>(w,f)||attempt<Array<double,2>>(w,f)||attempt<Array<double,3>>(w,f)||attempt<Array<double,4>>(w,f);
   assert(ok);
 }
}
struct Fields
{
 std::map<string,std::unique_ptr<WrapperBase>> wrappers;int count=0;
 template<class A,class...I>A add(string const & name,I...dims){auto w=std::make_unique<Wrapper<A>>(name);w->value.resize(count,dims...);A a=w->value;wrappers[name]=std::move(w);return a;}
 template<class A>A get(string const & name){return Wrapper<A>::cast(*wrappers.at(name)).value;}
 WrapperBase & getWrapperBase(string const & name){return *wrappers.at(name);}
 template<class F>void forWrappers(F f){for(auto & item:wrappers)f(*item.second);}
 void resize(int n){count=n;forWrappers([&](WrapperBase & w){w.resize(n);});}
};
namespace mpm
{
 enum class OversizedParticleTreatmentOption{None,Split,ResetDeformationGradient};
 enum class CPDIDomainScalingTypeOption{Homel,VPHencky};enum class DomainResetTypeOption{Isotropic,VPHencky};
 enum class SurfaceFlag{Interior=0,Surface=1,Cohesive=2,DamagedCohesive=3};
 template<class T>int toInteger(T x){return static_cast<int>(x);}
}
enum class ParticleType{SinglePoint,SinglePointBSpline,CPDI,CPTI};
struct ContinuumBase:Fields{string catalog="Solid";string getCatalogName()const{return catalog;}};
// FIELD_TAGS
struct ParticleSubRegion:Fields
{
 ParticleType type=ParticleType::SinglePointBSpline;ContinuumBase model;string material="solid";Array<int,1> active;
 struct Neighbors{int count=0;void resize(int n){count=n;}}neighbors;
 int size()const{return count;}ParticleType getParticleType()const{return type;}
 template<class F>typename F::Type getField(){return get<typename F::Type>(F::key());}
 template<class T>T const& getReference(string const&)const{static_assert(std::is_same_v<T,string>);return material;}
 void resize(int n){Fields::resize(n);model.resize(n);}
 void setActiveParticleIndices(){auto rank=get<Array<int,1>>("particleRank");active.resize(count);int n=0;for(int p=0;p<count;++p)if(rank[p]>=0)active[n++]=p;active.resize(n);}
 auto activeParticleIndices()const{return active;}
 auto & neighborList(){return neighbors;}
 auto getParticleID(){return get<Array<globalIndex,1>>("particleID");}
 auto getParticleSurfaceFlag(){return get<Array<int,1>>("particleSurfaceFlag");}
 auto getParticleDamage(){return get<Array<double,1>>("particleDamage");}
 auto getParticleVolume(){return get<Array<double,1>>("particleVolume");}
 auto getParticleGeometryShape(){return get<Array<int,1>>("particleGeometryShape");}
 auto getParticleGeometryCurrent(){return get<Array<double,3>>("particleGeometryCurrent");}
 auto getParticleCenter(){return get<Array<double,2>>("particleCenter");}
 auto getParticleVelocity(){return get<Array<double,2>>("particleVelocity");}
 auto getParticleSurfacePosition(){return get<Array<double,2>>("particleSurfacePosition");}
 auto getParticleSurfaceTraction(){return get<Array<double,2>>("particleSurfaceTraction");}
 struct GeometryView
 {
   Array<int,1> shapes;Array<double,3> current;
   auto particleShape(int p)const{return static_cast<particleGeometry::Shape>(shapes[p]);}
   auto currentDomain(int p,int dimensions)const{auto shape=particleShape(p);if(dimensions==2&&shape==particleGeometry::Shape::Parallelepiped)shape=particleGeometry::Shape::Parallelogram;if(dimensions==2&&shape==particleGeometry::Shape::TriangleWedge)shape=particleGeometry::Shape::Triangle;return particleGeometry::affine::Domain::fromGenerators(shape,current[p]);}
 };
 auto geometryViewOnHost(){return GeometryView{getParticleGeometryShape(),getParticleGeometryCurrent()};}
};
struct ParticleManager
{
 std::vector<std::unique_ptr<ParticleSubRegion>> regions;
 template<class F>void forParticleSubRegions(F f){for(auto & r:regions)f(*r);}
};
template<class T>T& getConstitutiveModel(ParticleSubRegion & r,string const&){return r.model;}
namespace MpiWrapper
{
 inline int rank=0;inline std::vector<int> counts;inline globalIndex globalMaximum=0;
 int commSize(){return counts.empty()?1:static_cast<int>(counts.size());}int commRank(int){return rank;}
 globalIndex max(globalIndex v){return std::max(v,globalMaximum);}int sum(int v){return v;}
 void allGather(int local,Array<int,1> out,int){if(counts.empty())out[0]=local;else{assert(counts[rank]==local);for(int i=0;i<out.size();++i)out[i]=counts[i];}}
}
