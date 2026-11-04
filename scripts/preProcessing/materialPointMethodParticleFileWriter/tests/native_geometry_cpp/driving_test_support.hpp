// CPU-only repository/view, reduction and MPI adapters for source-excerpt tests.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#define GEOS_HOST_DEVICE
#define GEOS_FORCE_INLINE inline
#define GEOS_MARK_FUNCTION
#define GEOS_UNUSED_PARAM(x) x [[maybe_unused]]
#define GEOS_UNUSED_VAR(x) (void)(x)
#define GEOS_ERROR_IF(c,m) do {if(c) throw std::runtime_error("rejected");}while(false)
#define GEOS_ERROR(m) throw std::runtime_error("error")
#define GEOS_LOG_RANK_0(m) do {}while(false)
using real64=double;using localIndex=int;using integer=int;
template<class T> struct View {
 T *data;int count;
 template<class V> View(V &v):data(v.data()),count(int(v.size())){}
 T &operator[](int i) const{return data[i];}int size()const{return count;}
};
template<class T> using arrayView1d=View<std::remove_const_t<T>>;
template<class T> struct Array:std::vector<T> {
 using std::vector<T>::vector;
 View<T> toViewConst(){return View<T>(*this);}
};
template<class T> using array1d=Array<T>;
struct Table {
 std::shared_ptr<std::vector<std::vector<double>>> data=std::make_shared<std::vector<std::vector<double>>>();
 Table()=default;
 Table(std::initializer_list<std::vector<double>> rows):data(std::make_shared<std::vector<std::vector<double>>>(rows)){}
 int size(int dim=0)const{return dim==0?int(data->size()):(data->empty()?0:int(data->front().size()));}
 std::vector<double> &operator[](int i)const{return (*data)[i];}
};
template<class T,int=0> using arrayView2d=Table;
template<class T> using arrayView3d=T(*)[1][3];
template<class T> using SortedArrayView=View<std::remove_const_t<T>>;
namespace LvArray {
 namespace math {using std::abs;using std::sqrt;using std::min;using std::max;using std::cos;using std::sin;using std::pow;using std::round;}
 template<class T> struct NumericLimits {static constexpr T max=std::numeric_limits<T>::max();};
 namespace tensorOps {
 template<int N,class A,class B> void copy(A &a,B const &b){for(int i=0;i<N;++i)a[i]=b[i];}
 template<int N,class A,class T> void fill(A &a,T v){for(int i=0;i<N;++i)a[i]=v;}
 template<int N,class A,class B> void subtract(A &a,B const &b){for(int i=0;i<N;++i)a[i]-=b[i];}
 template<int N,class A,class B> void add(A &a,B const &b){for(int i=0;i<N;++i)a[i]+=b[i];}
 template<int N,class A> void scale(A &a,double v){for(int i=0;i<N;++i)a[i]*=v;}
 }
}
bool isZero(double x){return std::abs(x)<1.e-14;}
bool isEqual(double a,double b){return isZero(a-b);}
double clampValue(double x,double lo,double hi){return std::clamp(x,lo,hi);}
struct serialPolicy{};struct parallelDeviceReduce{};
template<class P,class F> void forAll(int n,F f){for(int i=0;i<n;++i)f(i);}
namespace mpm {enum class UpdateMethodOption{FLIP,FMPM};enum class InterpolationOption{Linear,Cosine,Smoothstep};enum class BoundaryConditionOption{Outflow,Symmetry,Moving,Contact};}
namespace nodes {constexpr int REFERENCE_POSITION_USD=0;}
namespace RAJA {
 struct seq_multi_reduce{};
 template<class P,class T> struct ReduceSum {std::shared_ptr<T> v;explicit ReduceSum(T x):v(std::make_shared<T>(x)){} void operator+=(T x)const{*v+=x;}T get()const{return *v;}};
 template<class P,class T> struct MultiReduceSum {std::vector<ReduceSum<P,T>> v;MultiReduceSum(int n,T x):v(n,ReduceSum<P,T>(x)){for(auto &item:v)item=ReduceSum<P,T>(x);}auto const &operator[](int i)const{return v[i];}};
}
constexpr int MPI_DOUBLE=0,MPI_SUM=0,MPI_COMM_GEOS=0;
void MPI_Allreduce(double const *a,double *b,int,int,int,int){*b=*a;}
namespace fields {namespace mpm {struct particleStress{};}}
struct ParticleSubRegion {
 Array<double> volume;Table stress;Array<int> indices;
 arrayView1d<double const> getParticleVolume(){return volume;}
 template<class T> Table getField(){return stress;}
 SortedArrayView<int const> activeParticleIndices(){return indices;}
};
struct ParticleManager {double bulk=1000;std::vector<ParticleSubRegion> regions;template<class F>void forParticleSubRegions(F f){for(auto &r:regions)f(r);}};
struct SpatialPartition {double ratio[3]={1,1,1};void updateSizes(Array<double>const &rate,double dt){for(int i=0;i<3;++i)ratio[i]*=1+rate[i]*dt;}};
struct NodeManager {Table positions{{0,0,0},{1,1,0}};Table referencePosition(){return positions;}};
struct Wall {
 bool axisymmetric=true;int boundaryType=3,dir0=0,dir1=1,dir2=2,positiveNormal=1,face=1;
 int numDims=2,g=0,fieldIndex=0,useInternalForceAsFaceReaction=0;
 double dt=.1,smallMass=1e-12,faceReaction=0;
 double domainIncrementalStretch[3]={1,1,1},domainL[3]={};
 double gridPosition[1][3]={{1,1,0}},gridVelocity[1][1][3]={{{1,2,0}}},gridAcceleration[1][1][3]={};
 double gridCenterOfVolume[1][1][3]={{{.01,.01,0}}},gridSurfacePosition[1][1][3]={};
 double gridSurfaceFieldMass[1][1]={},gridMass[1][1]={{10}},gridInternalForce[1][1][3]={};int gridGhostRank[1]={-1};
 double boundaryFaceFrictionCoefficients[6]={.2,.2,.2,.2,.2,.2};
 int enablePrescribedBoundaryTransverseVelocities[6]={};double prescribedBoundaryTransverseVelocities[6][2]={};
 int targetKind[1][1][3]={},targetCount[1][1][3]={},targetDVelocityMode[1][1][3]={};
 double targetVelocity[1][1][3]={},targetAcceleration[1][1][3]={},targetDVelocity[1][1][3]={};
 void apply();
};
