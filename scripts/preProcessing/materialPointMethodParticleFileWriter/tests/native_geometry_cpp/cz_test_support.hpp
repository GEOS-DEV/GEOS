// CPU-only array and tensor adapters for extracted production CZ routines.
#include "physicsSolvers/solidMechanics/MPMGeometryMeasures.hpp"
#include "mesh/ParticleGeometry.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <random>
using namespace geos;
using globalIndex=long long;
template<class T> using arrayView2d=T(*)[2];
template<class T> using arrayView3d=T(*)[2][3];
template<class T> using arrayView4d=T(*)[2][3][3];
struct Serial {};
template<class P,class F> void forAll(int n,F f) {for(int i=0;i<n;++i) f(i);}
#define GEOS_UNUSED_VAR(x) (void)(x)
#define GEOS_ERROR_IF(x,msg) assert(!(x))
bool isZero(double x) {return std::abs(x)<1e-30;}
namespace mpm {enum class SurfaceFlag {Cohesive=3,DamagedCohesive=4};
int toInteger(SurfaceFlag f) {return static_cast<int>(f);} }
namespace LvArray {
template<class T,class U> T integerConversion(U x) {return static_cast<T>(x);}
namespace math {using std::abs; using std::floor;}
namespace tensorOps {
template<int N,class A> void fill(A&& a,double v) {for(int i=0;i<N;++i) a[i]=v;}
template<int N,class A,class B> void copy(A&& a,B const& b) {for(int i=0;i<N;++i) a[i]=b[i];}
template<int N,int M,class A,class B> void copy(A&& a,B const& b) {for(int i=0;i<N;++i) for(int j=0;j<M;++j) a[i][j]=b[i][j];}
template<int N,class A,class B> double AiBi(A const&a,B const&b) {double s=0;for(int i=0;i<N;++i)s+=a[i]*b[i];return s;}
template<int N,class A> double l2NormSquared(A const&a) {return AiBi<N>(a,a);}
template<int N,class A> double l2Norm(A const&a) {return std::sqrt(l2NormSquared<N>(a));}
template<int N,class A> void scale(A&&a,double s) {for(int i=0;i<N;++i)a[i]*=s;}
template<int N,class A,class B> void scaledCopy(A&&a,B const&b,double s) {for(int i=0;i<N;++i)a[i]=s*b[i];}
template<int N,class A,class B> void scaledAdd(A&&a,B const&b,double s) {for(int i=0;i<N;++i)a[i]+=s*b[i];}
template<int N,class A,class B> void subtract(A&&a,B const&b) {scaledAdd<N>(a,b,-1);}
template<int N,int M,class R,class A,class B> void Ri_eq_AijBj(R&&r,A const&a,B const&b) {for(int i=0;i<N;++i){r[i]=0;for(int j=0;j<M;++j)r[i]+=a[i][j]*b[j];}}
template<int N,int M,int K,class R,class A,class B> void Rij_eq_AikBkj(R&&r,A const&a,B const&b) {for(int i=0;i<N;++i)for(int j=0;j<M;++j){r[i][j]=0;for(int k=0;k<K;++k)r[i][j]+=a[i][k]*b[k][j];}}
template<int N,class R,class A> void cofactor(R&&r,A const&a) {static_assert(N==3);for(int i=0;i<3;++i)for(int j=0;j<3;++j)r[i][j]=a[(i+1)%3][(j+1)%3]*a[(i+2)%3][(j+2)%3]-a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3];}
template<int N,class A> double invert(A&&a) {double c[3][3];cofactor<N>(c,a);double d=AiBi<3>(a[0],c[0]);for(int i=0;i<3;++i)for(int j=0;j<3;++j)a[i][j]=c[j][i]/d;return d;}
}}
void close(double a,double b,double tol=2e-11) {if(std::abs(a-b)>2e-14+tol*std::max(std::abs(a),std::abs(b))){std::cerr<<a<<" != "<<b<<"\n";std::abort();}}
