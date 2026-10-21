// Test-only CPU adapters, not solver implementations.
#include "mesh/ParticleGeometry.hpp"
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <random>
using namespace geos;
template<class T> using arraySlice1d=T*;
template<class T> using arraySlice2d=T(*)[3];
template<class T,int USD=1> using arrayView2d=T(*)[3];
template<class T> using arrayView3d=T(*)[9][9];
namespace nodes { constexpr int REFERENCE_POSITION_USD=1; }
namespace mpm { enum class CPDIDomainScalingTypeOption { Homel,VPHencky }; }
namespace LvArray {
namespace math {
using std::abs; using std::min; using std::max; using std::pow; using std::sqrt;
using std::log; using std::exp; using std::atan2; using std::cos; using std::sin; using std::floor;
}
namespace tensorOps {
template<int N,int M,class A,class B> void copy(A&& a,B const& b)
{ for(int i=0;i<N;++i) for(int j=0;j<M;++j) a[i][j]=b[i][j]; }
template<int N,class A> double l2NormSquared(A const& a)
{ double s=0; for(int i=0;i<N;++i) s+=a[i]*a[i]; return s; }
template<int N,class A> double l2Norm(A const& a) {return std::sqrt(l2NormSquared<N>(a));}
template<int N,class A> void scale(A&& a,double s) {for(int i=0;i<N;++i) a[i]*=s;}
template<int N,class A> double determinant(A const& a) {static_assert(N==3);return particleGeometry::affine::determinant(a);}
template<int N,class A,class B> void denseToSymmetric(A&& a,B const& b)
{static_assert(N==3);a[0]=b[0][0];a[1]=b[1][1];a[2]=b[2][2];a[3]=b[1][2];a[4]=b[0][2];a[5]=b[0][1];}
template<int N,class E,class V,class A> void symEigenvectors(E&& e,V&& v,A const& a)
{
  static_assert(N==3);
  double m[3][3]={{a[0],a[5],a[4]},{a[5],a[1],a[3]},{a[4],a[3],a[2]}};
  double q[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  for(int it=0;it<60;++it)
  {
    int p=0,r=1;
    for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) if(std::abs(m[i][j])>std::abs(m[p][r])) {p=i;r=j;}
    if(std::abs(m[p][r])<1e-15*(std::abs(m[0][0])+std::abs(m[1][1])+std::abs(m[2][2]))) break;
    double theta=.5*std::atan2(2*m[p][r],m[r][r]-m[p][p]),c=std::cos(theta),s=std::sin(theta);
    double R[3][3]={{1,0,0},{0,1,0},{0,0,1}},tmp[3][3]={},next[3][3]={},nq[3][3]={};
    R[p][p]=c;R[r][r]=c;R[p][r]=s;R[r][p]=-s;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) for(int k=0;k<3;++k) {tmp[i][j]+=m[i][k]*R[k][j];nq[i][j]+=q[i][k]*R[k][j];}
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) for(int k=0;k<3;++k) next[i][j]+=R[k][i]*tmp[k][j];
    copy<3,3>(m,next);copy<3,3>(q,nq);
  }
  for(int i=0;i<3;++i) {e[i]=m[i][i];for(int j=0;j<3;++j) v[i][j]=q[j][i];}
}
}}
