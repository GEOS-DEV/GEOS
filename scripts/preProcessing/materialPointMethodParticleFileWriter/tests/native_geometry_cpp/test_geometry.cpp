// Standalone numerical tests of the production header. The test driver supplies
// only DataTypes typedefs; this is NOT a GEOS library or device build.
#include "mesh/ParticleReferenceGeometry.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <random>
using namespace geos;

int main()
{
  std::mt19937 rng(972);
  std::uniform_real_distribution<double> d(-0.25,0.25);
  double maxError=0;
  for(int trial=0;trial<3000;++trial)
  {
    double oldF[3][3],newF[3][3],reference[3][3],copy[3][3],before[3][3],after[3][3];
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
    {
      oldF[i][j]=(i==j?1.:0.)+d(rng);
      newF[i][j]=(i==j?1.:0.)+d(rng);
      reference[i][j]=(i==j?1.:0.)+d(rng);
      copy[i][j]=reference[i][j];
    }
    particleGeometry::deform(reference,oldF,before);
    assert(particleGeometry::rebase(oldF,newF,reference));
    particleGeometry::deform(reference,newF,after);
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
      maxError=std::max(maxError,std::abs(before[i][j]-after[i][j]));
    for(int shape=0;shape<3;++shape)
    {
      double c[3]={1,2,3},v[8][3],mean[3]={};
      particleGeometry::vertices(shape,c,before,v);
      int n=particleGeometry::vertexCount(shape);
      for(int k=0;k<n;++k) for(int j=0;j<3;++j) mean[j]+=v[k][j]/n;
      for(int j=0;j<3;++j) assert(std::abs(mean[j]-c[j])<3e-14);
      double expected=particleGeometry::determinant(oldF)*particleGeometry::signedMeasure(shape,copy);
      assert(std::abs(particleGeometry::signedMeasure(shape,before)-expected)<1e-12);
    }
  }
  assert(maxError<1e-13);
  double I[3][3]={{1,0,0},{0,1,0},{0,0,1}},G[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  double singular[3][3]={};
  assert(!particleGeometry::rebase(I,singular,G));
  assert(G[0][0]==1 && G[1][1]==1 && G[2][2]==1);
  double bad[3][3]={{1,0,0},{0,1,0},{0,0,std::numeric_limits<double>::infinity()}};
  assert(!particleGeometry::rebase(I,bad,G));
  double scale[3][3]={{1e-90,0,0},{0,1e-90,0},{0,0,1e-90}};
  assert(particleGeometry::rebase(I,scale,G));
  double current[3][3];particleGeometry::deform(G,scale,current);
  assert(std::abs(current[0][0]-1)<1e-14);
  int fromSiloWedge[6]={2,1,5,3,0,4},fromSiloTet[4]={1,0,2,3};
  for(int i=0;i<6;++i) assert(particleGeometry::siloVertex(1,fromSiloWedge[i])==particleGeometry::vtkVertex(1,i));
  for(int i=0;i<4;++i) assert(particleGeometry::siloVertex(2,fromSiloTet[i])==i);
  std::cout<<"3000 affine/reset trials passed; maximum generator error="<<maxError<<"\n";
}
