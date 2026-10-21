// Production math tested independently of GEOS registration, MPI and devices.
#include "mesh/ParticleGeometry.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <random>
#include <string>
using namespace geos;
using particleGeometry::Shape;
using particleGeometry::affine::Domain;
using particleGeometry::affine::InputError;
using particleGeometry::affine::validateFileGeometry;

static bool close( double a, double b )
{ return std::abs(a-b) < 2e-12*std::max({1.0,std::abs(a),std::abs(b)}); }

int main( int argc, char ** argv )
{
  // Also accept real PFW records on stdin. The Python regression drives this
  // path to check the actual C++ reader codec, rather than a Python replica.
  if( (argc == 2 || argc == 3) && std::string(argv[1]) == "--decode" )
  {
    double const thickness = argc == 3 ? std::stod(argv[2]) : 1.0;
    double shape, volume;
    while( std::cin >> shape >> volume )
    {
      double g[3][3] = {};
      int const dim = shape >= 3 ? 2 : 3;
      for( int a=0; a<dim; ++a ) for( int i=0; i<dim; ++i ) std::cin >> g[a][i];
      assert(std::cin);
      if(validateFileGeometry(shape,volume,g,thickness)!=InputError::None) return 2;
      std::cout << "valid\n";
    }
    return 0;
  }

  std::mt19937 rng(1941);
  std::uniform_real_distribution<double> random(-0.3,0.3);
  int childrenChecked=0;
  for( int trial=0; trial<1000; ++trial )
  {
    double g[3][3], F[3][3];
    for( int a=0; a<3; ++a ) for( int j=0; j<3; ++j )
    {
      g[a][j]=(a==j?1.0:0.0)+random(rng);
      F[a][j]=(a==j?1.0:0.0)+random(rng);
    }
    for( Shape shape : {Shape::Parallelepiped,Shape::TriangleWedge,Shape::Tetrahedron,
                       Shape::Parallelogram,Shape::Triangle} )
    {
      auto parent=Domain::fromGenerators(shape,g);
      // A true 2D reference lives in xy; 3D components remain available in the API.
      if( particleGeometry::affine::dimension(shape)==2 )
        for(int a=0;a<3;++a) {parent.generators[a][2]=0; parent.generators[2][a]=0;}
      double center[3]={1.2,-0.4,2.1}, vertices[8][3], mean[3]={};
      parent.vertices(center,vertices);
      int n=parent.verticesCount();
      double largest=0;
      for(int a=0;a<n;++a)
      {
        for(int j=0;j<3;++j) mean[j]+=vertices[a][j]/n;
        for(int b=a+1;b<n;++b)
        {
          double d2=0;
          for(int j=0;j<3;++j) d2+=(vertices[a][j]-vertices[b][j])*(vertices[a][j]-vertices[b][j]);
          largest=std::max(largest,d2);
        }
      }
      for(int j=0;j<3;++j) assert(close(mean[j],center[j]));
      assert(close(parent.diameter(),std::sqrt(largest)));
      if(shape==Shape::Parallelepiped)
      {
        int const order[8]={6,2,5,1,7,3,4,0};
        for(int k=0;k<8;++k) for(int i=0;i<3;++i)
        {
          double expected=center[i];
          for(int a=0;a<3;++a) expected+=(k&(1<<(2-a))?-0.5:0.5)*g[a][i];
          assert(close(vertices[order[k]][i],expected));
        }
      }

      // Affine covariance: transform the existing vertex, and independently
      // reconstruct the transformed domain. In 3D volume scales by det(F).
      Domain deformed=parent;
      particleGeometry::affine::deform(parent.generators,F,deformed.generators);
      double newCenter[3]={}, next[8][3];
      for(int i=0;i<3;++i) for(int j=0;j<3;++j) newCenter[i]+=F[i][j]*center[j];
      deformed.vertices(newCenter,next);
      for(int a=0;a<n;++a) for(int i=0;i<3;++i)
      {
        double expected=0;
        for(int j=0;j<3;++j) expected+=F[i][j]*vertices[a][j];
        assert(close(expected,next[a][i]));
      }
      if(particleGeometry::affine::dimension(shape)==3)
        assert(close(deformed.signedMeasure(),parent.signedMeasure()*particleGeometry::affine::determinant(F)));

      // A ray from the centroid through each vertex hits that vertex at
      // t=1, for every shape, shear and orientation. This exercises CZ anchors.
      for(int k=0;k<n;++k)
      {
        double direction[3]; parent.vertexOffset(k,direction);
        assert(close(parent.rayDistance(direction),1));
      }
      double zero[3]={}; assert(parent.rayDistance(zero)==0);

      // Reset policy and subsequent convection. FOnly leaves G and G0 intact;
      // FAndGeometry selects newF*G0. Neither policy rebases reference edges.
      for(bool resetGeometry : {false,true})
      {
        auto current=deformed;
        double newF[3][3]={{1.1,.03,0},{0,.9,0},{0,0,1}};
        particleGeometry::affine::reset(parent.generators,current.generators,newF,resetGeometry);
        auto expected=resetGeometry ? parent : deformed;
        if(resetGeometry)
          particleGeometry::affine::deform(parent.generators,newF,expected.generators);
        for(int step=0;step<3;++step)
        {
          double expectedNext[3][3];
          particleGeometry::affine::deform(expected.generators,F,expectedNext);
          particleGeometry::affine::advance(current.generators,F);
          for(int a=0;a<3;++a) for(int i=0;i<3;++i)
          {
            assert(close(current.generators[a][i],expectedNext[a][i]));
            expected.generators[a][i]=expectedNext[a][i];
          }
        }
      }

      // All nonempty direction masks: child volume and first moments sum to
      // the parent; every child's edges/offsets agree with legacy bisection.
      bool const box=shape==Shape::Parallelepiped || shape==Shape::Parallelogram;
      int const dimensions=particleGeometry::affine::dimension(shape);
      for(int mask=1;mask<(1<<dimensions);++mask)
      {
        double volume=0, moment[3]={};
        int directions=0;
        for(int a=0;a<dimensions;++a) directions+=(mask>>a)&1;
        for(int side=0;side<(1<<dimensions);++side) if((side & ~mask)==0)
        {
          Domain child; double offset[3];
          bool const valid=parent.boxChild(mask,side,child,offset);
          assert(valid==box);
          if(!box) continue;
          ++childrenChecked;
          assert(close(child.measure()*(1<<directions),parent.measure()));
          volume+=child.measure();
          for(int i=0;i<3;++i)
          {
            double expected=0;
            for(int a=0;a<dimensions;++a) if(mask&(1<<a))
              expected+=(side&(1<<a)?1.0:-1.0)*0.25*parent.generators[a][i];
            assert(close(offset[i],expected));
            moment[i]+=child.measure()*(center[i]+offset[i]);
          }
          for(int a=0;a<3;++a) for(int i=0;i<3;++i)
            assert(close(child.generators[a][i],parent.generators[a][i]*(mask&(1<<a)?0.5:1.0)));
        }
        if(box)
        {
          assert(close(volume,parent.measure()));
          for(int i=0;i<3;++i) assert(close(moment[i],parent.measure()*center[i]));
        }
      }
    }
  }

  double g[3][3]={{2,0,0},{0,3,0},{0,0,4}};
  assert(validateFileGeometry(0,24,g)==InputError::None);
  assert(validateFileGeometry(1,12,g)==InputError::None);
  assert(validateFileGeometry(2,4,g)==InputError::None);
  assert(validateFileGeometry(2.5,4,g)==InputError::Shape);
  assert(validateFileGeometry(5,4,g)==InputError::Shape);
  assert(validateFileGeometry(3,0.6,g,0.1)==InputError::None);
  assert(validateFileGeometry(4,0.3,g,0.1)==InputError::None);
  assert(validateFileGeometry(2,0,g)==InputError::Volume);
  assert(validateFileGeometry(2,5,g)==InputError::Measure);
  g[0][2]=0.1;
  assert(validateFileGeometry(1,12,g)==InputError::Wedge);
  g[0][2]=0; g[0][0]=-2;
  assert(validateFileGeometry(0,24,g)==InputError::Generators);
  g[0][0]=std::numeric_limits<double>::quiet_NaN();
  assert(validateFileGeometry(0,24,g)==InputError::Generators);
  // Backend connectivity must map the same canonical vertices to each format.
  int fromSiloWedge[6]={2,1,5,3,0,4},fromSiloTet[4]={1,0,2,3};
  for(int i=0;i<6;++i) assert(particleGeometry::siloVertex(1,fromSiloWedge[i])==particleGeometry::vtkVertex(1,i));
  for(int i=0;i<4;++i) assert(particleGeometry::siloVertex(2,fromSiloTet[i])==i);
  std::cout << "Unified geometry: 5000 domains; " << childrenChecked
            << " box children; affine, centroid, volume, diameter and input checks passed.\n";
}
