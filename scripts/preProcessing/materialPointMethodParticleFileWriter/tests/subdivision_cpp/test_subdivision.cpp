/* Standalone tests of the production geometry/subdivision headers. */
#include "mesh/ParticleSubdivision.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
using namespace geos;
using namespace particleGeometry;
using namespace particleGeometry::subdivision;
using affine::Domain;

void close( double a, double b, double scale = 1.0 )
{ assert( std::abs(a-b) < 2e-10 * std::max(scale,std::max(std::abs(a),std::abs(b))) ); }
void inverse( double const (&a)[3][3], double (&b)[3][3] )
{
  double det=affine::determinant(a);
  for(int i=0;i<3;++i) for(int j=0;j<3;++j)
    b[i][j]=(a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3])/det;
}
bool inside( Domain const & d, double const (&point)[3], double tolerance )
{
  double inv[3][3],q[3]={}; inverse(d.generators,inv);
  for(int a=0;a<3;++a) for(int i=0;i<3;++i) q[a]+=inv[i][a]*point[i];
  if(d.shape==Shape::Parallelepiped || d.shape==Shape::Parallelogram)
  { for(int a=0;a<affine::dimension(d.shape);++a) if(std::abs(q[a])>0.5+tolerance) return false; }
  else
  {
    int n=d.shape==Shape::Tetrahedron?3:2; double sum=0;
    for(int a=0;a<n;++a) { if(q[a]<-1.0/(n+1)-tolerance) return false; sum+=q[a]; }
    if(sum>1.0/(n+1)+tolerance) return false;
    if(d.shape==Shape::TriangleWedge && std::abs(q[2])>0.5+tolerance) return false;
  }
  return true;
}
std::mt19937 rng(271828);
double uniform() { return std::uniform_real_distribution<double>(0.00001,0.99999)(rng); }
Domain randomDomain( Shape shape, int dimensions )
{
  Domain d;d.shape=shape;
  double a=6.28*uniform(),b=dimensions==3?6.28*uniform():0;
  double rotation[3][3]={{std::cos(a),-std::sin(a)*std::cos(b),std::sin(a)*std::sin(b)},
                         {std::sin(a),std::cos(a)*std::cos(b),-std::cos(a)*std::sin(b)},
                         {0,std::sin(b),std::cos(b)}};
  double g[3][3]={{0.2+4*uniform(),2*uniform()-1,dimensions==3?uniform():0},
                 {0,0.2+3*uniform(),dimensions==3?uniform():0},{0,0,0.2+2*uniform()}};
  affine::deform(g,rotation,d.generators);
  return d;
}
void sample( Domain const & parent, double (&x)[3] )
{
  double q[3]={};
  if(parent.shape==Shape::Parallelepiped || parent.shape==Shape::Parallelogram)
    for(int a=0;a<affine::dimension(parent.shape);++a) q[a]=uniform()-0.5;
  else
  {
    int n=parent.shape==Shape::Tetrahedron?3:2;
    double bary[4]={},sum=0;
    for(int v=0;v<=n;++v) {bary[v]=-std::log(uniform());sum+=bary[v];}
    for(int a=0;a<n;++a) q[a]=bary[a+1]/sum-1.0/(n+1);
    if(parent.shape==Shape::TriangleWedge) q[2]=uniform()-0.5;
  }
  for(int i=0;i<3;++i) { x[i]=0;for(int a=0;a<3;++a) x[i]+=q[a]*parent.generators[a][i]; }
}
void checkPlan( Domain const & parent, Domain const & reference, Plan const & plan )
{
  double total=0,volume=0,referenceVolume=0,centroid[3]={};
  double increment[3][3]={{1.1,0.13,0},{0.02,0.9,0},{0,0,1}};
  for(int c=0;c<plan.count;++c)
  {
    Child const & map=plan.children[c];Domain d=map.geometry(parent),d0=map.geometry(reference);
    double offset[3];map.offset(parent,offset);
    assert(d.signedMeasure()>0 && d0.signedMeasure()>0);
    close(d.measure()/parent.measure(),map.fraction);
    close(d0.measure()/reference.measure(),map.fraction);
    total+=map.fraction;volume+=d.measure();referenceVolume+=d0.measure();
    for(int i=0;i<3;++i) centroid[i]+=map.fraction*offset[i];
    for(int v=0;v<d.verticesCount();++v)
    {
      double x[3];d.vertexOffset(v,x);for(int i=0;i<3;++i)x[i]+=offset[i];
      assert(inside(parent,x,2e-10));
    }
    // Split commutes with physical convection and FAndGeometry reset.
    Domain updatedParent=parent;affine::advance(updatedParent.generators,increment);
    Domain updatedChild=d;affine::advance(updatedChild.generators,increment);
    Domain expected=map.geometry(updatedParent);
    Domain resetParent=parent;affine::reset(reference.generators,resetParent.generators,increment,true);
    Domain resetChild=d;affine::reset(d0.generators,resetChild.generators,increment,true);
    Domain expectedReset=map.geometry(resetParent);
    Domain kept=d;affine::reset(d0.generators,kept.generators,increment,false);
    for(int a=0;a<3;++a)for(int i=0;i<3;++i)
    {close(updatedChild.generators[a][i],expected.generators[a][i]);close(resetChild.generators[a][i],expectedReset.generators[a][i]);close(kept.generators[a][i],d.generators[a][i]);}
  }
  close(total,1);close(volume,parent.measure());close(referenceVolume,reference.measure());
  for(double x:centroid)close(x,0);
  // Independent point ownership check: no overlap/gap inside the parent.
  for(int trial=0;trial<80;++trial)
  {
    double x[3];sample(parent,x);int owners=0;
    for(int c=0;c<plan.count;++c)
    {
      double dx[3],local[3];plan.children[c].offset(parent,dx);
      for(int i=0;i<3;++i)local[i]=x[i]-dx[i];
      if(inside(plan.children[c].geometry(parent),local,1e-12))++owners;
    }
    assert(owners==1);
  }
}
int main()
{
  int count=0;
  for(Shape shape:{Shape::Parallelepiped,Shape::Tetrahedron,Shape::TriangleWedge,Shape::Parallelogram,Shape::Triangle})
    for(int dimensions:{2,3})
    {
      if((affine::dimension(shape)==2 && dimensions!=2)||(shape==Shape::Tetrahedron && dimensions==2))continue;
      for(int iteration=0;iteration<300;++iteration)
      {
        Domain p=randomDomain(shape,dimensions),r=randomDomain(shape,dimensions);
        Plan plan=sizePlan(p,dimensions,diameter(p,dimensions)*(0.2+1.2*uniform()));
        assert(plan.count>=1 && plan.count<=(1<<dimensions));
        checkPlan(p,r,plan);++count;
        if(dimensions==2 && shape==Shape::TriangleWedge)
          for(int c=0;c<plan.count;++c)
          {auto child=plan.children[c].geometry(p);for(int i=0;i<3;++i)close(child.generators[2][i],p.generators[2][i]);}
      }
    }
  // Every canonical simplex edge, not only edges incident on vertex zero.
  for(int n:{2,3})for(int a=0;a<=n;++a)for(int b=a+1;b<=n;++b)
  {
    Plan plan;plan.count=2;for(int side=0;side<2;++side)plan.children[side]=simplexChild(n,a,b,side);
    Domain p=randomDomain(n==3?Shape::Tetrahedron:Shape::TriangleWedge,n);
    checkPlan(p,p,plan);
  }
  // Existing box masks and exact legacy child ordering.
  for(int dimensions:{2,3})for(int mask=1;mask<(1<<dimensions);++mask)
  {
    Domain p=randomDomain(Shape::Parallelepiped,dimensions);Plan plan=boxPlan(mask,dimensions);checkPlan(p,p,plan);
    for(int c=0;c<plan.count;++c)
    {
      int pos=0;for(int a=0;a<dimensions;++a)if(plan.children[c].center[a]>0)pos|=1<<a;
      Domain expected;double dx[3],actual[3];assert(p.boxChild(mask,pos,expected,dx));
      Domain actualDomain=plan.children[c].geometry(p);plan.children[c].offset(p,actual);
      for(int a=0;a<3;++a){close(dx[a],actual[a]);for(int i=0;i<3;++i)close(expected.generators[a][i],actualDomain.generators[a][i]);}
    }
  }
  // Repeated steps reach a diameter target without assuming each bisection does.
  for(Shape shape:{Shape::Tetrahedron,Shape::TriangleWedge})
  {
    Domain p=randomDomain(shape,3);double target=0.23*diameter(p,3);std::vector<Domain> leaves{p};
    bool done=false;
    for(int step=0;step<8 && !done;++step)
    {
      std::vector<Domain> next;done=true;
      for(Domain const & leaf:leaves){Plan plan=sizePlan(leaf,3,target);if(plan.unresolved)done=false;for(int c=0;c<plan.count;++c)next.push_back(plan.children[c].geometry(leaf));}
      leaves=next;
    }
    assert(done);for(Domain const & leaf:leaves)assert(diameter(leaf,3)<=target*(1+1e-10));
  }
  std::cout<<count<<" randomized geometry plans; partition, orientation, centroid, reset, box masks, and repeated bisection passed.\n";
}
