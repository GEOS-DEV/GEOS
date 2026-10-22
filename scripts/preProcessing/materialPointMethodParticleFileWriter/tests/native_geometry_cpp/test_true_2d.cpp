#include "physicsSolvers/solidMechanics/MPMMappingKernels.hpp"
#include "physicsSolvers/solidMechanics/MPMGeometryMeasures.hpp"
#include "mesh/ParticleSubdivision.hpp"
#include <cassert>
#include <iostream>
#include <limits>
#include <map>
using namespace geos;

static void near(double a,double b,double tol=2e-12) { assert(std::abs(a-b)<tol); }

template<int DIM, bool SIMPLEX, bool CUBIC=false>
void mapping(bool const pointOnly=false)
{
  int grid[9][9][9]={}; double coordinates[729][3]={}; int next=0;
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<(DIM==2?1:9);++k)
  { grid[i][j][k]=next; coordinates[next][0]=i; coordinates[next][1]=j; coordinates[next++][2]=k; }
  double const origin[3]={}, h[3]={1,1,1}; int cells[3]={8,8,DIM==2?0:8};
  double center[3]={3.1,3.3,DIM==2?0:3.2};
  double generators[DIM][DIM]={};
  for(int d=0;d<DIM;++d) generators[d][d]=0.6+0.1*d;
  generators[0][1]=0.17; generators[1][0]=-0.13;
  constexpr int count=SIMPLEX ? (DIM+1)*(1<<DIM) : (1<<(2*DIM));
  int nodes[count]; double values[count], gradients[count][3];
  bool ok;
  int used=count;
  if(pointOnly)
  {
    ok=mpmMapping::point<DIM,CUBIC>(center,grid,origin,h,cells,nodes,values,gradients);
    used=CUBIC ? (1<<(2*DIM)) : 1<<DIM;
  }
  else ok=mpmMapping::domain<DIM,SIMPLEX>(center,generators,grid,origin,h,cells,nodes,values,gradients);
  assert(ok);
  double sum=0, x[3]={}, gradient[3]={}, linearGradient[3][3]={};
  for(int n=0;n<used;++n)
  {
    assert(nodes[n]>=0 && nodes[n]<next); sum+=values[n];
    for(int d=0;d<3;++d)
    {
      x[d]+=values[n]*coordinates[nodes[n]][d]; gradient[d]+=gradients[n][d];
      for(int k=0;k<DIM;++k) linearGradient[d][k]+=coordinates[nodes[n]][d]*gradients[n][k];
    }
  }
  near(sum,1);
  for(int d=0;d<3;++d)
  {
    near(x[d],center[d]); near(gradient[d],0);
    for(int k=0;k<DIM;++k) near(linearGradient[d][k],d==k?1:0);
  }
  center[0]=std::numeric_limits<double>::quiet_NaN();
  assert((!mpmMapping::point<DIM,false>(center,grid,origin,h,cells,nodes,values,gradients)));
  center[0]=-1e-15; center[1]=4;
  assert((mpmMapping::point<DIM,false>(center,grid,origin,h,cells,nodes,values,gradients)));
  center[0]=-1e-6;
  assert((!mpmMapping::point<DIM,false>(center,grid,origin,h,cells,nodes,values,gradients)));
}

int main()
{
  mapping<2,false>(); mapping<2,true>(); mapping<2,false>(true); mapping<2,false,true>(true);
  mapping<3,false>(); mapping<3,true>(); mapping<3,false>(true); mapping<3,false,true>(true);
  double g[2][2]={{2,0.3},{0.2,1.4}}, original[2][2]={{2,0.3},{0.2,1.4}};
  double f[3][3]={{1.1,0.4,0},{-0.2,1,0},{0,0,1}};
  auto quad=particleGeometry::affine::Domain::fromGenerators(particleGeometry::Shape::Parallelogram,g);
  near(quad.measure(),2.74);
  assert(particleGeometry::affine::validateFileGeometry(3,2.74*0.07,g,0.07)==particleGeometry::affine::InputError::None);
  assert(particleGeometry::affine::validateFileGeometry(3,2.74,g,0.07)==particleGeometry::affine::InputError::Measure);
  particleGeometry::affine::advance<2>(g,f);
  near(g[0][0],2.32); near(g[0][1],-0.1);
  particleGeometry::affine::reset<2>(original,g,f,false); near(g[0][0],2.32);
  particleGeometry::affine::reset<2>(original,g,f,true); near(g[0][0],2.32);
  for(double scale: {1e-15,1.0,1e15})
  {
    double g0[2][2]={{2*scale,.2*scale},{-.1*scale,scale}}, current[2][2], recovered[3][3];
    particleGeometry::affine::deform<2>(g0,f,current);
    assert(particleGeometry::affine::recoverDeformationGradient2D(g0,current,recovered));
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) near(recovered[i][j],f[i][j]);
  }
  for(auto shape: {particleGeometry::Shape::Parallelogram,particleGeometry::Shape::Triangle})
  {
    auto parent=particleGeometry::affine::Domain::fromGenerators(shape,g);
    auto plan=particleGeometry::subdivision::sizePlan(parent,2,1.1);
    double sum=0, dx[3]={};
    assert(plan.count<=4);
    for(int c=0;c<plan.count;++c)
    {
      auto child=plan.children[c].geometry(parent);
      assert(child.signedMeasure()>0); sum+=child.measure();
      double offset[3];plan.children[c].offset(parent,offset);
      for(int d=0;d<3;++d) dx[d]+=plan.children[c].fraction*offset[d];
    }
    near(sum,parent.measure()); for(double d:dx) near(d,0);
  }
  double n[3]={1,0,0}, p[3]={};
  near(mpmMeasures::lineArea(2,3,0.07,n,p),3*0.07);
  p[0]=1; near(mpmMeasures::lineArea(2,3,0.07,n,p),1.5*0.07);
  p[0]=2; near(mpmMeasures::lineArea(2,3,0.07,n,p),0);
  p[0]=0;n[0]=n[1]=1;
  near(mpmMeasures::lineArea(1,1,0.07,n,p),2*std::sqrt(2.0)/3*0.07);
  assert(mpmMapping::dispatch(2,[](auto d){return decltype(d)::value;})==2);
  assert(mpmMapping::dispatch(3,[](auto d){return decltype(d)::value;})==3);
  std::cout<<"DIM=2/3 Q1, B-spline, CPDI/CPTI, compact geometry, splitting, thickness and line-area invariants passed\n";
}
