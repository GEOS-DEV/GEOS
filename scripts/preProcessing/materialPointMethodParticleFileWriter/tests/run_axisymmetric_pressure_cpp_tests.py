"""Compile production background-pressure and 2D force bodies with CPU adapters.

Checks spatial selectors, simultaneous inner/outer full-ring pressure loads,
and cancellation at bonded CPDI/CPTI interfaces. Does not build GEOS or test
restart I/O.
"""
from pathlib import Path
import subprocess
import tempfile
from axisymmetric_test_support import program_prefix

root, prefix = program_prefix()
source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
start = source.index('void SolidMechanicsMPM::updateGridBackgroundStress(')
stop = source.index('\n/**', start)
background = source[start:stop].replace('SolidMechanicsMPM::', '')
start = source.index('// 2D path. This preserves the original x-y component behavior')
stop = source.index('          RAJA::atomicAdd(', start)
force = source[start:stop]
# The mapping adapters do not use arrayView2d; the pressure body uses six-wide
# stress arrays. Both stress and coordinate storage below use padded rows.
prefix = prefix.replace('using arrayView2d=T(*)[3];', 'using arrayView2d=T(*)[6];')
adapters = r'''
#include <array>
#include <stdexcept>
namespace LvArray { namespace math { inline double square(double x) { return x*x; } }
namespace tensorOps {
template<int N,class A,class B> void copy(A&& a,B const& b)
{ for(int i=0;i<N;++i) a[i]=b[i]; }
}}
struct serialPolicy {};
template<class Policy,class F> void forAll(int n,F f) { for(int i=0;i<n;++i) f(i); }
template<class T> using array2d=T[729][6];
struct NodeManager {
  double stress[729][6]={},xyz[729][6]={};
  int size() const { return 729; }
  template<class T> T& getReference(char const*) {return stress;}
  double const (*referencePosition() const)[6] {return xyz;}
};
namespace viewKeyStruct {char const* gridBackgroundStressString() {return "stress";} }
#define GEOS_ERROR_IF(condition,message) do {if(condition) throw std::runtime_error(message);} while(false)
struct PressureKernel {
  bool m_axisymmetric=true;
  int m_uniformConfiningPressure=1,m_enableBoreholePressure=0,m_enableConfiningPressure=1,m_numDims=2;
  double m_boreholeRadius=1,m_boreholeStress[6]={};
  double m_confiningPressureSphereRadius=0,m_confiningStress[6]={-7,-7,-7,0,0,0};
  double m_confiningPressureBoxMin[3]={1,1,-1},m_confiningPressureBoxMax[3]={3,3,1};
'''+background+r'''
};
std::array<double,2> mappedForce(NodeManager const& nodes,int mappedNode,double pVolume,
                                 double const* grad,double stress) {
  auto const& gridBackgroundStress=nodes.stress;
  bool const m_axisymmetric=true;
  double const pStress0=stress,pStress1=stress,pStress2=stress,pStress5=0,supplementalPressure=0;
  double const grad0=grad[0],grad1=grad[1],grad2=grad[2];
'''+force+r'''
  return {internalForce0,internalForce1};
}
bool near(double a,double b) {return std::abs(a-b)<2.e-10*std::max({1.,std::abs(a),std::abs(b)});}
void clearStress(NodeManager& nodes) {
  // Production grid initialization does this before every background update.
  for(auto& row:nodes.stress) for(double& value:row)value=0;
}
int main() {
  NodeManager nodes;
  int map[9][9][9],cells[3]={8,8,0};
  double origin[3]={-1,-1,0},h[3]={1,1,1};
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<9;++k) {
    int n=(i*9+j)*9+k; map[i][j][k]=n;
    nodes.xyz[n][0]=i-1; nodes.xyz[n][1]=j-1;
  }
  PressureKernel pressure;
  pressure.updateGridBackgroundStress(nodes);
  for(auto const& row:nodes.stress) for(int d=0;d<6;++d) assert(near(row[d],d<3?-7.:0.));
  // Tile a radius-4, height-4 cylinder, with a CPDI core and CPTI outer layer.
  // At zero material stress, pressure must act only on the external boundary.
  // At sigma=-p I, every particle-to-node force cancels independently.
  Kernel mapping; mapping.m_axisymmetric=true;
  double forces[9][9][2]={},work=0,volume=0;
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) for(int tri=0;tri<(i==3?2:1);++tri) {
    double center[3]={i+.5,j+.5,0},g[3][3]={{1,0,0},{0,1,0},{0,0,1}},area=1;
    ParticleType type=ParticleType::CPDI;
    if(i==3) {
      type=ParticleType::CPTI;area=.5;
      if(tri==0) {center[0]=i+2./3;center[1]=j+1./3;g[1][0]=1;}
      else {center[0]=i+1./3;center[1]=j+2./3;g[0][1]=1;}
    }
    double V=axisymmetric::volume(area,center[0]);volume+=V;
    int mapped[16];double w[16],grad[16][3];
    assert(mapping.computeParticleMapping<2>(type,center,g,map,origin,h,cells,mapped,w,grad));
    for(int a=0;a<(type==ParticleType::CPTI?12:16);++a) {
      auto f=mappedForce(nodes,mapped[a],V,grad[a],0);
      auto balanced=mappedForce(nodes,mapped[a],V,grad[a],-7);
      assert(near(balanced[0],0)&&near(balanced[1],0));
      int ni=mapped[a]/81,nj=(mapped[a]/9)%9;
      forces[ni][nj][0]+=f[0];forces[ni][nj][1]+=f[1];
      // v_r=r, v_y=y gives divergence=3; pressure work is -3*p*V.
      work+=f[0]*nodes.xyz[mapped[a]][0]+f[1]*nodes.xyz[mapped[a]][1];
    }
  }
  double outer=0,top=0,bottom=0;
  for(int j=0;j<9;++j) outer+=forces[5][j][0];
  for(int i=0;i<9;++i) {bottom+=forces[i][1][1];top+=forces[i][5][1];}
  for(int i=1;i<5;++i) for(int j=1;j<=5;++j) assert(near(forces[i][j][0],0));
  for(int i=1;i<=5;++i) for(int j=2;j<5;++j) assert(near(forces[i][j][1],0));
  assert(near(outer,-7*axisymmetric::twoPi*4*4));
  assert(near(top,-7*.5*axisymmetric::twoPi*16));assert(near(bottom,-top));
  assert(near(work,-3*7*volume));
  // Independent inner/outer pressures on an annulus, CPDI bulk / CPTI skin.
  pressure.m_uniformConfiningPressure=0;
  pressure.m_enableBoreholePressure=1;pressure.m_boreholeRadius=1.25;
  for(int d=0;d<3;++d)pressure.m_boreholeStress[d]=-11;
  pressure.m_confiningPressureBoxMin[0]=-99;pressure.m_confiningPressureBoxMin[1]=-99;
  pressure.m_confiningPressureBoxMax[0]=3.75;pressure.m_confiningPressureBoxMax[1]=99;
  // z bounds intentionally exclude every node: axisymmetric boxes ignore z.
  pressure.m_confiningPressureBoxMin[2]=8;pressure.m_confiningPressureBoxMax[2]=9;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  for(auto const& row:nodes.stress) assert(near(row[0],row[2])); // hoop pressure retained
  double annulusForces[9][9][2]={},radialWork=0;
  for(int i=1;i<4;++i) for(int j=0;j<4;++j) for(int tri=0;tri<(i==3?2:1);++tri) {
    double center[3]={i+.5,j+.5,0},g[3][3]={{1,0,0},{0,1,0},{0,0,1}},area=1;
    ParticleType type=ParticleType::CPDI;
    if(i==3) {
      type=ParticleType::CPTI;area=.5;
      if(tri==0) {center[0]=i+2./3;center[1]=j+1./3;g[1][0]=1;}
      else {center[0]=i+1./3;center[1]=j+2./3;g[0][1]=1;}
    }
    double V=axisymmetric::volume(area,center[0]);
    int mapped[16];double w[16],grad[16][3];
    assert(mapping.computeParticleMapping<2>(type,center,g,map,origin,h,cells,mapped,w,grad));
    for(int a=0;a<(type==ParticleType::CPTI?12:16);++a) {
      auto f=mappedForce(nodes,mapped[a],V,grad[a],0);
      int ni=mapped[a]/81,nj=(mapped[a]/9)%9;
      annulusForces[ni][nj][0]+=f[0];annulusForces[ni][nj][1]+=f[1];
      radialWork+=f[0]*nodes.xyz[mapped[a]][0];
    }
  }
  double inner=0;outer=0;
  for(int j=0;j<9;++j) {inner+=annulusForces[2][j][0];outer+=annulusForces[5][j][0];}
  assert(near(inner,11*axisymmetric::twoPi*1*4));
  assert(near(outer,-7*axisymmetric::twoPi*4*4));
  for(int i=3;i<5;++i) for(int j=0;j<9;++j)assert(near(annulusForces[i][j][0],0));
  assert(near(radialWork,axisymmetric::twoPi*4*(11-7*16)));

  // Borehole selector is independent of axial y and excludes its cutoff.
  pressure.m_enableConfiningPressure=0;pressure.m_boreholeRadius=2;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][8][0]][0],-11)); // r=1, y=7
  assert(near(nodes.stress[map[0][8][0]][0],-11)); // reflected buffer r=-1
  assert(near(nodes.stress[map[3][1][0]][0],0));   // r=2 exactly
  // Selection is recomputed from current coordinates, not cached grid indices.
  nodes.xyz[map[2][8][0]][0]=2.5;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][8][0]][0],0));nodes.xyz[map[2][8][0]][0]=1;
  // Finite axial box bounds can also load end surfaces. At a corner equality
  // is inside. Confinement overwrites the borehole in an overlap.
  pressure.m_enableConfiningPressure=1;
  pressure.m_confiningPressureBoxMin[0]=0;pressure.m_confiningPressureBoxMin[1]=0;
  pressure.m_confiningPressureBoxMax[0]=3;pressure.m_confiningPressureBoxMax[1]=3;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][8][0]][0],-7));
  assert(near(nodes.stress[map[2][2][0]][0],-11));
  assert(near(nodes.stress[map[4][4][0]][0],0));
  pressure.m_uniformConfiningPressure=1;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][2][0]][0],-7));
  // Origin-centered sphere in the meridional plane, ignoring geometric z.
  pressure.m_uniformConfiningPressure=0;pressure.m_enableBoreholePressure=0;
  pressure.m_confiningPressureSphereRadius=2;
  nodes.xyz[map[2][2][0]][2]=10;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][2][0]][0],0));
  assert(near(nodes.stress[map[3][1][0]][0],0));
  assert(near(nodes.stress[map[2][3][0]][0],-7));

  // Cartesian borehole, box and sphere selectors retain their behavior.
  pressure.m_axisymmetric=false;pressure.m_enableBoreholePressure=1;pressure.m_enableConfiningPressure=0;
  clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][8][0]][0],0));
  assert(near(nodes.stress[map[2][2][0]][0],-11));
  pressure.m_enableBoreholePressure=0;pressure.m_enableConfiningPressure=1;
  pressure.m_confiningPressureSphereRadius=0;
  pressure.m_confiningPressureBoxMin[0]=1;pressure.m_confiningPressureBoxMin[1]=1;
  pressure.m_confiningPressureBoxMin[2]=-1;pressure.m_confiningPressureBoxMax[2]=1;
  clearStress(nodes);
  pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[3][3][0]][0],0));assert(near(nodes.stress[map[1][1][0]][0],-7));
  pressure.m_confiningPressureSphereRadius=2;
  for(auto& row:nodes.stress) for(double& value:row)value=0;
  pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[1][1][0]][0],0));assert(near(nodes.stress[map[4][4][0]][0],-7));
  assert(near(nodes.stress[map[2][2][0]][0],0)); // Cartesian 2D sphere ignores z
  pressure.m_numDims=3;clearStress(nodes);pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[2][2][0]][0],-7)); // Cartesian 3D sphere includes z
  std::cout<<"Production pressure selectors, hoop balance, CPDI/CPTI inner/outer resultants, pressure work and Cartesian selections passed.\n";
}
'''
with tempfile.TemporaryDirectory(prefix='geos-axis-pressure-') as folder:
    tmp=Path(folder);(tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\nnamespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    cpp=tmp/'test.cpp';cpp.write_text(prefix+adapters);exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer',
                    '-I'+str(tmp),'-I'+str(root/'src/coreComponents'),str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
