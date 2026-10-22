"""Compile production background-pressure and 2D force bodies with CPU adapters.

Checks full-ring boundary resultants and cancellation at bonded CPDI/CPTI
interfaces. This standalone test does not build GEOS or test restart I/O.
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
  // Existing Cartesian box and sphere modes retain their selection behavior.
  pressure.m_uniformConfiningPressure=0;
  bool rejected=false;try {pressure.updateGridBackgroundStress(nodes);}catch(std::runtime_error const&){rejected=true;}
  assert(rejected);pressure.m_axisymmetric=false;
  for(auto& row:nodes.stress) for(double& value:row)value=0;
  pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[3][3][0]][0],0));assert(near(nodes.stress[map[1][1][0]][0],-7));
  pressure.m_confiningPressureSphereRadius=2;
  for(auto& row:nodes.stress) for(double& value:row)value=0;
  pressure.updateGridBackgroundStress(nodes);
  assert(near(nodes.stress[map[1][1][0]][0],0));assert(near(nodes.stress[map[4][4][0]][0],-7));
  std::cout<<"Production uniform pressure, bonded CPDI/CPTI cancellation, cylindrical resultants, pressure work and Cartesian selections passed.\n";
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
