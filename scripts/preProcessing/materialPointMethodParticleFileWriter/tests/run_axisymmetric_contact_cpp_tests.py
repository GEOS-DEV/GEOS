"""Compile production mapping, normal assembly and both contact laws on the CPU.

This extracts the actual source bodies; small adapters replace GEOS containers.
It does not exercise a complete GEOS timestep, MPI synchronization or a GPU.
"""
from pathlib import Path
import os
import re
import subprocess
import tempfile
from axisymmetric_test_support import program_prefix

root, prefix = program_prefix()
source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
enums = (root/'src/coreComponents/physicsSolvers/solidMechanics/MPMSolverEnums.hpp').read_text()


def extract(text, signature):
    start = text.index(signature)
    # Braces in commented-out code must not terminate an extracted function.
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                    lambda match: ' '*len(match.group()), text)
    end = masked.index('{', start) + 1
    depth = 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return text[start:end].replace('SolidMechanicsMPM::', '')


definitions = '\nnamespace mpm {\n'
for name in ('ContactGapCorrectionOption', 'OverlapCorrectionOption'):
    definitions += extract(enums, 'enum struct '+name) + ';\n'
definitions += '}\n'
start = source.index('          real64 surfaceNormalContribution0 =')
end = source.index('          // Optional explicit particle surface normal:', start)
normal = source[start:end]
adapters = r'''
#include <array>
#include <stdexcept>
#define GEOS_UNUSED_PARAM(x) x [[maybe_unused]]
#define GEOS_ERROR(message) throw std::runtime_error(message)
namespace LvArray { namespace tensorOps {
template<int N,class A,class B> void copy(A&& a,B const& b)
{for(int i=0;i<N;++i) a[i]=b[i];}
template<int N,class A,class B> void subtract(A&& a,B const& b)
{for(int i=0;i<N;++i) a[i]-=b[i];}
template<int N,class A,class B> void scaledCopy(A&& a,B const& b,double s)
{for(int i=0;i<N;++i) a[i]=s*b[i];}
template<int N,class A,class B> double AiBi(A const& a,B const& b)
{double s=0;for(int i=0;i<N;++i)s+=a[i]*b[i];return s;}
}}
'''
for name in ('computeOrthonormalBasis', 'computePairwiseNodalContactForce',
             'computePairwiseNodalContactImpulse'):
    adapters += extract(source, 'void SolidMechanicsMPM::'+name+'(') + '\n'

adapters += r'''
using Vec=std::array<double,3>;
Vec normalContribution(bool m_axisymmetric,double pVolume,double const* grad) {
  double const grad0=grad[0],grad1=grad[1],grad2=grad[2];
''' + normal + r'''
  // The production P2G scatter only writes components < numDims.
  return {surfaceNormalContribution0,surfaceNormalContribution1,
          m_axisymmetric ? 0.0 : surfaceNormalContribution2};
}
bool near(double a,double b,double tol=2.e-10) {
  return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});
}
void close(double a,double b) {
  if(!near(a,b)) {std::cerr<<"Mismatch: "<<a<<" vs "<<b<<'\n';std::abort();}
}
using Gap=mpm::ContactGapCorrectionOption;
Vec impulse(bool forcePath,Vec n,Vec va,Vec vb,double ma,double mb,double mu,
            Gap mode,double gap,bool separable=true) {
  double const dt=.025,h[3]={1,1,1};
  double unitNormal[3]={n[0],n[1],n[2]};
  Vec qa{},qb{},sa{},sb{},ja{},jb{},zero{};
  for(int d=0;d<3;++d) {
    qa[d]=ma*va[d];qb[d]=mb*vb[d];sa[d]=-.5*gap*n[d];sb[d]=.5*gap*n[d];
  }
  if(forcePath) {
    computePairwiseNodalContactForce(mode,mpm::OverlapCorrectionOption::Off,
      1.01,1.2,100,h,1,1,1.e-30,1,0,separable,dt,mu,unitNormal,ma,mb,ma,mb,
      va.data(),vb.data(),qa.data(),qb.data(),ma,mb,sa.data(),sb.data(),
      sa.data(),sb.data(),zero.data(),zero.data(),ja.data(),jb.data());
    for(int d=0;d<3;++d) {ja[d]*=dt;jb[d]*=dt;}
  } else {
    computePairwiseNodalContactImpulse(mode,mpm::OverlapCorrectionOption::Off,
      1.01,1.2,100,h,1,1,1.e-30,1,0,separable,dt,mu,unitNormal,ma,mb,ma,mb,
      va.data(),vb.data(),qa[0],qa[1],qa[2],qb[0],qb[1],qb[2],ma,mb,
      sa.data(),sb.data(),sa.data(),sb.data(),ja.data(),jb.data());
  }
  double before=0,after=0;
  for(int d=0;d<3;++d) {
    close(ja[d],-jb[d]);
    double aa=va[d]+ja[d]/ma,bb=vb[d]+jb[d]/mb;
    before+=ma*va[d]*va[d]+mb*vb[d]*vb[d];
    after+=ma*aa*aa+mb*bb*bb;
  }
  assert(after<=before+1.e-11*std::max(before,1.e-30));
  close(ja[2],0);
  return ja;
}

void contactLaws() {
  int cases=0;
  for(double angle : {0.,.23,.8,1.5707963267948966,2.4})
    for(double radius : {1.e-6,.01,1.,20.})
      for(double ratio : {.01,1.,100.}) for(double mu : {0.,.2,10.}) {
    Vec n={std::cos(angle),std::sin(angle),0},t={-n[1],n[0],0},va{},vb{};
    for(int d=0;d<2;++d) {va[d]=n[d]+.7*t[d];vb[d]=-n[d]-.7*t[d];}
    double ma=axisymmetric::twoPi*radius,mb=ratio*ma,reduced=ma*mb/(ma+mb);
    for(auto mode : {Gap::Simple,Gap::Implicit,Gap::Softened}) {
      auto j=impulse(true,n,va,vb,ma,mb,mu,mode,-.01);
      auto k=impulse(false,n,va,vb,ma,mb,mu,mode,-.01);
      double jn=0,jt=0;
      for(int d=0;d<3;++d) {close(j[d],k[d]);jn+=j[d]*n[d];jt+=j[d]*t[d];}
      close(jn,-2*reduced);close(jt,-std::min(mu*2*reduced,1.4*reduced));
      // Reversing the normal velocity makes the pair recede: no attraction.
      auto separating=impulse(true,n,vb,va,ma,mb,mu,mode,-.01);
      auto separatingFMPM=impulse(false,n,vb,va,ma,mb,mu,mode,-.01);
      for(int d=0;d<3;++d) {close(separating[d],0);close(separatingFMPM[d],0);}
      ++cases;
    }
    for(bool forcePath : {false,true}) {
      auto open=impulse(forcePath,n,va,vb,ma,mb,mu,Gap::Implicit,.25);
      for(double value:open) close(value,0);
      auto soft=impulse(forcePath,n,va,vb,ma,mb,mu,Gap::Softened,.25);
      close(soft[0]*n[0]+soft[1]*n[1],-1.5*reduced);
      auto bonded=impulse(forcePath,n,va,vb,ma,mb,mu,Gap::Simple,.25,false);
      for(int d=0;d<2;++d) close(va[d]+bonded[d]/ma,vb[d]-bonded[d]/mb);
    }
  }
  std::cout<<cases<<" force/impulse comparisons plus separation, gap, friction and energy checks passed.\n";
}

void mappedInterfaces() {
  int map[9][9][9],cells[3]={8,8,0};
  double origin[3]={-1,-1,0},h[3]={1,1,1};
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<9;++k) map[i][j][k]=(i*9+j)*9+k;
  Kernel mapper;mapper.m_axisymmetric=true;
  // Axial, cylindrical and inclined interfaces. All include particles at r=0;
  // the inclined case is tiled entirely with CPTI, the others mix CPDI/CPTI.
  for(int orientation=0;orientation<3;++orientation) {
    double normals[9][9][2][3]={},mass[9][9][2]={};
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
      bool const simplex=orientation==2 || i>=2;
      for(int tri=0;tri<(simplex?2:1);++tri) {
        double center[3]={i+.5,j+.5,0},g[3][3]={{1,0,0},{0,1,0},{0,0,1}},area=1;
        auto type=ParticleType::CPDI;
        if(simplex) {
          type=ParticleType::CPTI;area=.5;
          if(tri==0) {center[0]=i+2./3;center[1]=j+1./3;g[1][0]=1;}
          else {center[0]=i+1./3;center[1]=j+2./3;g[0][1]=1;}
        }
        int field=orientation==0 ? int(j>=2) : (orientation==1 ? int(i>=2) : int(center[1]>center[0]));
        int mapped[16];double values[16],grad[16][3];
        assert(mapper.computeParticleMapping<2>(type,center,g,map,origin,h,cells,mapped,values,grad));
        double volume=axisymmetric::volume(area,center[0]);
        for(int a=0;a<(simplex?12:16);++a) {
          int ni=mapped[a]/81,nj=(mapped[a]/9)%9;
          auto n=normalContribution(true,volume,grad[a]);
          mass[ni][nj][field]+=volume*values[a];
          for(int d=0;d<3;++d) normals[ni][nj][field][d]+=n[d];
        }
      }
    }
    // Internal interface contributions cancel in the union, including the axis.
    for(int i=1;i<5;++i) for(int j=2;j<5;++j) {
      assert(std::abs(normals[i][j][0][0]+normals[i][j][1][0])<1.e-12);
      assert(std::abs(normals[i][j][0][1]+normals[i][j][1][1])<1.e-12);
    }
    for(int k=2;k<5;++k) {
      int i=orientation==1?3:k,j=orientation==0?3:k;
      double expected[3]={orientation==0?0.:1.,orientation==0?1.:0.,0};
      if(orientation==2) {expected[0]=-std::sqrt(.5);expected[1]=std::sqrt(.5);}
      auto const& n=normals[i][j][0];double norm=std::hypot(n[0],n[1]);
      assert(norm>1);
      for(int d=0;d<2;++d) assert(std::abs(n[d]/norm-expected[d])<1.e-12);
      Vec direction={n[0]/norm,n[1]/norm,0},va=direction,vb={-direction[0],-direction[1],0};
      auto impulseA=impulse(true,direction,va,vb,mass[i][j][0],mass[i][j][1],0,Gap::Implicit,-.01);
      for(int d=0;d<2;++d) close(va[d]+impulseA[d]/mass[i][j][0],vb[d]-impulseA[d]/mass[i][j][1]);
    }
    // For axial/cylindrical interfaces this support sees only the axis,
    // whose radial face has zero revolved area. An inclined interface can
    // enter the same support at r>0 and legitimately contribute there.
    if(orientation<2) for(int j=1;j<=5;++j) for(int field=0;field<2;++field)
      assert(std::abs(normals[1][j][field][0])<1.e-12);
  }
  // Cartesian assembly is untouched, and the hoop slot is never a 2D normal.
  double gradient[3]={2,-3,7};
  auto cartesian=normalContribution(false,5,gradient),axis=normalContribution(true,5,gradient);
  close(cartesian[0],10);close(cartesian[1],-15);close(cartesian[2],35);
  close(axis[0],45);close(axis[1],-15);close(axis[2],0);
  std::cout<<"Mixed CPDI/CPTI axial, cylindrical and inclined interfaces: normal directions, interior cancellation and axis checks passed.\n";
}
int main() {mappedInterfaces();contactLaws();}
'''

with tempfile.TemporaryDirectory(prefix='geos-axis-contact-') as directory:
    tmp = Path(directory)
    (tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text(
        '#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\n'
        'namespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    cpp = tmp/'test.cpp'
    cpp.write_text(prefix+definitions+adapters)
    exe = tmp/'test'
    subprocess.run(['c++', '-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Werror',
                    '-Wshadow', '-Wfloat-equal', '-pedantic', '-fsanitize=address,undefined',
                    '-fno-omit-frame-pointer', '-I'+str(tmp), '-I'+str(root/'src/coreComponents'),
                    str(cpp), '-o', str(exe)], check=True)
    subprocess.run([str(exe)], check=True, env=os.environ)
