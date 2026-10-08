#!/usr/bin/env python3
"""Compile the actual search/fit block extracted from SolidMechanicsMPM.cpp.

Minimal GEOS type/macro stubs isolate this block only. This is not a complete
GEOS adapter build. Tests modes 0, 1 and 2 through the production weighted fit.
"""
from pathlib import Path
import argparse, subprocess
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--cxx',default='g++')
parser.add_argument('--sanitize',action='store_true')
a=parser.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
sm=Path(__file__).resolve().parents[2]
source=(sm/'SolidMechanicsMPM.cpp').read_text()
start=source.index('    real64 const markerSearchRadius = surfaceFitRadius * cellDiagonal;')
stop=source.index('    if( timeMortar ) fittingSeconds += elapsed( fittingStart, mortarNow() );',start)
stop=source.index('\n',stop)
block=source[start:stop]
preamble=r'''
#include "MPMMortarSearch.hpp"
#include "MPMMortarIndexedSearch.hpp"
#include <chrono>
#include <cmath>
#include <set>
#include <memory>
#include <sstream>
#include <tuple>
#include <random>
#include <iostream>
using real64=double;using localIndex=int;using globalIndex=long long;
namespace mpm=geos::mpm;
namespace LvArray{namespace math{using std::min;using std::max;using std::floor;using std::sqrt;}}
#define GEOS_ERROR_IF(c,msg) do{if(c){std::ostringstream s;s<<msg;throw std::runtime_error(s.str());}}while(false)
struct WeakInterfaceMortarNodeKey{long long i,j,k;bool operator<(WeakInterfaceMortarNodeKey const&b)const{return std::tie(i,j,k)<std::tie(b.i,b.j,b.k);}};
struct WeakInterfaceTraceMarker{long long particleID;int group;std::array<double,3>point,normal;};
using Fit=std::array<double,7>;
std::vector<Fit> run(int mode,int numDims,int periodicMask)
{
 bool fastMarkerSearch=mode!=0,validateMarkerSearch=mode==2,timeMortar=true;
 int m_planeStrain=numDims==2,planeNodeK=1;
 double physicalMidZ=.5,cellDiagonal=std::sqrt(3.)*.1,surfaceFitRadius=2.5;
 int surfaceFitNeighbors=12;
 double surfaceFitRadiusSquared=surfaceFitRadius*surfaceFitRadius*cellDiagonal*cellDiagonal;
 double distanceWeightFloorSquared=.000625,minNormalCoherence=0,maxFitResidual=1.e6,explicitSurfaceNormalTolerance=1.e-12;
 long long markerSearchMismatches=0;
 double fittingSeconds=0,searchIndexSeconds=0,markerTreeNodes=0,markerTieReplays=0,markerCandidatesVisited=0,maximumMarkerCandidates=0,markerDistanceEvaluations=0,markerQueries=0;
 using Clock=std::chrono::steady_clock;
 auto mortarNow=[](){return Clock::now();};
 auto elapsed=[](Clock::time_point a,Clock::time_point b){return std::chrono::duration<double>(b-a).count();};
 std::array<double,3>gridGlobalMin{{-.1,-.1,-.1}},gridGlobalExtent{{1.2,1.2,1.2}};
 std::array<int,3>periodic{{periodicMask&1,(periodicMask>>1)&1,(periodicMask>>2)&1}};
 std::array<long long,3>globalCellCount{{12,12,12}};double hEl[3]={.1,.1,.1};int nEl[3]={12,12,12};
 auto positiveModulo=[](long long a,long long b){return (a%b+b)%b;};
 auto canonicalCoordinate=[&](double x,int d){if(!periodic[d])return x;double a=x-gridGlobalMin[d];a-=std::floor(a/gridGlobalExtent[d])*gridGlobalExtent[d];return gridGlobalMin[d]+a;};
 auto canonicalPoint=[&](std::array<double,3>x){for(int d=0;d<3;++d)x[d]=canonicalCoordinate(x[d],d);return x;};
 auto minimumImageDelta=[&](double x,int d){return periodic[d]?x-std::round(x/gridGlobalExtent[d])*gridGlobalExtent[d]:x;};
 auto cellKeyFromPoint=[&](std::array<double,3>const&x){return WeakInterfaceMortarNodeKey{static_cast<long long>(std::floor((x[0]+.1)/.1)),static_cast<long long>(std::floor((x[1]+.1)/.1)),static_cast<long long>(std::floor((x[2]+.1)/.1))};};
 std::vector<std::vector<std::vector<int>>>ijkMap(13,std::vector<std::vector<int>>(13,std::vector<int>(13)));
 std::vector<std::array<double,3>>gridPosition;
 for(int i=0;i<13;++i)for(int j=0;j<13;++j)for(int k=0;k<13;++k){ijkMap[i][j][k]=gridPosition.size();gridPosition.push_back({{-.1+.1*i,-.1+.1*j,-.1+.1*k}});}
 auto mortarGridPoint=[&](int i){return gridPosition[i];};
 std::mt19937_64 rng(782);std::uniform_real_distribution<double>u(-1,1);
 std::vector<WeakInterfaceTraceMarker>pairMarkers;
 for(int i=0;i<2200;++i){std::array<double,3>n{{u(rng),u(rng),u(rng)}};double l=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);std::array<double,3>p{};for(int d=0;d<3;++d){n[d]/=l;p[d]=.28*n[d];}pairMarkers.push_back({i,i%2,canonicalPoint(p),n});}
'''
post=r'''
 if(markerSearchMismatches)throw std::runtime_error("Production indexed search mismatch");
 std::vector<Fit> out;
 for(auto const & f:nodalFit)out.push_back({{f.phi,f.coherence,f.residual,double(f.near),double(f.hasScalar),double(f.qualityAccepted),double(f.valid)}});
 return out;
}
int main(){try{int checked=0;for(int dim:{2,3})for(int periodic:{0,1,3}){auto reference=run(0,dim,periodic);for(int mode:{1,2}){auto result=run(mode,dim,periodic);if(result!=reference)throw std::runtime_error("Production nodal fit changed");checked+=result.size();}}std::cout<<"PASS extracted production search/fit: exact nodal scalar/quality/availability comparisons="<<checked<<"\n";}catch(std::exception const&e){std::cerr<<e.what()<<"\n";return 1;}}
'''
cpp=a.output/'productionSearch.cpp';cpp.write_text(preamble+'\n'+block+'\n'+post)
flags=['-std=c++14','-O3']
if a.sanitize: flags=['-std=c++14','-O1','-g','-fsanitize=address,undefined','-fno-pie','-no-pie']
cmd=[a.cxx,*flags,'-I'+str(sm/'kernels'),str(cpp),'-o',str(a.output/'productionSearch')]
subprocess.run(cmd,check=True,timeout=60)
subprocess.run([str(a.output/'productionSearch')],check=True,timeout=90)
