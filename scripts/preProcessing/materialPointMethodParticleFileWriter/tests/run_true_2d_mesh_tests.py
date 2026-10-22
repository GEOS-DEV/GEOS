"""Execute production native connectivity/area and partition routines with adapters.

The Cartesian-rank adapter checks coordinate dimensionality; it is not MPI.
"""
from pathlib import Path
import re
import subprocess
import tempfile
root=next(p for p in Path(__file__).resolve().parents if (p/'src/coreComponents/mesh').exists())
mesh=root/'src/coreComponents/mesh'
def function(s,signature):
    a=s.index(signature);b=s.index('{',a)+1;depth=1
    while depth:
        depth+=(s[b]=='{')-(s[b]=='}');b+=1
    return s[a:b]
internal=(mesh/'generators/InternalMeshGenerator.cpp').read_text()
relation=function(internal,'static void getElemToNodesRelationInBox(')
native_box=function(internal,'if( elementType == ElementType::Quadrilateral || elementType == ElementType::Triangle )')
util=(mesh/'generators/CellBlockUtilities.cpp').read_text()
a=util.index('    case ElementType::Quadrilateral:',util.index('localIndex getFaceNodes('));b=util.index('    case ElementType::Hexahedron:',a)
face=util[a:b]
source=(mesh/'CellElementSubRegion.cpp').read_text()
a=source.index('    case ElementType::Quadrilateral:',source.index('calculateElementCenterAndVolume'));b=source.index('    case ElementType::Hexahedron:',a)
area=source[a:b]
partition=(mesh/'mpiCommunications/SpatialPartition.cpp').read_text()
methods='\n'.join(function(partition,s).replace('SpatialPartition::','') for s in (
    'void SpatialPartition::addNeighbors(', 'bool SpatialPartition::isCoordInPartition(',
    'void SpatialPartition::setSpatialDimension(', 'void SpatialPartition::setPlaneStrainThickness('))
program=r'''
#include <array>
#include <vector>
#include <set>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iostream>
using real64=double;using integer=int;using localIndex=int;using MPI_Comm=int;
#define GEOS_ERROR_IF(c,...) do{if(c)throw std::runtime_error("invalid input");}while(false)
#define GEOS_ERROR(...) throw std::runtime_error("invalid topology")
namespace LvArray{template<class T,class U>T integerConversion(U x){return static_cast<T>(x);}}
enum class ElementType {Hexahedron,Wedge,Quadrilateral,Triangle,Tetrahedron};
template<class T> struct Vec:std::array<T,3>{T &operator()(int i){return (*this)[i];}T const &operator()(int i)const{return (*this)[i];}};
struct NeighborCommunicator{int rank;explicit NeighborCommunicator(int r):rank(r){}};
namespace MpiWrapper {int dimension=3;int cartRank(MPI_Comm,int const *coords){int r=0;for(int d=0;d<dimension;++d){assert(coords[d]>=0 && coords[d]<3);r=3*r+coords[d];}return r;}}
double MapValueToRange(double x,double lo,double hi){double y=std::fmod(x-lo,hi-lo);return lo+(y<0?y+hi-lo:y);}
struct Partition{
 int m_spatialDimension=3;double m_planeStrainThickness=1;
 Vec<int> m_partitions,m_coords,m_periodic;Vec<double> m_min,m_max,m_gridMin,m_gridMax;
 std::vector<NeighborCommunicator> m_neighbors;
'''+methods+'\n};\n'+relation+r'''
void cellNodes(int firstNodeIndex,int const (&numNodesInDir)[3],int (&nodeOfBox)[8],ElementType elementType)
{
'''+native_box+r'''
}
int faceNodes(ElementType elementType,int elemIdx,int faceNumber,int const (*elementToNodes)[8],std::array<int,2>& faceNodes)
{switch(elementType){
'''+face+r'''
 default:throw std::runtime_error("not native");}}
double cellArea(ElementType m_elementType,double const (&input)[4][3])
{
 double m_elementVolume[1]={};int k=0;
 auto getElementCoordinatesaAndComputeElementCenter=[&](auto& Xlocal){for(int a=0;a<4;++a)for(int d=0;d<3;++d)Xlocal[a][d]=input[a][d];};
 switch(m_elementType){
'''+area+r'''
 default:throw std::runtime_error("not native");}
 return m_elementVolume[0];
}
int main()
{
 // Independent Euler/count/area checks on a non-square 7x5 structured mesh.
 int const nx=7,ny=5,nodes[3]={nx+1,ny+1,1};std::set<std::pair<int,int>> edges;int count=0;
 for(int j=0;j<ny;++j)for(int i=0;i<nx;++i)
 {
   int box[8]={},index[3]={i,j,0},order[8]={},connect[1][8]={};
   cellNodes(i+(nx+1)*j,nodes,box,ElementType::Quadrilateral);
   getElemToNodesRelationInBox(ElementType::Quadrilateral,0,index,0,order,4);
   double positions[4][3]={};
   for(int a=0;a<4;++a)
   {
     int id=box[order[a]];assert(id>=0 && id<(nx+1)*(ny+1));connect[0][a]=id;
     positions[a][0]=1e6+.2*(id%(nx+1));positions[a][1]=-2e6+.3*(id/(nx+1));
   }
   assert(std::abs(cellArea(ElementType::Quadrilateral,positions)-.06)<1e-9);
   for(int f=0;f<4;++f){std::array<int,2> e;assert(faceNodes(ElementType::Quadrilateral,0,f,connect,e)==2);if(e[0]>e[1])std::swap(e[0],e[1]);edges.emplace(e[0],e[1]);}
   double sum=0;
   for(int t=0;t<2;++t)
   {
     getElemToNodesRelationInBox(ElementType::Triangle,0,index,t,order,3);
     for(int a=0;a<3;++a){int id=box[order[a]];positions[a][0]=.2*(id%(nx+1));positions[a][1]=.3*(id/(nx+1));connect[0][a]=id;}
     sum+=cellArea(ElementType::Triangle,positions);
     for(int f=0;f<3;++f){std::array<int,2> e;assert(faceNodes(ElementType::Triangle,0,f,connect,e)==2);assert(e[0]!=e[1]);}
   }
   assert(std::abs(sum-.06)<1e-12);++count;
 }
 assert(count==nx*ny);assert(static_cast<int>(edges.size())==nx*(ny+1)+ny*(nx+1));
 for(int dim:{2,3})for(int periodic:{0,1})for(int center:{0,1})
 {
   Partition p;for(int d=0;d<3;++d){p.m_partitions[d]=3;p.m_coords[d]=center;p.m_periodic[d]=periodic;p.m_min[d]=1;p.m_max[d]=2;p.m_gridMin[d]=0;p.m_gridMax[d]=3;}
   if(dim==2)p.m_partitions[2]=1;
   p.setSpatialDimension(dim);MpiWrapper::dimension=dim;int c[3]={},comm=0;p.addNeighbors(0,comm,c);
   int expected=(periodic||center)?(dim==2?8:26):(dim==2?3:7);
   assert(static_cast<int>(p.m_neighbors.size())==expected);
   assert(p.isCoordInPartition(1.5,0));assert(!p.isCoordInPartition(2.1,0));
   if(dim==2){assert(c[2]==0);assert(p.isCoordInPartition(1e100,2));}
   p.setPlaneStrainThickness(.03);assert(std::abs(p.m_planeStrainThickness-.03)<1e-15);
   bool bad=false;try{p.setPlaneStrainThickness(0);}catch(std::runtime_error const&){bad=true;}assert(bad);
 }
 std::cout<<"Production mesh fragments: quad/tri connectivity, edge uniqueness, positive areas, 2D/3D neighbor topology and XY ownership passed.\n";
}
'''
with tempfile.TemporaryDirectory(prefix='geos-native-mesh-') as folder:
    folder=Path(folder);cpp=folder/'test.cpp';cpp.write_text(program);exe=folder/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer',str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
