// Standalone tests of the production header; build with C++14 or newer.
#include "../kernels/MPMBoundaryMortar.hpp"
#include <iostream>
#include <random>
#include <sstream>

namespace bm = geos::mpm::boundaryMortar;
int checks = 0;
void check( bool c, char const * label ) { if( !c ) throw std::runtime_error( label ); ++checks; }
void near( double a, double b, double tol=2e-11 ) { check( std::abs(a-b)<tol, "numeric comparison" ); }
template< class F > void fails( F f ) { bool failed=false; try{f();}catch(std::runtime_error const&){failed=true;}check(failed,"expected a checked failure"); }
bm::Row row( std::initializer_list<bm::PairCoefficient> c, double area=1 ) { return bm::Row{c,area}; }

int main( int argc, char ** argv )
{
  // Tiny text protocol for an independent NumPy/SciPy SVD comparison runner.
  if( argc > 1 && std::string(argv[1]) == "--solve" )
  {
    std::size_t n,m,k,f; std::cin>>n>>m>>k>>f;
    std::vector<double> mass(n),v(n),cut(n,1e-30);
    for(auto& x:mass)std::cin>>x;
    for(auto& x:v)std::cin>>x;
    std::vector<bm::Relation> bc(k);
    for(auto& r:bc){long long b;std::cin>>r.a>>b>>r.sign>>r.offset;r.b=b<0?bm::invalidIndex:std::size_t(b);}
    std::vector<std::pair<std::size_t,std::size_t>> ties(f);
    for(auto& p:ties)std::cin>>p.first>>p.second;
    std::vector<bm::Row> rows(m);
    for(auto& r:rows){std::size_t count;std::cin>>r.measure>>count;r.coefficients.resize(count);for(auto& c:r.coefficients)std::cin>>c.a>>c.b>>c.value;}
    bm::Options opt; opt.relativeTolerance=1e-9;opt.rankTolerance=1e-13;
    try{
      auto r=bm::project(mass,v,cut,bc,ties,rows,opt);
      std::cout.precision(17);
      for(double x:r.velocity)std::cout<<x<<" ";
      std::cout<<"\n";
      for(double x:r.reactionImpulse)std::cout<<x<<" ";
      std::cout<<"\n";
      std::cout<<r.maximumJump<<" "<<r.boundaryResidual<<" "<<r.rank<<"\n";
    }catch(std::exception const& e){std::cerr<<e.what()<<"\n";return 2;}
    return 0;
  }
  bm::Options o;
  {
    auto r=bm::project({2,3},{2,-1},{0,0},{},{},{row({{0,1,1}})},o);
    near(r.velocity[0],.2);near(r.velocity[1],.2);
    near(2*r.velocity[0]+3*r.velocity[1],1);
    near(r.reactionImpulse[0],0);near(r.reactionImpulse[1],0);
  }
  {
    // Redundant mortar rows, arbitrary row scaling: same mass-minimum answer.
    auto r=bm::project({2,3},{2,-1},{0,0},{},{},
      {row({{0,1,1}}),row({{0,1,7}},7),row({{0,1,.001}},.001)},o);
    near(r.velocity[0],.2);check(r.rank==1,"rank of redundant constraints");
  }
  {
    // Positive-mass donor with a massless dependent buffer remains represented.
    auto r=bm::project({3,0,4,0},{1,-1,-2,2},{0,0,0,0},
      {{1,0,-1,0},{3,2,-1,0}},{},{row({{1,3,.5}})},o);
    near(r.velocity[0],-5.0/7);near(r.velocity[1],5.0/7);
    near(r.velocity[2],-5.0/7);check(r.accepted[0]==1,"fold support before cutoff");
  }
  {
    // Moving buffer: q, U, 2U-q. Both donor and buffer carry physical mass.
    auto r=bm::project({3,2,7},{1,.4,-.2},{0,0,0},
      {{1,bm::invalidIndex,1,.4},{2,0,-1,.8}},{},{},o);
    near(r.velocity[0],1);near(r.velocity[2],-.2);near(r.velocity[1],.4);
  }
  {
    // Common-velocity fallback participates in the solve, not an overwrite.
    auto r=bm::project({1,2,3,4},{1,-2,3,-4},{0,0,0,0},{},{{0,1}},
      {row({{0,1,.4},{2,3,.6}})},o);
    near(r.velocity[0],r.velocity[1]);near(r.velocity[2],r.velocity[3]);
    check(r.fallbackDof[0]==1,"fallback retained");
  }
  {
    // A completely fixed, already satisfied row has no multiplier unknown.
    auto r=bm::project({1,1},{7,9},{0,0},
      {{0,bm::invalidIndex,1,.25},{1,bm::invalidIndex,1,.25}},{},{row({{0,1,1}})},o);
    near(r.velocity[0],.25);check(r.rank==0,"redundant fully prescribed row");
    fails([&]{bm::project({1,1},{0,0},{0,0},
      {{0,bm::invalidIndex,1,0},{1,bm::invalidIndex,1,1}},{},{row({{0,1,1}})},o);});
  }
  {
    bm::AffineMap map(3,1e-12);
    map.relate(0,1,-1,2);map.relate(1,2,-1,4);map.fix(2,3);
    auto r=bm::reduce(map,{1,1,1},{0,0,0});near(r.value[0],1);near(r.value[1],1);near(r.value[2],3);
    fails([&]{map.fix(0,5);});
    bm::AffineMap cycle(2,1e-12);cycle.relate(0,1,-1,2);cycle.relate(1,0,1,0);
    auto c=bm::reduce(cycle,{1,1},{0,0});near(c.value[0],1);near(c.value[1],1);
  }
  {
    bm::Options tiny=o;tiny.maximumMatrixEntries=1;
    fails([&]{bm::project({1,1},{0,1},{0,0},{},{},{row({{0,1,1}})},tiny);});
    fails([&]{bm::project({-1,1},{0,1},{0,0},{},{},{row({{0,1,1}})},o);});
  }
  {
    bm::Faces faces;faces.cells={{10,10,10}};faces.low={{0,0,0}};faces.high={{1,1,1}};
    for(int d=0;d<3;++d){faces.kind[2*d]=1;faces.controlled[2*d][d]=true;}
    std::set<bm::LatticeKey> keys{{{2,2,2}}};bm::closeBoundaryDependencies(keys,faces);
    check(keys.size()==8,"three-plane ghost/donor closure");
    std::vector<bm::LatticeKey> list(keys.begin(),keys.end());
    for(int d=0;d<3;++d)
    {
      auto rel=bm::makeBoundaryRelations(list,2,d,faces);
      check(rel.size()==8,"component-wise independent fields");
      std::vector<double> m(16,1),v(16,2),cut(16,0);
      auto r=bm::project(m,v,cut,rel,{}, {},o);
      near(r.boundaryResidual,0);
      for(auto const& c:rel)near(r.velocity[c.a],-r.velocity[c.b]);
    }
    auto img=bm::symmetryImages({{{.2,.3,.4}},{{.2,.3,.4}}},faces,1,1e-12);
    check(img.size()==8,"octant symmetry image closure");
    std::array<double,3> sum{{0,0,0}};
    for(auto const& x:img)for(int d=0;d<3;++d)sum[d]+=x.normal[d];
    for(double x:sum)near(x,0);
    check(faces.boundaryMask({{0,.1,.2}},1e-12)==1,"physical x=0 mask");
    check(!faces.inside({{-.01,.1,.2}}),"exclude padded integration cells");
    check(faces.inside({{0,.1,.2}}),"retain terminal facets");
    faces.periodic[0]=true;
    check(faces.boundaryMask({{0,.1,.2}},1e-12)==0,"periodic face not physical boundary");
  }
  {
    // Admissible affine weak jump: vA-vB=x+.3y-.2 on the plane phi=0.
    // Coefficients integrate four samples with weighted phi exactly zero.
    std::vector<double> mass(8,1),cut(8,0),trial{-.2,0,.8,0,-.2,0,.8,0};
    auto r=bm::project(mass,trial,cut,{}, {},
      {row({{0,1,.4},{2,3,.1},{4,5,.4},{6,7,.1}})},o);
    for(std::size_t i=0;i<trial.size();++i)near(r.velocity[i],trial[i]);
  }
  {
    // Actual cubic B-spline boundary trace: two fields, three symmetry planes,
    // three opposite Moving planes, and overlapping mortar constraints. The
    // normal trace must remain prescribed after the coupled correction, even
    // though the adjacent material-side nodal velocities are independent.
    bm::Faces f; f.cells={{6,6,6}}; f.low={{0,0,0}}; f.high={{1,1,1}};
    for(int d=0; d<3; ++d)
    {
      f.kind[2*d]=1; f.kind[2*d+1]=2;
      f.controlled[2*d][d]=f.controlled[2*d+1][d]=true;
      f.target[2*d+1][d]=.01*(d+1);
    }
    std::vector<bm::LatticeKey> keys;
    std::map<bm::LatticeKey,std::size_t> index;
    for(int i=0;i<=6;++i)for(int j=0;j<=6;++j)for(int k=0;k<=6;++k)
    { bm::LatticeKey x{{i,j,k}}; index[x]=keys.size(); keys.push_back(x); }
    std::vector<double> mass(2*keys.size()), trial(mass.size()), cut(mass.size(),0);
    for(std::size_t i=0;i<mass.size();++i){mass[i]=1+.1*(i%7);trial[i]=.03*std::sin(double(i));}
    std::vector<bm::Row> rows(5);
    for(std::size_t a=0;a<rows.size();++a)
    {
      rows[a].measure=1;
      for(std::size_t i=a;i<keys.size();i+=7)
        rows[a].coefficients.push_back({2*i,2*i+1,1./49});
    }
    for(int d=0;d<3;++d)
    {
      auto bc=bm::makeBoundaryRelations(keys,2,d,f);
      auto result=bm::project(mass,trial,cut,bc,{},rows,o);
      near(result.maximumJump,0);
      for(int side=0;side<2;++side)for(int j=0;j<=6;++j)for(int k=0;k<=6;++k)
      {
        bm::LatticeKey boundary{{0,0,0}};
        boundary[d]=side?5:1; boundary[(d+1)%3]=j; boundary[(d+2)%3]=k;
        auto ghost=boundary, donor=boundary;
        ghost[d]=side?6:0; donor[d]=side?4:2;
        for(int field=0;field<2;++field)
        {
          double const trace=(result.velocity[2*index.at(ghost)+field]
            +4*result.velocity[2*index.at(boundary)+field]
            +result.velocity[2*index.at(donor)+field])/6;
          near(trace,f.target[2*d+side][d]);
        }
      }
    }
    // No coefficient clipping: GEOS' partition-of-unity closure can yield a
    // tiny negative coefficient in an otherwise nonnegative trace row.
    auto q=bm::project({1,1,1,1},{1,2,3,4},{0,0,0,0},{},{},
      {row({{0,1,1.0},{2,3,-2e-17}})},o);
    near(q.maximumJump,0);
    fails([&]{bm::project({1,1},{0,1},{0,0},{},{{0,99}},{},o);});
  }
  std::cout<<"PASS: "<<checks<<" standalone production-header checks\n";
}
