// Tests of production reconstruction, field mapping, separability and contact.
double scalarField(std::vector<Sample> const & samples,double r,double y,double R) {
  double value=0,weight=0;
  for(auto const &p:samples) for(int sign:{-1,1}) {
    double w=p.area*kernelDevice(std::hypot(r-sign*p.r,y-p.y),R,1);
    weight+=w;value+=w*(markSurfaceAsDamage(p.flag)?1.:p.d);
  }
  return weight>0?value/weight:0;
}
void gradients() {
  std::vector<Sample> samples;
  for(int i=0;i<18;++i) for(int j=-12;j<=12;++j)
    samples.push_back({.05+i*.1,j*.1,.01,.4+.2*std::sin(j*.1)+.1*std::cos(i*.1)});
  // Finite differences of the existing normalized scalar kernel provide an
  // independent derivative check at the axis, near it and away from it.
  int cases=0;
  for(double r:{0.,1e-8,.05,.35,1.1}) for(double y:{-.21,0.,.17}) {
    auto g=reconstruct(samples,r,y,.47);
    double const e=1e-7; // C1 kernel cutoff gives O(e) error for samples exactly at support radius.
    assert(std::abs(g[0]-(scalarField(samples,r+e,y,.47)-scalarField(samples,r-e,y,.47))/(2*e))<2e-8);
    double const finiteY=(scalarField(samples,r,y+e,.47)-scalarField(samples,r,y-e,.47))/(2*e);
    if(std::abs(g[1]-finiteY)>=2e-8) std::cerr<<"r="<<r<<" y="<<y<<" gradient="<<g[1]<<" finite difference="<<finiteY<<'\n';
    assert(std::abs(g[1]-finiteY)<2e-8);
    close(g[2],0);++cases;
  }
  close(reconstruct(samples,0,.12,.47)[0],0);
  for(auto &p:samples) p.d=.37;
  for(double r:{0.,.05,.4,1.72}) for(double y:{-1.19,0.,1.19}) {
    auto g=reconstruct(samples,r,y,.47);for(double v:g) close(v,0);
  }
  // Radius-independent planar pattern: translating away from the axis changes
  // ring volumes but must not rotate the recovered meridional crack normal.
  std::vector<Sample> plane;
  for(int i=-6;i<=6;++i) for(int j=-6;j<=6;++j)
    plane.push_back({3.+i*.1,j*.1,.01,.5+.12*i*.1+.2*j*.1});
  auto g=reconstruct(plane,3.,0.,.51);
  close(g[0]/g[1],.6);
  for(auto &p:plane) {p.r+=20;p.area*=.31;}
  auto moved=reconstruct(plane,23.,0.,.51);
  for(int d=0;d<3;++d) close(g[d],moved[d]);
  // The existing surface-as-damage convention remains in use.
  samples.push_back({.3,.1,100.,1.,0,2});
  auto surface=reconstruct(samples,.3,0.,.47);assert(surface[1]>0);
  auto empty=reconstruct({},.3,0.,.47);for(double v:empty) close(v,0);
  std::cout<<cases<<" normalized-kernel derivative cases; constant fields, axis parity, radius invariance and surface flags passed.\n";
}

void fieldPairs() {
  Vec zero{},reference{0,1,0},plus{0,1,0},minus{0,-1,0};
  constexpr int groups=3;
  // Arbitrary body IDs survive DFG's side*groups+group mapping. Explicit
  // meridional normals take precedence over the damage gradient as before.
  for(int group=0;group<groups;++group) {
    assert(partitionField(groups,0,group,minus.data(),zero.data(),reference.data())==group);
    assert(partitionField(groups,1,group,plus.data(),zero.data(),reference.data())==group);
    assert(partitionField(groups,1,group,minus.data(),zero.data(),reference.data())==group+groups);
    assert(partitionField(groups,1,group,plus.data(),minus.data(),reference.data())==group+groups);
  }
  double a[6]={0,2,0,0,0,0},b[6]={0,5,0,0,0,0};
  close(computeSurfaceQualityFromMappingNormalTensor(1,a,b,2,5),1);
  auto separates=[&](int A,int B,double d,double peak,double quality=1.,double single=0.,double threshold=-1.) {
    return evaluateSeparabilityCriterion(groups,threshold,.5,DBL_MAX,1.,.8,A,B,
      d,d,peak,peak,zero.data(),plus.data(),quality,single);
  };
  for(int group=0;group<groups;++group) {
    assert(!separates(group,group+groups,0.,0.));
    assert(!separates(group,group+groups,.8,.9));
    assert(separates(group,group+groups,.8,1.));
    assert(!separates(group,group+groups,.8,1.,.1));
    assert(!separates(group,group+groups,.8,1.,1.,.8,.5));
    // Damage-triggered opening: intact fields tie; separable receding fields
    // receive no attractive impulse, on both force and FMPM impulse paths.
    for(bool forcePath:{false,true}) {
      auto tied=impulse(forcePath,plus,minus,plus,2,5,.2,Gap::Simple,0,false);
      close(minus[1]+tied[1]/2,plus[1]-tied[1]/5);
      auto opened=impulse(forcePath,plus,minus,plus,2,5,.2,Gap::Simple,0,true);
      for(double value:opened) close(value,0);
    }
  }
  for(int A=0;A<6;++A) for(int B=A+1;B<6;++B)
    assert(separates(A,B,0.,0.)==(A%groups!=B%groups));
  // Three contacting bodies at one node: accumulate all pair impulses using
  // the production law. Total momentum is conserved for unequal ring masses.
  for(bool forcePath:{false,true}) {
    double mass[3]={1,2,4};Vec v[3]={{0,1,0},{0,0,0},{0,-1,0}},sum[3]={};
    for(int A=0;A<3;++A) for(int B=A+1;B<3;++B) {
      auto j=impulse(forcePath,plus,v[A],v[B],mass[A],mass[B],0,Gap::Simple,0);
      for(int d=0;d<3;++d) {sum[A][d]+=j[d];sum[B][d]-=j[d];}
    }
    for(int d=0;d<3;++d) close(sum[0][d]+sum[1][d]+sum[2][d],0);
  }
  std::cout<<"Three groups / six DFG fields: partitioning, damage/quality gates, intact tie, crack opening and contact momentum passed.\n";
}
// Integration adapter for the generated example's central band. Uses actual
// production mapping, reconstruction, partition and separability functions;
// a serial scatter reproduces the documented mass/volume-weighted sums.
void exampleBand(char const *path) {
  std::ifstream file(path);assert(file);
  std::string line,name;std::getline(file,line);std::istringstream header(line);
  std::vector<std::string> names;while(header>>name) names.push_back(name);
  auto column=[&](char const *key) {auto it=std::find(names.begin(),names.end(),key);assert(it!=names.end());return it-names.begin();};
  struct Particle {Sample sample;double g[3][3]={};ParticleType type;};
  std::vector<Particle> particles;std::vector<Sample> samples;
  int groups=1;bool damaged=false;
  while(std::getline(file,line)) {
    if(line.empty()) continue;
    std::istringstream row(line);std::vector<double> value(names.size());
    for(auto &v:value) {row>>v;assert(row);}
    Particle p;auto &s=p.sample;
    s.r=value[column("PositionX")];s.y=value[column("PositionY")];
    s.area=value[column("ParticleVolume")]/(axisymmetric::twoPi*s.r);
    s.d=value[column("Damage")];s.group=int(value[column("ContactGroup")]);
    s.flag=int(value[column("SurfaceFlag")]);assert(s.flag==0);
    p.g[0][0]=value[column("GeometryXX")];p.g[0][1]=value[column("GeometryXY")];
    p.g[1][0]=value[column("GeometryYX")];p.g[1][1]=value[column("GeometryYY")];p.g[2][2]=1;
    p.type=int(value[column("ParticleType")])==3 ? ParticleType::CPTI : ParticleType::CPDI;
    groups=std::max(groups,s.group+1);damaged=damaged || s.d>.99;
    particles.push_back(p);samples.push_back(s);
  }
  assert(!particles.empty());
  int map[9][9][9],cells[3]={8,8,0};
  double origin[3]={-.001,-.004,0},h[3]={.001,.001,1};
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<9;++k) map[i][j][k]=(i*9+j)*9+k;
  Kernel mapper;mapper.m_axisymmetric=true;
  struct Contribution {int node,group;double weight,damage;Vec gradient;};
  std::vector<Contribution> contributions;
  std::array<Vec,729> reference{};
  for(auto const &p:particles) if(std::abs(p.sample.y)<.0026) {
    auto const &s=p.sample;auto gradient=reconstruct(samples,s.r,s.y,.0015);
    int mapped[16];double values[16],grad[16][3],center[3]={s.r,s.y,0};
    assert(mapper.computeParticleMapping<2>(p.type,center,p.g,map,origin,h,cells,mapped,values,grad));
    double volume=axisymmetric::volume(s.area,s.r);
    for(int a=0;a<(p.type==ParticleType::CPTI?12:16);++a) {
      int node=mapped[a];
      auto norm=[](Vec const &v){return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];};
      if(norm(gradient)>norm(reference[node]) || (!(norm(gradient)<norm(reference[node])) && gradient>reference[node])) reference[node]=gradient;
      contributions.push_back({node,s.group,volume*values[a],s.d,gradient});
    }
  }
  struct Field {double volume=0,damage=0,peak=0,tensor[6]={};};
  std::vector<std::array<Field,6>> fields(729);Vec zero{};
  for(auto const &c:contributions) {
    int side=partitionField(groups,1,c.group,c.gradient.data(),zero.data(),reference[c.node].data());
    auto &f=fields[c.node][side];f.volume+=c.weight;f.damage+=c.weight*c.damage;f.peak=std::max(f.peak,c.damage);
    double length=std::hypot(c.gradient[0],c.gradient[1]);
    if(length>.0001/h[0]) {
      double x=c.gradient[0]/length,y=c.gradient[1]/length;
      f.tensor[0]+=c.weight*x*x;f.tensor[1]+=c.weight*y*y;f.tensor[3]+=c.weight*x*y;
    }
  }
  int separated=0,centralGroup=groups==1?0:1;
  for(int i=1;i<=7;++i) {
    auto const &a=fields[map[i][4][0]][centralGroup],&b=fields[map[i][4][0]][centralGroup+groups];
    if(a.volume<=1e-30 || b.volume<=1e-30) continue;
    double quality=computeSurfaceQualityFromMappingNormalTensor(1,a.tensor,b.tensor,a.volume,b.volume);
    bool separate=evaluateSeparabilityCriterion(groups,-1.,.5,DBL_MAX,.0015,.5,centralGroup,centralGroup+groups,
      a.damage/a.volume,b.damage/b.volume,a.peak,b.peak,zero.data(),zero.data(),quality,0.);
    if(separate) ++separated;
  }
  assert(damaged ? separated==7 : separated==0);
  std::cout<<path<<": "<<separated<<" of 7 central-plane nodes permit opening under the unchanged criterion.\n";
}
int main(int argc,char **argv) {gradients();fieldPairs();for(int i=1;i<argc;++i) exampleBand(argv[i]);}
