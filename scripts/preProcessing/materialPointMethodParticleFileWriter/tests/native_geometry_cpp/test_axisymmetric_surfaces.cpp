// Header kernels: exact measures, inverse geometry and meridional regression.
bool near(double a,double b,double tolerance=2e-10)
{return std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)});}
void testMeasures()
{
  using namespace axisymmetricSurface;
  double lower[2]={0,-100},upper[2]={100,100};bool symmetry[4]={true,false,false,false};
  std::mt19937 rng(314159);std::uniform_real_distribution<double> uniform(0,1);
  for(int k=0;k<150;++k)
  {
    double radius=4*uniform(rng),hx=.3+uniform(rng),hy=.3+uniform(rng),angle=6.28*uniform(rng);
    auto s=support(hx,hy,radius,0,lower,upper,symmetry);
    double n[3]={std::cos(angle),std::sin(angle),0},opposite[3]={-n[0],-n[1],0};
    double bound=std::abs(n[0])*hx+std::abs(n[1])*hy,offset=(uniform(rng)*1.6-.8)*bound;
    double v=halfSpaceVolume(s,n,offset),complement=halfSpaceVolume(s,opposite,-offset);
    assert(near(v+complement,s.volume()));
    assert(v>=0 && v<=s.volume()*(1+1e-12));
    if(v>1e-8*s.volume() && v<(1-1e-8)*s.volume())
      assert(near(offsetFromVolume(s,n,v),offset));
    double point[3]={offset*n[0],offset*n[1],0},eps=1e-5*std::min(hx,hy);
    double derivative=(halfSpaceVolume(s,n,offset+eps)-halfSpaceVolume(s,n,offset-eps))/(2*eps);
    assert(near(derivative,lineArea(s,n,point),2e-8));
  }
  // Axis support is finite, although the area of the revolution axis is zero.
  auto axis=support(1,1,0,0,lower,upper,symmetry);
  assert(near(axis.volume(),twoPi/6));
  double radial[3]={1,0,0},zero[3]={};
  assert(near(lineArea(axis,radial,zero),0));
  assert(near(halfSpaceVolume(axis,radial,0),0)); // never assume a half-full axis node
  // A radius-2 cylinder between axial symmetry planes has area 2*pi*R*H.
  double area=0,h=.25,H=3,R=2;
  lower[1]=0;upper[1]=H;symmetry[2]=symmetry[3]=true;
  for(int i=0;i<=12;++i) for(int j=0;j<=12;++j)
  {
    auto s=support(h,h,i*h,j*h,lower,upper,symmetry);
    double position[3]={R-i*h,0,0};area+=lineArea(s,radial,position);
  }
  assert(near(area,twoPi*R*H));
  // Inclined line revolved about y: exact frustum area, including r variation.
  area=0;R=1;double slope=.25,norm=std::hypot(1.,slope),n[3]={1/norm,-slope/norm,0};
  for(int i=0;i<=12;++i) for(int j=0;j<=12;++j)
  {
    auto s=support(h,h,i*h,j*h,lower,upper,symmetry);
    double d=(R+slope*j*h-i*h)/norm,position[3]={d*n[0],d*n[1],0};area+=lineArea(s,n,position);
  }
  assert(near(area,.5*twoPi*(2*R+slope*H)*H*norm));
  std::cout<<"150 exact clipped-volume/area derivative/inverse checks; axis, cylinder and frustum areas passed.\n";
}

void testRegression()
{
  using namespace axisymmetricSurface;
  for(double length:{1e-6,1.,1e6}) for(double angle:{0.,.4,1.5707963267948966})
  {
    double reference[2]={std::cos(angle),std::sin(angle)},n[3]={1,1,0};
    // Two separated strips, reflected/rotated; each particle samples equal area.
    auto visit=[&](auto add){for(int i=-6;i<=6;++i) for(int j=1;j<=5;++j) for(int sign:{-1,1})
    {
      double a=sign*(.12+.08*j),b=i*.15;
      add(a*reference[0]-b*reference[1],a*reference[1]+b*reference[0],double(sign),.01);
    }};
    assert(logisticNormal(visit,length,length,80,1e-7,n));
    assert(n[0]*reference[0]+n[1]*reference[1]>1-1e-8);
    double refined[3]={1,1,0};
    auto splitVisit=[&](auto add){visit([&](double x,double y,double label,double weight)
      {for(int k=0;k<4;++k)add(x,y,label,weight/4);});};
    assert(logisticNormal(splitVisit,length,length,80,1e-7,refined));
    assert(near(n[0],refined[0])&&near(n[1],refined[1]));
  }
  // Axis images make a horizontal interface normal purely axial.
  double n[3]={.2,1,0};
  auto visit=[&](auto add){for(int i=0;i<8;++i)for(int j=-4;j<=4;++j)for(int image:{-1,1})
  {if(j!=0)add(image*(i+.5)/8,j*.15,j<0?-1.:1.,.01);}};
  assert(logisticNormal(visit,.01,.01,80,1e-7,n));assert(std::abs(n[0])<2e-5&&n[1]>.999999);
  double unchanged[3]={.6,.8,0};
  auto oneSide=[](auto add){add(0.,0.,1.,1.);};
  assert(!logisticNormal(oneSide,1,1,50,1e-5,unchanged));assert(near(unchanged[0],.6));
  std::cout<<"Meridional logistic normals: length scaling, refinement weights, axis images and one-sided fallback passed.\n";
}

void testSurfaceStress()
{
  for(auto type:{ParticleType::CPDI,ParticleType::CPTI,ParticleType::SinglePointBSpline,ParticleType::SinglePoint})
  {
    auto full=std::make_unique<SurfaceKernel>(12);full->addRectangle(12,1,false,type);full->project();
    for(auto const &p:full->particles)full->assemble(p,0,-1,1);
    double maxForce=0;
    for(auto const &f:full->gridSurfaceTensionForce)for(auto const &v:f)for(auto x:v)maxForce=std::max(maxForce,std::abs(x));
    std::cout<<"Uniform phase max force: "<<maxForce<<"\n";assert(maxForce<1e-9);
    for(int n:{12,24})
    {
      auto cylinder=std::make_unique<SurfaceKernel>(n);cylinder->m_boundaryConditionTypes[1]=0;
      cylinder->addRectangle(n,.5,false,type);cylinder->project();
      for(auto const &p:cylinder->particles)cylinder->assemble(p,0,-1,1);
      double reconstructedArea=0;
      for(int i=3;i<=n+3;++i)for(int j=3;j<=n+3;++j)
      {
        int node=cylinder->m_ijkMap[i][j][0];auto support=cylinder->axisymmetricSurfaceSupport(cylinder->gridPosition[node]);
        double volume=cylinder->gridMaterialVolume[node][0],cap=support.volume();
        double normal[3]={i==3?0:cylinder->boundary[node][0],(j==3||j==n+3)?0:cylinder->boundary[node][1],0};
        double norm=std::hypot(normal[0],normal[1]);
        if(volume<1e-12*cap||volume>(1-1e-12)*cap||norm<1e-12*cap/cylinder->m_hEl[0])continue;
        normal[0]/=norm;normal[1]/=norm;
        double d=axisymmetricSurface::offsetFromVolume(support,normal,volume),offset[3]={d*normal[0],d*normal[1],0};
        reconstructedArea+=axisymmetricSurface::lineArea(support,normal,offset);
      }
      assert(near(reconstructedArea,axisymmetricSurface::twoPi*.5));
      double area=axisymmetricSurface::twoPi*.5,ratio=-cylinder->work(false)/area;
      std::cout<<"Cylinder void capillary work / gamma*area (n="<<n<<"): "<<ratio<<"\n";
      assert(std::abs(ratio-1)<.15);
      double axialForce=0;for(auto const &f:cylinder->gridSurfaceTensionForce)for(auto const &v:f){axialForce+=v[1];assert(std::abs(v[2])<1e-15);}
      assert(std::abs(axialForce)<1e-9);
      auto pair=std::make_unique<SurfaceKernel>(n);pair->addRectangle(n,1,true,type);pair->project();
      for(auto const &p:pair->particles)pair->assemble(p,0,1,1);
      ratio=-pair->work(false)/area;
      std::cout<<"Cylinder phase pair capillary work / gamma*area: "<<ratio<<"\n";assert(std::abs(ratio-1)<.05);
    }
  }
}

void testGatherAndSphere()
{
  auto kernel=std::make_unique<SurfaceKernel>(24);
  double center[3]={.012,.52,0},geometry[3][3]={{.02,0,0},{.006,.023,0},{0,0,0}};
  int nodes[16];double weights[16],gradient[16][3];
  assert(kernel->computeParticleMapping<2>(ParticleType::CPDI,center,geometry,kernel->m_ijkMap,
    kernel->m_xLocalMin,kernel->m_hEl,kernel->m_nEl,nodes,weights,gradient));
  double normal[3]={.6,.8,0},distance=.031,result[3]={},old[3]={};
  for(int a=0;a<16;++a)
  {
    auto const &x=kernel->gridPosition[nodes[a]];
    double d=distance+normal[0]*(center[0]-x[0])+normal[1]*(center[1]-x[1]);
    double offset[3]={d*normal[0],d*normal[1],0};
    axisymmetricSurface::gatherOffset(weights[a],x,center,normal,offset,result);
    for(int k=0;k<2;++k)old[k]+=weights[a]*offset[k];
  }
  assert(near(result[0],distance*normal[0])&&near(result[1],distance*normal[1]));
  assert(std::hypot(old[0]-result[0],old[1]-result[1])>1e-5);
  std::cout<<"Particle-relative plane gather corrects cylindrical centroid bias.\n";
  for(auto type:{ParticleType::CPDI,ParticleType::CPTI,ParticleType::SinglePointBSpline})
  {
    double previous=1;
    for(int n:{16,32})
    {
      auto sphere=std::make_unique<SurfaceKernel>(n);
      sphere->m_boundaryConditionTypes[1]=sphere->m_boundaryConditionTypes[2]=sphere->m_boundaryConditionTypes[3]=0;
      sphere->addRectangle(n,1,false,type);
      auto &ps=sphere->particles;
      ps.erase(std::remove_if(ps.begin(),ps.end(),[](auto const&p){return std::hypot(p.x[0],p.x[1]-.5)>.35;}),ps.end());
      sphere->project();for(auto const&p:ps)sphere->assemble(p,0,-1,1);
      double area=2*axisymmetricSurface::twoPi*.35*.35,ratio=-sphere->work(true)/(2*area),error=std::abs(ratio-1);
      std::cout<<"Sphere capillary dilation work / (2*gamma*area), n="<<n<<": "<<ratio<<"\n";
      assert(error<.2&&error<previous);previous=error;
    }
  }
}
