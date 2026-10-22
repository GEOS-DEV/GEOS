int main()
{
  constexpr double pi=3.14159265358979323846;
  // Exact annular/cylindrical area, nodal radial first moment and axis limit.
  double ny[3]={0,1,0}, nr[3]={1,0,0}, zero[3]={};
  for(double h: {1e-4,.01,1.}) for(double r: {0.,.25*h,h,10*h})
  {
    double const expected=r>=h ? 2*pi*r*h : 2*pi*(h*h/6+r*h/2+r*r/2-r*r*r/(6*h));
    close(mpmMeasures::axisymmetricLineArea(h,h,r,ny,zero),expected);
    close(mpmMeasures::axisymmetricNodalVolume(h,h,r),expected*h);
    close(mpmMeasures::lineArea(h,h,.3,ny,zero),.3*h);
    if(r>0) close(mpmMeasures::axisymmetricLineArea(h,h,r,nr,zero),2*pi*r*h);
    else close(mpmMeasures::axisymmetricLineArea(h,h,r,nr,zero),0);
    double pos[3]={.2*h,0,0};
    close(mpmMeasures::axisymmetricLineArea(h,h,r,nr,pos),2*pi*(r+.2*h)*h*.8);
  }
  // Independently integrate random clipped inclined lines with a fine midpoint
  // rule; covers normal reversal/rescaling, offset, radial support crossing.
  std::mt19937 gen(723); std::uniform_real_distribution<double> rand(-1.,1.);
  for(int c=0;c<60;++c)
  {
    double angle=3*rand(gen),n[3]={std::cos(angle),std::sin(angle),0};
    double pos[3]={.7*rand(gen),.5*rand(gen),0},r=.2+.7*(rand(gen)+1);
    double distance=n[0]*pos[0]+n[1]*pos[1],sum=0;
    constexpr int samples=200000;
    for(int k=0;k<samples;++k)
    {
      double t=-2.+4.*(k+.5)/samples;
      double x=distance*n[0]-t*n[1],y=distance*n[1]+t*n[0];
      if(std::abs(x)<1 && std::abs(y)<.7 && r+x>0)sum+=2*pi*(r+x)*(1-std::abs(x))*(1-std::abs(y)/.7)*4/samples;
    }
    double area=mpmMeasures::axisymmetricLineArea(1,.7,r,n,pos);
    close(area,sum,2e-8);
    double flipped[3]={-3*n[0],-3*n[1],0};
    close(area,mpmMeasures::axisymmetricLineArea(1,.7,r,flipped,pos));
  }
  // Delay CZ activation until a non-identity F; rebase with noncommuting
  // meridional shear/stretch and an independent hoop stretch.
  for(int c=0;c<120;++c)
  {
    double ref[1][3][3]={{{1.3,.2,0},{-.1,.9,0},{0,0,1.7}}};
    double F[1][3][3]; LvArray::tensorOps::copy<3,3>(F[0],ref[0]);
    double cof[3][3];relativeArea(F,ref,cof);
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) close(cof[i][j],i==j?1:0,1e-10);
    double rel[3][3]={{1.2+.1*rand(gen),.3*rand(gen),0},{.3*rand(gen),.8+.1*rand(gen),0},{0,0,1.1+.2*rand(gen)}};
    LvArray::tensorOps::Rij_eq_AikBkj<3,3,3>(F[0],rel,ref[0]);
    relativeArea(F,ref,cof);
    double expected[3][3]; LvArray::tensorOps::cofactor<3>(expected,rel);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)close(cof[i][j],expected[i][j]);
    double fresh[3][3]={{.9,.1,0},{0,1.4,0},{0,0,.8}};
    rebase(ref,F[0],fresh); LvArray::tensorOps::copy<3,3>(F[0],fresh);
    relativeArea(F,ref,cof);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)close(cof[i][j],expected[i][j]);
  }
  // TypeB uses material-domain ray intersections for quads AND triangles.
  // Barycentric/edge coordinates independently identify the boundary anchor.
  for(int shape: {3,4}) for(double angle: {.1,.7,1.7,3.2,4.9})
  {
    double g[2][2]={{.8,.15},{.2,.6}},n[3]={std::cos(angle),std::sin(angle),0},x[3];
    anchor(shape,g,n,x);
    double det=g[0][0]*g[1][1]-g[0][1]*g[1][0];
    double u=(x[0]*g[1][1]-x[1]*g[1][0])/det,v=(g[0][0]*x[1]-g[0][1]*x[0])/det;
    if(shape==3)close(std::max(std::abs(u),std::abs(v)),.5);
    else {assert(u>=-1./3-1e-12 && v>=-1./3-1e-12 && u+v<=1./3+1e-12);close(std::min({u+1./3,v+1./3,1./3-u-v}),0);}
    close(x[0]*n[1]-x[1]*n[0],0);
  }
  // Actual production state kernel and actual Uncoupled law: mixed-mode
  // forces with full 3D cofactor, hoop area growth, damage and compression.
  int cases=0;
  for(double angle: {0.,.35,pi/2}) for(double r: {.001,.1,10.})
  for(double hoop: {.7,1.,1.5}) for(double opening: {-.02,0.,.02,.12})
  for(double slip: {0.,.01,.12})
  {
    double damage=0; Law law{200,100,.1,.1,&damage};
    double masses[1][2]={{1,3}},refArea[1][2]={{2*pi*r*.03,2*pi*r*.03}};
    double normals0[1][2][3]={{{std::cos(angle),std::sin(angle),0},{-std::cos(angle),-std::sin(angle),0}}};
    double normals[1][2][3]={},cof[1][2][3][3]={},disp[1][2][3]={},force[1][2][3]={};
    double F[3][3]={{1.2,.2,0},{0,.9,0},{0,0,hoop}};
    for(int side=0;side<2;++side) {
      LvArray::tensorOps::cofactor<3>(cof[0][side],F);
      LvArray::tensorOps::Ri_eq_AijBj<3,3>(normals[0][side],cof[0][side],normals0[0][side]);
      double scale=LvArray::tensorOps::l2Norm<3>(normals[0][side]);
      LvArray::tensorOps::scale<3>(normals[0][side],1/scale);
      double n[3]={normals[0][0][0],normals[0][0][1],0};
      for(int d=0;d<2;++d)disp[0][side][d]=(side==0?-.5:.5)*opening*n[d]+(side==0?.5:-.5)*slip*(d==0?-n[1]:n[0]);
    }
    CohesiveZoneStateUpdateKernel::launch<Serial>(1,law,.001,1,1e-20,0,0,1,0,0,0,0,0,0,masses,disp,normals,cof,force,normals0,refArea);
    double areaVector[3];LvArray::tensorOps::Ri_eq_AijBj<3,3>(areaVector,cof[0][0],normals0[0][0]);
    double area=refArea[0][0]*LvArray::tensorOps::l2Norm<3>(areaVector);
    bool failed=opening>=.1 || slip>=.1;
    close(damage,failed?1:0);
    for(int d=0;d<3;++d) {
      close(force[0][0][d]+force[0][1][d],0);
      double n=normals[0][0][d],t=d==0?-normals[0][0][1]:(d==1?normals[0][0][0]:0);
      close(force[0][0][d],failed?0:area*(200*opening*n-100*slip*t));
    }
    ++cases;
  }
  // Compression may be delegated to ordinary contact without suppressing
  // cohesive shear. Damage is irreversible when the interface recloses.
  for(int prevent: {0,1})
  {
    double damage=0; Law law{200,100,.1,.1,&damage};
    double masses[1][2]={{1,1}},areas[1][2]={{2,2}};
    double normals[1][2][3]={{{0,1,0},{0,-1,0}}},disp[1][2][3]={{{.01,.02,0},{0,0,0}}};
    double cof[1][2][3][3]={{{{1,0,0},{0,1,0},{0,0,1}},{{1,0,0},{0,1,0},{0,0,1}}}};
    double force[1][2][3]={};
    CohesiveZoneStateUpdateKernel::launch<Serial>(1,law,.001,1,1e-20,prevent,0,1,0,0,0,0,0,0,masses,disp,normals,cof,force,normals,areas);
    close(force[0][0][0],-2);close(force[0][0][1],prevent?0:-8);
    double normalStress=0,shearStress=0;
    law.jumpDisplacementUpdate(0,.12,0,normalStress,shearStress);
    law.jumpDisplacementUpdate(0,0,0,normalStress,shearStress);
    close(damage,1);close(normalStress,0);
  }
  // Cached mappings survive a partial failure; complete failure releases the
  // particle to ordinary contact. Node ordering must not affect this result.
  for(int failed=0;failed<2;++failed) {
    double damage[2]={};damage[failed]=1;int flag;double force;
    partialGather(damage,flag,force);assert(flag==1);close(force,failed==0?3:2);
  }
  {double damage[2]={1,1};int flag;double force;partialGather(damage,flag,force);assert(flag==0);close(force,0);}
  std::cout<<"Axisymmetric CZ: ring integrals, 120 reset/activation cases, "<<cases<<" production-law cases, partial failure passed\n";
}
