void near(double a,double b){assert(std::abs(a-b)<3e-10*std::max(1.0,std::max(std::abs(a),std::abs(b))));}
ParticleSubRegion & addRegion(ParticleManager & pm,int n,ParticleType type=ParticleType::SinglePointBSpline,int dimensions=3)
{
 pm.regions.push_back(std::make_unique<ParticleSubRegion>());auto & r=*pm.regions.back();r.count=n;r.type=type;r.model.count=n;
 // CREATE_FIELDS
 r.add<Array<globalIndex,1>>("particleID");
 for(string name:{"particleRank","particleSurfaceFlag","particleGeometryShape","localToGlobalMap","ghostRank"})r.add<Array<int,1>>(name);
 for(string name:{"particleDamage","particleVolume","particleInternalEnergy","importedScalar"})r.add<Array<double,1>>(name);
 for(string name:{"particleCenter","particleVelocity","particleSurfacePosition","particleSurfaceTraction"})r.add<Array<double,2>>(name,3);
 r.add<Array<double,3>>("particleGeometryCurrent",dimensions,dimensions);
 r.getField<fields::mpm::particleGeometryReference>().resize(n,dimensions,dimensions);
 r.model.add<Array<double,4>>("history",2,2,3);
 r.model.add<Array<double,1>>("constantTable");r.model.wrappers.at("constantTable")->sized=false;
 r.neighbors.count=n;
 for(int p=0;p<n;++p)
 {
   r.getParticleID()[p]=100+p;r.get<Array<int,1>>("particleRank")[p]=0;r.getParticleGeometryShape()[p]=2;
   r.getField<fields::mpm::particleMass>()[p]=7+p;r.getParticleVolume()[p]=13+p;
   r.getField<fields::mpm::particleReferenceVolume>()[p]=17+p;r.get<Array<double,1>>("particleInternalEnergy")[p]=27+p;
   r.getField<fields::mpm::particleMeltFlag>()[p]=0.7;r.get<Array<double,1>>("importedScalar")[p]=123+p;
   auto hist=r.model.get<Array<double,4>>("history");for(int a=0;a<2;++a)for(int b=0;b<2;++b)for(int c=0;c<3;++c)hist[p][a][b][c]=1000+p*100+a*20+b*5+c;
   for(int i=0;i<3;++i)
   {
     r.getParticleCenter()[p][i]=1+i;r.getField<fields::mpm::particleReferencePosition>()[p][i]=4+i;
     r.getParticleVelocity()[p][i]=2-i;r.getParticleSurfacePosition()[p][i]=0.3+i;
     r.getField<fields::mpm::particleEstimatedSurfacePosition>()[p][i]=0.1+i;
     r.getParticleSurfaceTraction()[p][i]=3+i;r.getField<fields::mpm::particleReferenceSurfaceTraction>()[p][i]=6+i;
     for(int a=0;a<3;++a)
     {
       if(a<dimensions && i<dimensions)
       {
         r.getParticleGeometryCurrent()[p][a][i]=a==i?(a==0?2.3:0.6):0.0;
         r.getField<fields::mpm::particleGeometryReference>()[p][a][i]=a==i?0.7:0;
       }
       r.getField<fields::mpm::particleDeformationGradient>()[p][a][i]=a==i?1.2:0;
       r.getField<fields::mpm::particleVelocityGradient>()[p][a][i]=(a==i?0.2:0.1)*(a+1);
     }
     r.getField<fields::mpm::particleReferenceSurfacePosition>()[p][i]=r.getParticleSurfacePosition()[p][i]/1.2;
   }
   r.getParticleSurfaceFlag()[p]=1;
 }
 r.setActiveParticleIndices();return r;
}
void checkTransfer(ParticleSubRegion & r,int shape,int affineVelocity)
{
 assert(r.size()>1);double mass=0,volume=0,v0=0,momentum[3]={},center[3]={},referenceCenter[3]={},force[3]={},referenceForce[3]={};std::set<globalIndex> ids;
 double totalKE=0,expectedIncrease=0;
 for(int p=0;p<r.size();++p)
 {
   ids.insert(r.getParticleID()[p]);double m=r.getField<fields::mpm::particleMass>()[p];double fraction=m/7;
   mass+=m;volume+=r.getParticleVolume()[p];v0+=r.getField<fields::mpm::particleReferenceVolume>()[p];
   assert(r.getParticleGeometryShape()[p]==shape);near(r.getField<fields::mpm::particleMeltFlag>()[p],0.7);
   near(r.get<Array<double,1>>("particleInternalEnergy")[p],27);near(r.get<Array<double,1>>("importedScalar")[p],123);
   auto hist=r.model.get<Array<double,4>>("history");for(int a=0;a<2;++a)for(int b=0;b<2;++b)for(int c=0;c<3;++c)near(hist[p][a][b][c],1000+a*20+b*5+c);
   double speed2=0;
   for(int i=0;i<3;++i)
   {
     double velocity=r.getParticleVelocity()[p][i];speed2+=velocity*velocity;
     double increment=0;for(int a=0;a<3;++a)increment+=r.getField<fields::mpm::particleVelocityGradient>()[p][i][a]*(r.getParticleCenter()[p][a]-(1+a));
     near(velocity,(2-i)+(affineVelocity?increment:0));expectedIncrease+=0.5*m*increment*increment;
     momentum[i]+=m*velocity;center[i]+=fraction*r.getParticleCenter()[p][i];referenceCenter[i]+=fraction*r.getField<fields::mpm::particleReferencePosition>()[p][i];
     force[i]+=r.getParticleSurfaceTraction()[p][i];referenceForce[i]+=r.getField<fields::mpm::particleReferenceSurfaceTraction>()[p][i];
     near(r.getParticleCenter()[p][i]+r.getParticleSurfacePosition()[p][i],1.3+2*i);
     near(1.2*r.getField<fields::mpm::particleReferenceSurfacePosition>()[p][i],r.getParticleSurfacePosition()[p][i]);
     near(r.getField<fields::mpm::particleCohesiveForce>()[p][i],0);
     for(int a=0;a<3;++a)near(r.getField<fields::mpm::particleDeformationGradient>()[p][a][i],a==i?1.2:0);
   }
   near(r.getField<fields::mpm::particleKineticEnergy>()[p],0.5*m*speed2);totalKE+=0.5*m*speed2;
   assert(r.getField<fields::mpm::particleCopyFlag>()[p]==-1);assert(r.getField<fields::mpm::particleSubdivideFlag>()[p]==0);
 }
 near(mass,7);near(volume,13);near(v0,17);assert(static_cast<int>(ids.size())==r.size());assert(ids.count(100)==1);
 for(int i=0;i<3;++i){near(momentum[i],7*(2-i));near(center[i],1+i);near(referenceCenter[i],4+i);near(force[i],3+i);near(referenceForce[i],6+i);}
 near(totalKE,17.5+(affineVelocity?expectedIncrease:0));
}
int main()
{
 int cases=0;
 for(int shape:{0,1,2})for(int dimensions:{3})for(int correction:{0,1})
 {
   ParticleManager pm;auto & r=addRegion(pm,1);r.getParticleGeometryShape()[0]=shape;addRegion(pm,0);
   Solver solver;solver.initialize(2);solver.m_numDims=dimensions;solver.m_planeStrain=dimensions==2;solver.m_applyVelocityGradientWhenSplittingParticles=correction;
   assert(solver.subdivideParticles(pm,1));checkTransfer(r,shape,correction);assert(r.neighbors.count==0);assert(solver.m_mappedNodes[0].size()==0);++cases;
 }
 // Native XY storage is genuinely 2x2; three-wide framework vectors stay planar.
 for(int shape:{3,4}) for(int correction:{0,1})
 {
   ParticleManager pm;auto & r=addRegion(pm,1,shape==4?ParticleType::CPTI:ParticleType::CPDI,2);
   r.getParticleGeometryShape()[0]=shape;
   for(auto name:{"particleCenter","particleVelocity","particleSurfacePosition","particleSurfaceTraction"})
     r.get<Array<double,2>>(name)[0][2]=0;
   auto L=r.getField<fields::mpm::particleVelocityGradient>();
   for(int a=0;a<3;++a)L[0][2][a]=L[0][a][2]=0;
   r.getField<fields::mpm::particleDeformationGradient>()[0][2][2]=1;
   Solver solver;solver.initialize(1);solver.m_numDims=2;solver.m_planeStrain=1;
   solver.m_applyVelocityGradientWhenSplittingParticles=correction;
   assert(solver.subdivideParticles(pm,1));
   double mass=0,volume=0,area=0,center[2]={},momentum[2]={};
   auto G=r.getParticleGeometryCurrent();
   assert(G.size(1)==2 && G.size(2)==2);
   for(int p=0;p<r.size();++p)
   {
     double m=r.getField<fields::mpm::particleMass>()[p];
     mass+=m;volume+=r.getParticleVolume()[p];
     area+=particleGeometry::affine::Domain::fromGenerators(static_cast<particleGeometry::Shape>(shape),G[p]).measure();
     for(int d=0;d<2;++d){center[d]+=m*r.getParticleCenter()[p][d];momentum[d]+=m*r.getParticleVelocity()[p][d];}
     near(r.getParticleCenter()[p][2],0);near(r.getParticleVelocity()[p][2],0);
     near(r.getField<fields::mpm::particleDeformationGradient>()[p][2][2],1);
   }
   near(mass,7);near(volume,13);near(area,2.3*.6/(shape==4?2:1));
   near(center[0],7);near(center[1],14);near(momentum[0],14);near(momentum[1],7);++cases;
 }
 // No split, active CZ, damaged CZ, delete flag, and inactive ghost exclusions.
 for(int reason=0;reason<5;++reason)
 {
   ParticleManager pm;auto & r=addRegion(pm,1);Solver solver;solver.initialize(1);
   if(reason==0)solver.m_subdivideParticles=0;
   if(reason==1)r.getField<fields::mpm::particleCohesiveZoneFlag>()[0]=1;
   if(reason==2)r.getParticleSurfaceFlag()[0]=mpm::toInteger(mpm::SurfaceFlag::DamagedCohesive);
   if(reason==3)r.getField<fields::mpm::particleDeleteFlag>()[0]=1;
   if(reason==4){r.get<Array<int,1>>("particleRank")[0]=-1;r.setActiveParticleIndices();}
   assert(!solver.subdivideParticles(pm,1));assert(r.size()==1);++cases;
 }
 // Native explicit splitting is independent of both box VP-Hencky options.
 for(bool scaling:{false,true})for(bool reset:{false,true})
 {
   ParticleManager pm;auto & r=addRegion(pm,1);Solver solver;solver.initialize(1);solver.m_cpdiDomainScaling=scaling;
   solver.m_cpdiDomainScalingType=mpm::CPDIDomainScalingTypeOption::VPHencky;
   solver.m_domainResetType=reset?mpm::DomainResetTypeOption::VPHencky:mpm::DomainResetTypeOption::Isotropic;
   assert(solver.subdivideParticles(pm,1));checkTransfer(r,2,1);++cases;
 }
 // Automatic mode selects CPDI only; oversized split independently selects native geometry.
 for(int oversized:{0,1})
 {
   ParticleManager pm;auto & r=addRegion(pm,1);Solver solver;solver.initialize(1);solver.m_subdivideParticlesAutomatic=1;
   if(oversized){solver.m_subdivideParticles=0;solver.m_oversizedParticleTreatment=mpm::OversizedParticleTreatmentOption::Split;solver.m_defGradResetMaxParticleDomainToGridCellRatio=0.5;}
   assert(solver.subdivideParticles(pm,1)==(oversized!=0));if(oversized)checkTransfer(r,2,1);++cases;
 }
 // CPTI uses native geometry selection, even with box VP-Hencky controls enabled.
 for(int policy:{0,1,2})
 {
   ParticleManager pm;auto & r=addRegion(pm,1,ParticleType::CPTI);Solver solver;solver.initialize(1);
   solver.m_cpdiDomainScaling=1;solver.m_cpdiDomainScalingType=mpm::CPDIDomainScalingTypeOption::VPHencky;
   solver.m_domainResetType=mpm::DomainResetTypeOption::VPHencky;
   if(policy==1)solver.m_subdivideParticlesAutomatic=1;
   if(policy==2){solver.m_subdivideParticles=0;solver.m_oversizedParticleTreatment=mpm::OversizedParticleTreatmentOption::Split;solver.m_defGradResetMaxParticleDomainToGridCellRatio=.5;}
   assert(solver.subdivideParticles(pm,1)==(policy!=1));
   if(policy!=1){assert(r.getParticleType()==ParticleType::CPTI);checkTransfer(r,2,1);}
   ++cases;
 }
 // Existing CPDI volume-feasibility branch and rank-count reset fallback.
 {
   ParticleManager pm;auto & r=addRegion(pm,1,ParticleType::CPDI);r.getParticleGeometryShape()[0]=0;
   Solver solver;solver.initialize(1);solver.m_cpdiDomainScaling=1;solver.m_cpdiDomainScalingType=mpm::CPDIDomainScalingTypeOption::VPHencky;solver.m_subdivideParticlesAutomatic=1;
   assert(solver.subdivideParticles(pm,1));checkTransfer(r,0,1);++cases;
 }
 {
   ParticleManager pm;auto & r=addRegion(pm,1);Solver solver;solver.initialize(1);solver.m_subdivideParticles=0;
   solver.m_oversizedParticleTreatment=mpm::OversizedParticleTreatmentOption::Split;solver.m_defGradResetMaxParticleDomainToGridCellRatio=0.5;solver.m_oversizedParticleResetRankParticleCountThreshold=0;
   assert(!solver.subdivideParticles(pm,1));assert(r.size()==1);++cases;
 }
 // Mixed geometries in ONE subregion: no assumption of a subregion-wide shape.
 {
   ParticleManager pm;auto & r=addRegion(pm,3);Solver solver;solver.initialize(1);
   for(int p=0;p<3;++p)r.getParticleGeometryShape()[p]=p;
   assert(solver.subdivideParticles(pm,3));double sums[3]={};std::set<globalIndex> ids;
   for(int p=0;p<r.size();++p){int shape=r.getParticleGeometryShape()[p];sums[shape]+=r.getField<fields::mpm::particleMass>()[p];near(r.get<Array<double,1>>("importedScalar")[p],123+shape);ids.insert(r.getParticleID()[p]);}
   for(int shape=0;shape<3;++shape){near(sums[shape],7+shape);}
   assert(static_cast<int>(ids.size())==r.size());++cases;
 }
 // Simulated MPI collectives: an empty rank participates in a remote split;
 // a populated rank allocates IDs after a remote rank's reserved range.
 {
   ParticleManager pm;addRegion(pm,0);Solver solver;solver.initialize(1);MpiWrapper::counts={0,3};MpiWrapper::rank=0;MpiWrapper::globalMaximum=500;
   assert(solver.subdivideParticles(pm,0));assert(pm.regions[0]->size()==0);++cases;
   MpiWrapper::counts.clear();MpiWrapper::globalMaximum=0;
 }
 {
   ParticleManager pm;auto & r=addRegion(pm,1);Solver solver;solver.initialize(1);
   auto d=particleGeometry::affine::Domain::fromGenerators(particleGeometry::Shape::Tetrahedron,r.getParticleGeometryCurrent()[0]);
   int extra=particleGeometry::subdivision::sizePlan(d,3,0.99998).count-1;
   MpiWrapper::counts={3,extra};MpiWrapper::rank=1;MpiWrapper::globalMaximum=500;
   assert(solver.subdivideParticles(pm,1));for(int p=1;p<r.size();++p)assert(r.getParticleID()[p]==503+p);++cases;
   MpiWrapper::counts.clear();MpiWrapper::rank=0;MpiWrapper::globalMaximum=0;
 }
 std::cout<<cases<<" production splitter cases passed with host wrapper/collective adapters.\n";
}
