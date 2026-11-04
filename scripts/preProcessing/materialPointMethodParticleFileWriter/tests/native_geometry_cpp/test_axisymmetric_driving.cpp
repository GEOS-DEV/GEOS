void close(double a,double b,double tol=1.e-11) {if(std::abs(a-b)>tol*std::max({1.,std::abs(a),std::abs(b)}))throw std::runtime_error("numeric mismatch");}
template<class F> void rejects(F f){try{f();}catch(std::runtime_error const&){return;}throw std::runtime_error("missing validation");}
ParticleManager material(double volume,double rr,double yy,double hoop) {
 ParticleManager p;ParticleSubRegion r;r.volume={volume};r.indices={0};r.stress=Table{{rr,yy,hoop,0,0,0}};p.regions.push_back(r);return p;
}
int main() {
 double const pi=std::acos(-1.);
 // The same exact finite wall map works for linear, cosine and quintic
 // interpolation, crossing knots, holding the final value, and continuation.
 for(auto type:{mpm::InterpolationOption::Linear,mpm::InterpolationOption::Cosine,mpm::InterpolationOption::Smoothstep}) {
  Driving d;d.m_fTableInterpType=type;d.validateAxisymmetricBoundaryControls();
  auto p=material(pi,0,0,0);SpatialPartition partition;NodeManager nodes;
  double t=0;
  for(double dt:{.13,.21,.19,.3,.31,.2}) {
   double const oldF=d.m_domainF[0];d.updatePrescribedKinematicsForExplicitStep(dt,t,p,partition);
   double const frac=std::min(t+dt,1.);
   double const blend=type==mpm::InterpolationOption::Linear?frac:(type==mpm::InterpolationOption::Cosine?.5*(1-std::cos(pi*frac)):10*std::pow(frac,3)-15*std::pow(frac,4)+6*std::pow(frac,5));
   close(d.m_domainF[0],1-.1*blend);close(d.m_domainF[1],1-.2*blend);
   close(d.m_domainIncrementalStretch[0],d.m_domainF[0]/oldF);
   d.resizeGrid(partition,nodes,dt);
   close(nodes.positions[1][0],d.m_domainF[0]);close(nodes.positions[1][1],d.m_domainF[1]);close(nodes.positions[1][2],0);
   close(d.physicalDomainVolume(),pi*std::pow(d.m_domainF[0],2)*d.m_domainF[1]);
   close(partition.ratio[0],d.m_domainF[0]);close(d.m_hEl[2],1);
   close(d.m_xGlobalMin[0],0);t+=dt;
  }
 }
 // Geometry versus material averages: annular box and unequal ring volumes.
 Driving a;a.m_xGlobalMin[0]=2;a.m_xGlobalMax[0]=3;a.m_domainExtent[0]=1;
 auto p=material(pi,10,20,99);auto second=material(4*pi,30,40,199);p.regions.push_back(second.regions[0]);
 Array<double> mean(3);double volume=0;a.computeBoxMetrics(p,mean,volume);
 close(mean[0],26);close(mean[1],36);close(mean[2],179);close(volume,5*pi);
 p.regions.pop_back();a.computeBoxMetrics(p,mean,volume);close(mean[0],2); // void is included in denominator
 // Legacy rr feedback, full-ring density, and correct annular absolute-radius CFL cap.
 a.m_prescribedBoundaryFTable=0;a.m_stressControl[0]=1;a.m_domainStress[0]=-1e10;
 SpatialPartition part;a.stressControl(.1,p,part);close(a.m_domainL[0],-.2*.1/.1/3);
 Driving d;d.m_prescribedBoundaryFTable=0;d.m_stressControl[0]=1;d.m_domainStress[0]=-100;
 p=material(pi,0,0,12345);d.stressControl(.1,p,part);close(d.m_domainL[0],-.1);close(d.m_domainL[1],0);close(d.m_domainL[2],0);
 // Stress feedback is sigma_rr, independent of the hoop value.
 p.regions[0].stress[0][2]=-999;d.stressControl(.1,p,part);close(d.m_domainL[0],-.1);
 // Closed-loop elastic surrogate includes physical radial/hoop coupling.
 // This tests controller signs and convergence, not full solver dynamics.
 double rr=0,ey=0;
 for(int n=0;n<200;++n) {
  p=material(pi,rr,ey,rr);d.stressControl(.01,p,part);
  rr+=1600*d.m_domainL[0]*.01;ey+=800*d.m_domainL[0]*.01;
 }
 close(rr,-100,1e-8);
 // Mixed Fyy prescription and rr control advance one grid map consistently.
 Driving mixed;mixed.m_stressControl[0]=1;mixed.m_stressTable=Table{{0,-100,0,0}};
 p=material(pi,0,0,0);mixed.updatePrescribedKinematicsForExplicitStep(.1,0,p,part);
 close(mixed.m_domainF[0],.99);close(mixed.m_domainF[1],.98);close(mixed.m_domainIncrementalStretch[0],.99);
 // A one-sided cap can stop a prescribed compression without moving the grid.
 Driving limited;limited.m_stressControl[0]=2;limited.m_stressTable=Table{{0,-100,0,0}};
 p=material(pi,-200,0,0);limited.updatePrescribedKinematicsForExplicitStep(.1,0,p,part);
 close(limited.m_domainL[0],0);close(limited.m_domainF[0],1);close(limited.m_domainIncrementalStretch[0],1);
 // Adaptive seek and active feedback run the complete production controller.
 Driving seek;seek.m_prescribedBoundaryFTable=0;seek.m_stressControl[0]=1;seek.m_stressControlAdaptiveSeek=1;
 seek.m_stressControlMinStrainRate=1;seek.m_stressControlMaxSeekStrainIncrement=0;
 seek.m_xGlobalMin[0]=9;seek.m_xGlobalMax[0]=10;seek.m_domainExtent[0]=1;
 seek.m_domainStress[0]=-100;p=material(19*pi,0,0,0);seek.stressControl(.1,p,part);
 assert(seek.m_domainL[0]<0);assert(std::abs(seek.m_domainL[0])<=.020000000001);close(seek.m_domainL[2],0);
 seek.m_stressControlAdaptiveInitialized=0;p=material(19*pi,-90,0,-90);seek.stressControl(.1,p,part);
 assert(seek.m_domainL[0]<0);assert(std::abs(seek.m_domainL[0])<=.020000000001);
 // A BC table supplies initial types and may release walls after the F ramp.
 Driving release;release.m_prescribedBcTable=1;release.m_boundaryConditionTypes[1]=0;
 release.m_bcTable=Table{{0,1,2,1,2},{1,1,0,1,0}};release.validateAxisymmetricBoundaryControls();
 // Axis preservation at BC-table switches and validation of future rows.
 Driving bc;bc.m_prescribedBcTable=1;bc.m_bcTable=Table{{0,1,2,1,2},{1,1,3,1,2}};
 bc.validateAxisymmetricBoundaryControls();bc.boundaryConditionUpdate(.1,1);assert(bc.m_boundaryConditionTypes[1]==3);assert(bc.m_boundaryConditionTypes[0]==1);
 bc.m_bcTable[1][1]=0;rejects([&]{bc.validateAxisymmetricBoundaryControls();});
 Driving invalid;invalid.m_boundaryConditionTypes[1]=1;rejects([&]{invalid.validateAxisymmetricBoundaryControls();});
 invalid=Driving{};invalid.m_updateMethod=mpm::UpdateMethodOption::FMPM;invalid.m_boundaryConditionTypes[1]=3;rejects([&]{invalid.validateAxisymmetricBoundaryControls();});
 invalid=Driving{};invalid.m_enablePrescribedBoundaryTransverseVelocities[1]=1;invalid.m_prescribedBoundaryTransverseVelocities[1][1]=1;rejects([&]{invalid.validateAxisymmetricBoundaryControls();});
 invalid=Driving{};invalid.m_enablePrescribedBoundaryTransverseVelocities[3]=1;invalid.m_prescribedBoundaryTransverseVelocities[3][1]=1;rejects([&]{invalid.validateAxisymmetricBoundaryControls();});
 // Actual per-node wall assembly: axisymmetric friction never adds slip energy,
 // has no impulse in separation, and scales reaction with full-ring mass.
 for(int normal=0;normal<2;++normal)for(int side=0;side<2;++side)for(double mu:{0.,.2,5.}) {
  Wall w;w.dir0=normal;w.dir1=(normal+1)%3;w.dir2=(normal+2)%3;w.face=2*normal+side;w.positiveNormal=side;
  int const tangent=1-normal;double const sign=side?1.:-1.;w.gridVelocity[0][0][normal]=sign;w.gridVelocity[0][0][tangent]=2;
  w.gridCenterOfVolume[0][0][normal]=sign*.01;w.boundaryFaceFrictionCoefficients[w.face]=mu;w.apply();
  close(w.targetVelocity[0][0][normal],0);close(w.faceReaction,-sign*100);
  if(mu>0) {close(w.targetVelocity[0][0][tangent],2-std::min(2.,mu));close(w.targetAcceleration[0][0][tangent],-std::min(2.,mu)/.1);}
  else assert(w.targetCount[0][0][tangent]==0);
 }
 Wall away;away.gridVelocity[0][0][0]=-1;away.apply();close(away.faceReaction,0);assert(away.targetCount[0][0][1]==0);
 Wall gap;gap.gridCenterOfVolume[0][0][0]=-1;gap.apply();close(gap.faceReaction,0);assert(gap.targetCount[0][0][1]==0);
 Wall moving;moving.domainL[0]=-.1;moving.domainIncrementalStretch[0]=.99;moving.gridVelocity[0][0][0]=0;moving.apply();close(moving.targetVelocity[0][0][0],-.099);assert(moving.faceReaction<0);
 Wall prescribed;prescribed.boundaryType=2;prescribed.domainL[0]=-.1;prescribed.domainIncrementalStretch[0]=.9;prescribed.apply();close(prescribed.targetVelocity[0][0][0],-.09);
 std::cout<<"PASS: production table/grid maps, ring stress metrics, PID/adaptive/one-sided controls, axis validation, moving/contact walls\n";
}
