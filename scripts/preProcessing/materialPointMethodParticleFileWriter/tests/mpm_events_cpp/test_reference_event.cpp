using namespace geos;
static void require(bool value, char const *message)
{
  if(!value) throw std::runtime_error(message);
}
static void configure(ReferenceCohesiveZonesMPMEvent & event)
{
  const_cast<string_array &>(event.getRegionNames())={"bond","second"};
}
int main()
{
  // User run: 1 ns startup step followed by a 2.4975 us CFL step. The old
  // half-step window admits BOTH calls, and resets the reference twice.
  for(bool axisymmetric : {false,true})
  {
    Runtime runtime; runtime.m_axisymmetric=axisymmetric;
    ReferenceCohesiveZonesMPMEvent event("activateBond",nullptr); configure(event);
    auto const &wrapper=event.wrappers.at("referenceApplied");
    require(wrapper.input==dataRepository::InputFlags::FALSE,"Latch must not be XML input");
    require(wrapper.restart==dataRepository::RestartFlags::WRITE_AND_READ,"Latch must survive restart");
    runtime.trigger(event,0.,1e-9);
    for(auto & entry:runtime.czManager.regions) entry.second.initialized=1;
    runtime.trigger(event,1e-9,2.4975024975024975e-6);
    for(auto const & entry:runtime.czManager.regions)
      require(entry.second.referenceRequests==1 && entry.second.initialized==1,
              "BUG: growing dt requests a second cohesive reference");
    require(event.referenceApplied()==1 && !event.isComplete(),"Activation must retain event end handling");

    ReferenceCohesiveZonesMPMEvent resumed("activateBond",nullptr); configure(resumed);
    restore(resumed,save(event));
    runtime.trigger(resumed,1e-9,1e-4);
    require(runtime.czManager.regions.at("bond").referenceRequests==1,"Restart reactivated the reference");
    resumed.setEndTime(1e-3);
    runtime.trigger(resumed,1e-3,1e-6);
    require(resumed.isComplete() && runtime.czManager.regions.at("bond").enabled==0,
            "Finite end time must disable the region");
    runtime.trigger(resumed,1e-3,1.);
    require(runtime.czManager.regions.at("bond").referenceRequests==1,"Completed event reactivated");

    ReferenceCohesiveZonesMPMEvent later("newReference",nullptr); configure(later);
    later.setStartTime(.01);
    runtime.trigger(later,.009,1e-6);
    runtime.trigger(later,.01,1e-6);
    require(runtime.czManager.regions.at("bond").referenceRequests==2,
            "A distinct event must still be able to request a new reference");
  }
  // hasStarted is already true while dependent events wait out their delay.
  // Save and restore in that interval, then activate once at the scheduled time.
  Runtime dependent;
  ReferenceCohesiveZonesMPMEvent delayed("delayed",nullptr); configure(delayed);
  delayed.dependencies={"ready"}; delayed.delay=.001; delayed.duration=.002;
  dependent.trigger(delayed,0.,1e-9);
  require(!delayed.hasStarted() && !delayed.referenceApplied(),"Unmet dependency must wait");
  dependent.eventManager.dependency.setIsComplete(1);
  dependent.trigger(delayed,.002,1e-6);
  require(delayed.hasStarted() && !delayed.referenceApplied(),"Scheduling must not count as activation");
  ReferenceCohesiveZonesMPMEvent pendingRestart("delayed",nullptr); configure(pendingRestart);
  pendingRestart.dependencies={"ready"}; restore(pendingRestart,save(delayed));
  dependent.trigger(pendingRestart,.0025,1e-6);
  require(!pendingRestart.referenceApplied(),"Delayed activation occurred too soon");
  dependent.trigger(pendingRestart,.003,1e-9);
  dependent.trigger(pendingRestart,.003+1e-9,1e-4);
  require(dependent.czManager.regions.at("bond").referenceRequests==1,"Delayed reference did not activate exactly once");
  dependent.trigger(pendingRestart,.005,1e-6);
  require(pendingRestart.isComplete(),"Dependent event did not finish");

  // First active step may fall after the old narrow start window when dt
  // shrinks across the scheduled time. Activation must still occur.
  Runtime skipped;
  ReferenceCohesiveZonesMPMEvent missed("missed",nullptr); configure(missed);
  missed.setStartTime(.001);
  skipped.trigger(missed,.0009,.00015);
  require(!missed.referenceApplied(),"Inactive event started early");
  skipped.trigger(missed,.00105,1e-6);
  require(missed.referenceApplied() && skipped.czManager.regions.at("bond").referenceRequests==1,
          "Shrinking dt skipped activation");
  std::cout<<"CZ event: adaptive dt, delayed/dependent activation, restart-state persistence, finite end and distinct events passed.\n";
}
