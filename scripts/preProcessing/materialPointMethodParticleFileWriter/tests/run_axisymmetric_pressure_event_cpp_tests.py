"""Compile actual pressure validation/handler bodies with CPU-only adapters.

Tests active/completed gating, ramp endpoints and retained state. Restart
registration is checked statically; this does not perform GEOS restart I/O.
"""
from pathlib import Path
import re
import subprocess
import tempfile

root = Path(__file__).resolve().parents[4]
source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()

def block(text, marker):
    start = text.index(marker)
    brace = text.index('{', start)
    depth = 1
    end = brace + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[start:end]

def validation(name):
    text = (root/f'src/coreComponents/events/mpmEvents/{name}MPMEvent.cpp').read_text()
    method = block(text, f'void {name}MPMEvent::postInputInitialization()')
    return method.replace(f'{name}MPMEvent::postInputInitialization', 'validate'+name)

event_source = source[source.index('void SolidMechanicsMPM::triggerEvents('):]
start = event_source.index('    GEOS_ERROR_IF(m_axisymmetric &&')
end = event_source.index('    // 2. If this event depends', start)
guard = event_source[start:end]
start = event_source.index('    real64 const startTime = event.getStartTime();')
end = event_source.index('    // Event specific handlers', start)
active = event_source[start:end]
handlers = block(event_source, '    if( event.getCatalogName() == "BoreholePressure" )')
handlers += block(event_source, '    if( event.getCatalogName() == "ConfiningPressure" )')
interpolate = block(source, 'void SolidMechanicsMPM::interpolateValueInRange(').replace('SolidMechanicsMPM::', '')
completion = source[source.index('void SolidMechanicsMPM::checkEventCompletion('):]
completion = block(completion, '    if( event.isComplete() || !event.hasStarted() )') + block(completion, '    if( time_n > event.getEndTime() )')

for name in ('boreholeRadius', 'boreholeStress', 'enableBoreholePressure',
             'enableConfiningPressure', 'confiningStress', 'uniformConfiningPressure',
             'confiningPressureSphereRadius', 'confiningPressureBoxMin', 'confiningPressureBoxMax'):
    wrapper = re.search(r'registerWrapper\( "'+name+r'".*?;', source, re.S)
    assert wrapper and 'RestartFlags::WRITE_AND_READ' in wrapper.group(), name

program = r'''
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using real64=double;
#define GEOS_MARK_FUNCTION
#define GEOS_ERROR_IF(condition,message) do {if(condition) throw std::runtime_error("invalid input");} while(false)
#define GEOS_ERROR(message) throw std::runtime_error("invalid interpolation")
#define GEOS_LOG_RANK_0(message) do {} while(false)
#define GEOS_LOG_RANK_0_IF(condition,message) do {} while(false)
namespace LvArray { namespace math {using std::pow;} }
namespace mpm {enum class InterpolationOption {Linear,Cosine,Smoothstep};}
struct MPMEventBase {void postInputInitialization() {}}
;
struct PressureEvent:MPMEventBase {
  std::string catalog="BoreholePressure";
  double m_startTime=0,m_endTime=1,m_startPressure=0,m_endPressure=10,m_boreholeRadius=2;
  double m_confiningPressureSphereRadius=0;
  int m_uniformPressure=0,m_interpType=0;
  std::vector<double> m_confiningPressureBoxMin={-1,-10,0},m_confiningPressureBoxMax={3,10,0};
  bool complete=false,started=false;
  std::string getCatalogName() const{return catalog;}
  bool isComplete() const{return complete;} bool hasStarted() const{return started;}
  void setIsComplete(int v){complete=v!=0;} void setHasStarted(int v){started=v!=0;}
  double getStartTime() const{return m_startTime;} double getEndTime() const{return m_endTime;}
  double getTimeInterval() const{return m_endTime-m_startTime;}
  double getStartPressure() const{return m_startPressure;} double getEndPressure() const{return m_endPressure;}
  double getBoreholeRadius() const{return m_boreholeRadius;} int getInterpType() const{return m_interpType;}
  int getUniformPressure() const{return m_uniformPressure;}
  double getConfiningPressureSphereRadius() const{return m_confiningPressureSphereRadius;}
  auto const& getConfiningPressureBoxMin() const{return m_confiningPressureBoxMin;}
  auto const& getConfiningPressureBoxMax() const{return m_confiningPressureBoxMax;}
''' + validation('BoreholePressure') + validation('ConfiningPressure') + r'''
};
using BoreholePressureMPMEvent=PressureEvent;
using ConfiningPressureMPMEvent=PressureEvent;
template<class T> T dynamicCast(PressureEvent& event){return static_cast<T>(event);}
struct Harness {
  bool m_axisymmetric=true;
  int m_enableBoreholePressure=0,m_enableConfiningPressure=0,m_uniformConfiningPressure=0;
  double m_boreholeRadius=0,m_confiningPressureSphereRadius=0,m_boreholeStress[6]={},m_confiningStress[6]={};
  std::vector<double> m_confiningPressureBoxMin,m_confiningPressureBoxMax;
''' + interpolate + r'''
  void trigger(PressureEvent& event,double time_n,double dt) {
''' + guard + active + handlers + r'''
  }
  void complete(PressureEvent& event,double time_n) {
''' + completion + r'''
  }
};
bool near(double a,double b){return std::abs(a-b)<1e-12*std::max({1.,std::abs(a),std::abs(b)});}
template<class F> void rejected(F fn){bool fail=false;try{fn();}catch(std::runtime_error const&){fail=true;}assert(fail);}
int main() {
  // Actual event post-input validation, including omitted boxes in sphere/uniform modes.
  PressureEvent valid;valid.validateBoreholePressure();valid.validateConfiningPressure();
  for(double radius:{0.,-1.,std::numeric_limits<double>::infinity()}) {
    auto bad=valid;bad.m_boreholeRadius=radius;rejected([&]{bad.validateBoreholePressure();});
  }
  for(int interp:{-1,3}) {auto bad=valid;bad.m_interpType=interp;
    rejected([&]{bad.validateBoreholePressure();});rejected([&]{bad.validateConfiningPressure();});}
  auto bad=valid;bad.m_startPressure=std::numeric_limits<double>::quiet_NaN();
  rejected([&]{bad.validateBoreholePressure();});rejected([&]{bad.validateConfiningPressure();});
  bad=valid;bad.m_confiningPressureBoxMin={};rejected([&]{bad.validateConfiningPressure();});
  bad=valid;bad.m_confiningPressureBoxMax={1,2};rejected([&]{bad.validateConfiningPressure();});
  bad=valid;bad.m_confiningPressureBoxMax[0]=-2;rejected([&]{bad.validateConfiningPressure();});
  bad=valid;bad.m_uniformPressure=1;rejected([&]{bad.validateConfiningPressure();});
  bad=valid;bad.m_confiningPressureSphereRadius=std::numeric_limits<double>::quiet_NaN();
  rejected([&]{bad.validateConfiningPressure();});
  for(int uniform:{0,1}) {auto good=valid;good.m_uniformPressure=uniform;
    good.m_confiningPressureSphereRadius=uniform?0:3;
    good.m_confiningPressureBoxMin.clear();good.m_confiningPressureBoxMax.clear();good.validateConfiningPressure();}

  for(std::string name:{"BoreholePressure","ConfiningPressure"}) for(int interp=0;interp<3;++interp) {
    Harness h;PressureEvent e;e.catalog=name;e.m_interpType=interp;
    double const* stress=name=="BoreholePressure"?h.m_boreholeStress:h.m_confiningStress;
    h.trigger(e,-.2,.1);assert(!e.started&&near(stress[0],0));
    h.trigger(e,0,.1);assert(e.started&&near(stress[0],0));
    h.trigger(e,.25,.1);
    double expected=0;h.interpolateValueInRange(.25,0,1,0,10,expected,static_cast<mpm::InterpolationOption>(interp));
    assert(near(stress[0],-expected));
    h.trigger(e,.9,.15); // next step overshoots the end; persist the exact endpoint
    for(int d=0;d<6;++d)assert(near(stress[d],d<3?-10:0));
    h.complete(e,1.05);assert(e.complete);
    e.m_endPressure=99;h.trigger(e,1.05,.15);assert(near(stress[0],-10));
    h.trigger(e,0,.1);assert(near(stress[0],-10)); // completed events cannot reactivate
    e=PressureEvent{};e.catalog=name;e.m_startTime=.5;e.m_endTime=.5;
    h.trigger(e,.5,.1);assert(near(stress[0],-10)); // instantaneous pressure step
    h.m_axisymmetric=false;e=PressureEvent{};e.catalog=name;e.m_interpType=interp;
    h.trigger(e,.9,.15);h.interpolateValueInRange(.9,0,1,0,10,expected,static_cast<mpm::InterpolationOption>(interp));
    assert(near(stress[0],-expected)); // retain Cartesian time sampling
  }
  Harness h;PressureEvent inner,outer;outer.catalog="ConfiningPressure";
  outer.m_startTime=.5;outer.m_endTime=2;outer.m_endPressure=3;
  h.trigger(inner,.9,.15);h.complete(inner,1.05);h.trigger(outer,1.05,.15);
  assert(near(h.m_boreholeStress[0],-10));assert(h.m_confiningStress[0]<0&&h.m_confiningStress[0]>-3);
  h.trigger(outer,1.9,.15);assert(near(h.m_confiningStress[0],-3)&&near(h.m_boreholeStress[0],-10));
  std::cout<<"Production pressure input validation, event guards, ramp endpoints and held loads passed; restart wrappers verified.\n";
}
'''
with tempfile.TemporaryDirectory(prefix='geos-axis-pressure-events-') as directory:
    tmp=Path(directory);cpp=tmp/'test.cpp';cpp.write_text(program);exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer',str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
