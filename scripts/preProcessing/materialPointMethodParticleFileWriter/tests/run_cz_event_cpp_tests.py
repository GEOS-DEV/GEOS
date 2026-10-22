"""Compile the production CZ event class, scheduler excerpt and event handler.

Repository adapters replace GEOS registration/serialization and region storage.
This exercises event sequencing; it does not execute a full GEOS restart/run.
"""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--handler-source', type=Path, help='Optional pre-fix solver source for demonstrating the regression')
args = parser.parse_args()
tests = Path(__file__).resolve().parent
root = next(p for p in tests.parents if (p/'src/coreComponents').is_dir())
source = (args.handler_source or root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()


def block(text, start):
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"', lambda match: ' '*len(match[0]), text)
    begin = masked.index('{', start)
    end, depth = begin+1, 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return text[begin:end]


start = source.index('void SolidMechanicsMPM::triggerEvents(')
schedule = source.index('    // 1. Do not run completed events.', start)
schedule_end = source.index('    // Event specific handlers', schedule)
handler = source.index('    if( event.getCatalogName() == "ReferenceCohesiveZones" )', schedule_end)
program = '''#include "ReferenceCohesiveZonesMPMEvent.hpp"
#include <iostream>
namespace geos {
struct Runtime {
  bool m_axisymmetric=false;
  int m_eventReporting=0, m_computeParticleSurfaceNormalsAndPositions=0;
  mpm::NormalsAndPositionsMethodOption m_normalAndPositionMethod=mpm::NormalsAndPositionsMethodOption::LogisticRegression;
  EventManager eventManager;
  CohesiveZoneManager czManager;
  struct groupKeyStruct { static string cohesiveZoneManagerString() { return "CZ"; } };
  template<class T> T & getGroup(string const &) { return czManager; }
  void trigger(MPMEventBase & event, real64 time_n, real64 dt) {
'''+source[schedule:schedule_end]+'''
    (void)timeInterval;
'''+block(source, handler)+'''\n  }
};
}
'''+(tests/'mpm_events_cpp/test_reference_event.cpp').read_text()

with tempfile.TemporaryDirectory(prefix='geos-cz-event-') as directory:
    tmp = Path(directory)
    events = root/'src/coreComponents/events/mpmEvents'
    for name in ('ReferenceCohesiveZonesMPMEvent.hpp', 'ReferenceCohesiveZonesMPMEvent.cpp'):
        shutil.copy2(events/name, tmp/name)
    shutil.copy2(tests/'mpm_events_cpp/reference_event_test_support.hpp', tmp/'MPMEventBase.hpp')
    for name in ('physicsSolvers/solidMechanics/MPMSolverEnums.hpp', 'constitutive/cohesiveZone/CohesiveZoneBase.hpp'):
        stub = tmp/name
        stub.parent.mkdir(parents=True, exist_ok=True)
        stub.write_text('#pragma once\n')
    cpp = tmp/'test.cpp'
    cpp.write_text(program)
    executable = tmp/'test'
    subprocess.run(['c++', '-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Werror', '-Wshadow',
                    '-Wfloat-equal', '-pedantic', '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                    '-I'+str(tmp), str(cpp), str(tmp/'ReferenceCohesiveZonesMPMEvent.cpp'), '-o', str(executable)], check=True)
    env = dict(os.environ)
    env.setdefault('ASAN_OPTIONS', 'detect_leaks=0')
    subprocess.run([str(executable)], check=True, env=env)
