# -*- coding: utf-8 -*-

# Copy this file to userDefs_<LC username>.py, or run setupMPM to generate it.
# Values left as CHANGEME have not been configured on that machine yet.
import os
import platform

_CHANGEME = 'CHANGEME'
_machines = ('dane', 'tuolumne')

_node = platform.node().lower()
node = platform.node()
dane = 'dane' in _node
lassen = 'lassen' in _node
tuolumne = ('tuolumne' in _node) or _node.startswith('tuo')
rzhound = 'rzhound' in _node
tioga = 'tioga' in _node

# Optional escape hatch for login nodes or unusual node names:
#   export GEOS_MPM_MACHINE=dane
#   export GEOS_MPM_MACHINE=tuolumne
_machine = os.environ.get('GEOS_MPM_MACHINE', '').strip().lower()
if not _machine:
  if dane:
    _machine = 'dane'
  elif tuolumne:
    _machine = 'tuolumne'

if _machine not in _machines:
  raise RuntimeError(
    "Could not identify Dane or Tuolumne from host '%s'. "
    "Set GEOS_MPM_MACHINE=dane or GEOS_MPM_MACHINE=tuolumne."
    % node
  )

_pfwPath = {machine: _CHANGEME for machine in _machines}
_geosPath = {machine: _CHANGEME for machine in _machines}
_testRunDirectory = {machine: _CHANGEME for machine in _machines}
_defaultRunDirectory = {machine: _CHANGEME for machine in _machines}
_defaultBank = {machine: _CHANGEME for machine in _machines}
_defaultPythonCommand = {machine: _CHANGEME for machine in _machines}
_defaultScheduler = {machine: _CHANGEME for machine in _machines}
_defaultSlurmPartition = {machine: _CHANGEME for machine in _machines}
_defaultFluxQueue = {machine: _CHANGEME for machine in _machines}
_defaultDeviceRanksPerNode = {machine: _CHANGEME for machine in _machines}
_defaultMpiRanksPerNode = {machine: _CHANGEME for machine in _machines}

# Fill entries by running setupMPM on each machine, or edit them manually.
_pfwPath["dane"] = _CHANGEME
_pfwPath["tuolumne"] = _CHANGEME

_geosPath["dane"] = _CHANGEME
_geosPath["tuolumne"] = _CHANGEME

_testRunDirectory["dane"] = _CHANGEME
_testRunDirectory["tuolumne"] = _CHANGEME

_defaultRunDirectory["dane"] = _CHANGEME
_defaultRunDirectory["tuolumne"] = _CHANGEME

_defaultBank["dane"] = _CHANGEME
_defaultBank["tuolumne"] = _CHANGEME

_defaultPythonCommand["dane"] = _CHANGEME
_defaultPythonCommand["tuolumne"] = _CHANGEME

_defaultScheduler["dane"] = _CHANGEME
_defaultScheduler["tuolumne"] = _CHANGEME

_defaultSlurmPartition["dane"] = _CHANGEME
_defaultSlurmPartition["tuolumne"] = _CHANGEME

_defaultFluxQueue["dane"] = _CHANGEME
_defaultFluxQueue["tuolumne"] = _CHANGEME

_defaultDeviceRanksPerNode["dane"] = _CHANGEME
_defaultDeviceRanksPerNode["tuolumne"] = _CHANGEME

_defaultMpiRanksPerNode["dane"] = _CHANGEME
_defaultMpiRanksPerNode["tuolumne"] = _CHANGEME

# Public names consumed by the particle-file-writer scripts.
pfwPath=_pfwPath[_machine]
geosPath=_geosPath[_machine]
testRunDirectory=_testRunDirectory[_machine]
defaultRunDirectory=_defaultRunDirectory[_machine]
defaultBank=_defaultBank[_machine]
defaultPythonCommand=_defaultPythonCommand[_machine]
pythonCommand=defaultPythonCommand
defaultPython=defaultPythonCommand  # backward-compatible alias

# Convenience defaults. Existing PFW inputs can ignore these.
defaultScheduler=_defaultScheduler[_machine]
defaultSlurmPartition=_defaultSlurmPartition[_machine]
defaultFluxQueue=_defaultFluxQueue[_machine]
defaultDeviceRanksPerNode=_defaultDeviceRanksPerNode[_machine]
defaultMpiRanksPerNode=_defaultMpiRanksPerNode[_machine]

_required = {
  'pfwPath': pfwPath,
  'geosPath': geosPath,
  'testRunDirectory': testRunDirectory,
  'defaultRunDirectory': defaultRunDirectory,
}
_missing = [name for name, value in _required.items() if value == _CHANGEME]
if _missing:
  raise RuntimeError(
    "The following userDefs values are still CHANGEME for %s: %s. "
    "Run setupMPM on that machine or edit this file."
    % (_machine, ', '.join(_missing))
  )
