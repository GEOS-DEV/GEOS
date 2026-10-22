"""Exercise the production cell-quadrature callback and generated launch scripts.

Framework adapters isolate the ProblemManager callback from GEOS dependencies.
The recursive FE dispatch and solver capability methods are extracted unchanged
from production source. This does not run GEOS, MPI, Slurm, or Flux.
"""
from pathlib import Path
import ast
import os
import subprocess
import tempfile

root = next(p for p in Path(__file__).resolve().parents if (p / 'src/coreComponents').exists())
core = root / 'src/coreComponents'
pfw = Path(__file__).resolve().parents[1]


def body(source, signature):
    start = source.index(signature)
    opening = source.index('{', start)
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[opening:end]


problem = (core / 'mainInterface/ProblemManager.cpp').read_text()
callback = body(problem, 'elemRegion.forElementSubRegions< CellElementSubRegion >( [&]( auto & subRegion )')
base = (core / 'physicsSolvers/PhysicsSolverBase.hpp').read_text()
mpm = (core / 'physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp').read_text()
base_method = body(base, 'virtual bool requiresCellShapeFunctionGradients()')
mpm_method = body(mpm, 'bool requiresCellShapeFunctionGradients()')
dispatch = (core / 'finiteElement/FiniteElementDispatch.hpp').read_text()
dispatch = dispatch[dispatch.index('template< typename ... FE_TYPES >'):dispatch.index('template< typename LAMBDA >\nvoid\ndispatchlowOrder3D')]

program = r'''
#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
using localIndex = int;
#define GEOS_ERROR(...) throw std::runtime_error("unsupported FE dispatch")
#define GEOS_UNUSED_VAR(x) (void)(x)
#define GEOS_UNUSED_PARAM(x)
#define TYPEOFREF(x) std::remove_reference_t<decltype(x)>
namespace LvArray { namespace system { template<class T> std::string demangleType(T const *) { return "FE"; } } }
struct PhysicsSolverBase {
  virtual ~PhysicsSolverBase() = default;
  virtual bool requiresCellShapeFunctionGradients() const ''' + base_method + r'''
};
struct SolidMechanicsMPM: PhysicsSolverBase {
  bool requiresCellShapeFunctionGradients() const override ''' + mpm_method + r'''
};
int initialized = 0, gradients = 0;
namespace finiteElement {
struct FiniteElementBase {
  virtual ~FiniteElementBase() = default;
  virtual int getNumQuadraturePoints() const = 0;
  template<class> struct MeshData {};
  template<class FE, class SUB> static void initialize(int, int, int, SUB &, MeshData<SUB> &) { ++initialized; }
};
template<int Q> struct Element: FiniteElementBase {
  static constexpr int numQuadraturePoints = Q;
  int getNumQuadraturePoints() const override { return numQuadraturePoints; }
};
using Quad = Element<4>; using Triangle = Element<1>; using Hex = Element<8>;
''' + dispatch + r'''
}
#define ALL_FE_TYPES finiteElement::Hex
enum class ElementType { Quad, Triangle, Hex };
struct Discretization {
  std::unique_ptr<finiteElement::FiniteElementBase> factory(ElementType type) const {
    if(type == ElementType::Quad) return std::make_unique<finiteElement::Quad>();
    if(type == ElementType::Triangle) return std::make_unique<finiteElement::Triangle>();
    return std::make_unique<finiteElement::Hex>();
  }
  template<class SUB, class FE> void calculateShapeFunctionGradients(int, SUB *, typename FE::template MeshData<SUB>, FE &) const { ++gradients; }
};
enum class RestartFlags { NO_WRITE };
struct Wrapper {
  std::unique_ptr<finiteElement::FiniteElementBase> value;
  Wrapper &setRestartFlags(RestartFlags) { return *this; }
  finiteElement::FiniteElementBase &reference() { return *value; }
};
struct Subregion {
  ElementType type;
  Wrapper wrapper;
  bool excluded = false;
  ElementType getElementType() const { return type; }
  std::string getName() const { return "cells"; }
  template<class T> Wrapper &registerWrapper(std::string const &, std::unique_ptr<T> value) { wrapper.value = std::move(value); return wrapper; }
  void excludeWrappersFromPacking(std::initializer_list<std::string>) { excluded = true; }
};
struct MeshLevel { std::string getName() const { return "level"; } };
void check(PhysicsSolverBase const *solver, ElementType type, int expected, bool cellGradients, int previous = 0)
{
  std::string const discretizationName = "FE1", meshBodyName = "grid", regionName = "CellRegion1";
  int const nodeManager = 0, edgeManager = 0, faceManager = 0, X = 0;
  MeshLevel meshLevel;
  Subregion cells{type, {}, false};
  Discretization const discretization;
  Discretization const *feDiscretization = &discretization;
  std::map<std::tuple<std::string, std::string, std::string, std::string>, int> regionQuadrature;
  auto const key = std::make_tuple(meshBodyName, meshLevel.getName(), regionName, cells.getName());
  regionQuadrature[key] = previous;
  initialized = gradients = 0;
  auto initializeCell = [&](auto &subRegion)
''' + callback + r''';
  initializeCell(cells);
  assert(regionQuadrature.at(key) == std::max(previous, expected));
  assert(cells.wrapper.value && cells.excluded);
  assert(initialized == (cellGradients ? 1 : 0));
  assert(gradients == (cellGradients ? 1 : 0));
}
int main()
{
  PhysicsSolverBase fem;
  SolidMechanicsMPM mpmSolver;
  check(&mpmSolver, ElementType::Quad, 4, false);
  check(&mpmSolver, ElementType::Triangle, 1, false);
  check(&mpmSolver, ElementType::Hex, 8, false);
  check(&mpmSolver, ElementType::Quad, 4, false, 8);
  check(&fem, ElementType::Hex, 8, true);
  check(&fem, ElementType::Hex, 8, true, 14);
  bool unsupported = false;
  try { check(&fem, ElementType::Quad, 4, true); }
  catch(std::runtime_error const &) { unsupported = true; }
  assert(unsupported); // This fix does not add native-2D FEM support.
  std::cout << "Production startup callback: MPM quad/tri/hex, quadrature maximum, and unchanged FEM dispatch passed.\n";
}
'''

with tempfile.TemporaryDirectory(prefix='geos-mpm-startup-') as folder:
    folder = Path(folder)
    cpp, exe = folder / 'test.cpp', folder / 'test'
    cpp.write_text(program)
    subprocess.run(['c++', '-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Werror',
                    '-Wshadow', '-Wfloat-equal', '-pedantic', '-fsanitize=address,undefined',
                    '-fno-omit-frame-pointer', str(cpp), '-o', str(exe)], check=True)
    subprocess.run([str(exe)], check=True)

    # Render the actual PFW f-string templates, without importing/running PFW.
    source = (pfw / 'particleFileWriter.py').read_text()
    templates = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.JoinedStr):
            continue
        text = ast.get_source_segment(source, node.value)
        if 'srun --export=ALL,MV2_SMP_USE_CMA' in text or 'flux run --env=*' in text:
            templates.append(node.value)
    assert len(templates) == 2
    variables = dict(jobName='test', mNodes=1, mWallTimeMinutes=5, mWallTime='00:05:00',
                     mBank='test', mPartition='test', mCores=1, geosPath='/bin/true',
                     geosInputFileName='test.xml', xpar=1, ypar=1, zpar=1, restartStr='')
    for name in ['srun', 'flux']:
        command = folder / name
        command.write_text('#!/usr/bin/env bash\nexit "$MOCK_LAUNCH_STATUS"\n')
        command.chmod(0o755)
    env = dict(os.environ, PATH=str(folder) + os.pathsep + os.environ['PATH'])
    for template in templates:
        script = folder / 'launch.sh'
        script.write_text(eval(compile(ast.Expression(template), '<PFW launch template>', 'eval'), variables))
        subprocess.run(['bash', '-n', str(script)], check=True)
        for status in [0, 42]:
            result = subprocess.run(['bash', str(script)], env=dict(env, MOCK_LAUNCH_STATUS=str(status)),
                                    capture_output=True, text=True)
            assert result.returncode == status, result
            assert ('has completed' in result.stdout) == (status == 0), result
    print('Production Slurm/Flux launch templates: success/failure status and completion reporting passed (mock launchers).')
