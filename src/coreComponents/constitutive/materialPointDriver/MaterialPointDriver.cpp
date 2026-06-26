/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MaterialPointDriver.cpp
 */

#include "MaterialPointDriver.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/SolidUtilities.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "LvArray/src/tensorOps.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace geos
{
namespace constitutive
{
namespace materialPointDriver
{

namespace
{

constexpr real64 smallNorm = 1.0e-30;

localIndex const voigtToI[6] = { 0, 1, 2, 1, 0, 0 };
localIndex const voigtToJ[6] = { 0, 1, 2, 2, 2, 1 };
localIndex const tensorToVoigt[3][3] = { { 0, 5, 4 }, { 5, 1, 3 }, { 4, 3, 2 } };

void throwError( string const & message )
{
  throw std::runtime_error( message );
}

real64 dot( Vector3 const & a, Vector3 const & b )
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector3 cross( Vector3 const & a, Vector3 const & b )
{
  return { { a[1] * b[2] - a[2] * b[1],
             a[2] * b[0] - a[0] * b[2],
             a[0] * b[1] - a[1] * b[0] } };
}

real64 norm( Vector3 const & a )
{
  return std::sqrt( dot( a, a ) );
}

Vector3 normalize( Vector3 const & a, string const & context )
{
  real64 const n = norm( a );
  if( n <= smallNorm )
  {
    throwError( "Cannot normalize near-zero vector while updating material frame: " + context );
  }
  return { { a[0] / n, a[1] / n, a[2] / n } };
}

Vector3 row( Matrix3 const & A, localIndex const i )
{
  return { { A[i][0], A[i][1], A[i][2] } };
}

void setRow( Matrix3 & A, localIndex const i, Vector3 const & a )
{
  A[i][0] = a[0];
  A[i][1] = a[1];
  A[i][2] = a[2];
}

Matrix3 identity()
{
  Matrix3 A{};
  for( localIndex i = 0; i < 3; ++i )
  {
    A[i][i] = 1.0;
  }
  return A;
}

Matrix3 multiply( Matrix3 const & A, Matrix3 const & B )
{
  Matrix3 C{};
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      for( localIndex k = 0; k < 3; ++k )
      {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return C;
}

Vector3 matVec( Matrix3 const & A, Vector3 const & x )
{
  return { { A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2],
             A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2],
             A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2] } };
}

real64 determinant( Matrix3 const & A )
{
  return A[0][0] * ( A[1][1] * A[2][2] - A[1][2] * A[2][1] )
       - A[0][1] * ( A[1][0] * A[2][2] - A[1][2] * A[2][0] )
       + A[0][2] * ( A[1][0] * A[2][1] - A[1][1] * A[2][0] );
}

Matrix3 cofactor( Matrix3 const & A )
{
  Matrix3 C{};
  C[0][0] =  A[1][1] * A[2][2] - A[1][2] * A[2][1];
  C[0][1] = -A[1][0] * A[2][2] + A[1][2] * A[2][0];
  C[0][2] =  A[1][0] * A[2][1] - A[1][1] * A[2][0];
  C[1][0] = -A[0][1] * A[2][2] + A[0][2] * A[2][1];
  C[1][1] =  A[0][0] * A[2][2] - A[0][2] * A[2][0];
  C[1][2] = -A[0][0] * A[2][1] + A[0][1] * A[2][0];
  C[2][0] =  A[0][1] * A[1][2] - A[0][2] * A[1][1];
  C[2][1] = -A[0][0] * A[1][2] + A[0][2] * A[1][0];
  C[2][2] =  A[0][0] * A[1][1] - A[0][1] * A[1][0];
  return C;
}

void copyToCArray( Matrix3 const & source, real64 ( & dest )[3][3] )
{
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      dest[i][j] = source[i][j];
    }
  }
}

void copyFromCArray( real64 const ( & source )[3][3], Matrix3 & dest )
{
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      dest[i][j] = source[i][j];
    }
  }
}

void computePolarRotation( Matrix3 const & F, Matrix3 & R )
{
  real64 f[3][3] = {};
  real64 r[3][3] = {};
  copyToCArray( F, f );
  LvArray::tensorOps::polarDecomposition< 3 >( r, f );
  copyFromCArray( r, R );
}

real64 voigtNorm( std::vector< real64 > const & values )
{
  real64 sum = 0.0;
  for( real64 const v : values )
  {
    sum += v * v;
  }
  return std::sqrt( sum );
}

void setSymmetricVelocityComponent( Matrix3 & L,
                                    localIndex const component,
                                    real64 const engineeringStrainRate )
{
  localIndex const i = voigtToI[component];
  localIndex const j = voigtToJ[component];
  if( i == j )
  {
    L[i][j] = engineeringStrainRate;
  }
  else
  {
    L[i][j] = 0.5 * engineeringStrainRate;
    L[j][i] = 0.5 * engineeringStrainRate;
  }
}

std::vector< real64 > solveLinearSystem( std::vector< std::vector< real64 > > A,
                                         std::vector< real64 > b )
{
  localIndex const n = static_cast< localIndex >( b.size() );
  for( localIndex col = 0; col < n; ++col )
  {
    localIndex pivot = col;
    real64 pivotAbs = std::abs( A[col][col] );
    for( localIndex rowIndex = col + 1; rowIndex < n; ++rowIndex )
    {
      real64 const valueAbs = std::abs( A[rowIndex][col] );
      if( valueAbs > pivotAbs )
      {
        pivotAbs = valueAbs;
        pivot = rowIndex;
      }
    }

    if( pivotAbs <= std::numeric_limits< real64 >::epsilon() )
    {
      throwError( "Singular stress-control Jacobian in material-point driver" );
    }

    if( pivot != col )
    {
      std::swap( A[pivot], A[col] );
      std::swap( b[pivot], b[col] );
    }

    real64 const diag = A[col][col];
    for( localIndex j = col; j < n; ++j )
    {
      A[col][j] /= diag;
    }
    b[col] /= diag;

    for( localIndex rowIndex = 0; rowIndex < n; ++rowIndex )
    {
      if( rowIndex == col )
      {
        continue;
      }
      real64 const factor = A[rowIndex][col];
      for( localIndex j = col; j < n; ++j )
      {
        A[rowIndex][j] -= factor * A[col][j];
      }
      b[rowIndex] -= factor * b[col];
    }
  }
  return b;
}

std::vector< string > mutableRealArrayNames()
{
  return { "density",
           "wavespeed",
           "velocityGradient",
           "deformationGradient",
           "materialDirection",
           "temperature",
           "temperatureRate",
           "jacobian",
           "lengthScale",
           "strengthScale",
           "damage",
           "basalPlaneDamage",
           "comminutionDamage",
           "plasticStrain",
           "relaxation",
           "plasticWork",
           "basalPlanePlasticWork",
           "effectiveBulkModulus",
           "effectiveShearModulus",
           "specificInternalEnergy",
           "supplementalPressure",
           "crackTipStressConcentration",
           "distanceToCrackTip",
           "crackTipDistance",
           "volume" };
}

} // namespace

Matrix3 MaterialFrameController::identityFrame()
{
  return identity();
}

Matrix3 MaterialFrameController::fromNormal( Vector3 const & normal,
                                             Vector3 tangentHint )
{
  Vector3 const n = normalize( normal, "input normal" );

  tangentHint = { { tangentHint[0] - dot( tangentHint, n ) * n[0],
                    tangentHint[1] - dot( tangentHint, n ) * n[1],
                    tangentHint[2] - dot( tangentHint, n ) * n[2] } };

  if( norm( tangentHint ) <= 1.0e-12 )
  {
    Vector3 fallback = { { 1.0, 0.0, 0.0 } };
    if( std::abs( dot( fallback, n ) ) > 0.9 )
    {
      fallback = { { 0.0, 1.0, 0.0 } };
    }
    tangentHint = { { fallback[0] - dot( fallback, n ) * n[0],
                      fallback[1] - dot( fallback, n ) * n[1],
                      fallback[2] - dot( fallback, n ) * n[2] } };
  }

  Vector3 const t1 = normalize( tangentHint, "input tangentHint" );
  Vector3 const t2 = normalize( cross( n, t1 ), "constructed tangent" );

  Matrix3 frame{};
  setRow( frame, 0, n );
  setRow( frame, 1, t1 );
  setRow( frame, 2, t2 );
  return frame;
}

Matrix3 MaterialFrameController::normalizedFrame( Matrix3 const & frame )
{
  Vector3 const e0 = normalize( row( frame, 0 ), "row 0" );
  Vector3 t1 = row( frame, 1 );
  t1 = { { t1[0] - dot( t1, e0 ) * e0[0],
           t1[1] - dot( t1, e0 ) * e0[1],
           t1[2] - dot( t1, e0 ) * e0[2] } };
  Vector3 const e1 = normalize( t1, "row 1" );
  Vector3 const e2 = normalize( cross( e0, e1 ), "row 2" );

  Matrix3 result{};
  setRow( result, 0, e0 );
  setRow( result, 1, e1 );
  setRow( result, 2, e2 );
  return result;
}

Matrix3 MaterialFrameController::update( Matrix3 const & referenceFrame,
                                         Matrix3 const & deformationGradient,
                                         Matrix3 const & polarRotation,
                                         MaterialFrameUpdateMode const mode )
{
  if( mode == MaterialFrameUpdateMode::fixed )
  {
    return referenceFrame;
  }

  Matrix3 result{};
  Matrix3 const cofF = cofactor( deformationGradient );

  if( mode == MaterialFrameUpdateMode::rotation )
  {
    for( localIndex i = 0; i < 3; ++i )
    {
      setRow( result, i, normalize( matVec( polarRotation, row( referenceFrame, i ) ), "rotation update" ) );
    }
    return result;
  }

  if( mode == MaterialFrameUpdateMode::fiber || mode == MaterialFrameUpdateMode::mpmCofactor )
  {
    Matrix3 const & map = ( mode == MaterialFrameUpdateMode::fiber ) ? deformationGradient : cofF;
    for( localIndex i = 0; i < 3; ++i )
    {
      setRow( result, i, normalize( matVec( map, row( referenceFrame, i ) ), "basis transport" ) );
    }
    return result;
  }

  Vector3 const n = normalize( matVec( cofF, row( referenceFrame, 0 ) ), "normal/cofactor update" );
  Vector3 tangentSeed = row( referenceFrame, 1 );
  if( mode == MaterialFrameUpdateMode::graphite )
  {
    tangentSeed = matVec( deformationGradient, tangentSeed );
  }
  tangentSeed = { { tangentSeed[0] - dot( tangentSeed, n ) * n[0],
                    tangentSeed[1] - dot( tangentSeed, n ) * n[1],
                    tangentSeed[2] - dot( tangentSeed, n ) * n[2] } };
  if( norm( tangentSeed ) <= 1.0e-12 )
  {
    tangentSeed = row( referenceFrame, 2 );
    if( mode == MaterialFrameUpdateMode::graphite )
    {
      tangentSeed = matVec( deformationGradient, tangentSeed );
    }
    tangentSeed = { { tangentSeed[0] - dot( tangentSeed, n ) * n[0],
                      tangentSeed[1] - dot( tangentSeed, n ) * n[1],
                      tangentSeed[2] - dot( tangentSeed, n ) * n[2] } };
  }

  Vector3 const t1 = normalize( tangentSeed, "normal/graphite tangent update" );
  Vector3 const t2 = normalize( cross( n, t1 ), "normal/graphite cross product" );
  setRow( result, 0, n );
  setRow( result, 1, t1 );
  setRow( result, 2, t2 );
  return result;
}

real64 EnergyIntegrator::stressPower( Voigt6 const & stress,
                                      Matrix3 const & velocityGradient )
{
  real64 power = 0.0;
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      power += stress[tensorToVoigt[i][j]] * velocityGradient[i][j];
    }
  }
  return power;
}

real64 EnergyIntegrator::integrateStressPowerTrapezoid( MaterialPointState & state,
                                                        Voigt6 const & stressBeginning,
                                                        Voigt6 const & stressEnd,
                                                        EnergyOptions const & options )
{
  if( options.mode != EnergyMode::stressPower )
  {
    return 0.0;
  }

  real64 const rho = std::max( state.density, smallNorm );
  real64 const p0 = stressPower( stressBeginning, state.velocityGradient );
  real64 const p1 = stressPower( stressEnd, state.velocityGradient );
  real64 const mechanicalIncrement = 0.5 * state.dt * ( p0 + p1 ) / rho;

  state.accumulatedStressPower += mechanicalIncrement;
  state.specificInternalEnergy += options.retentionFactor * mechanicalIncrement;

  if( options.temperatureMode == TemperatureMode::adiabatic )
  {
    real64 const oldTemperature = state.temperature;
    real64 const heatCapacity = std::max( options.heatCapacity, smallNorm );
    state.temperature += options.retentionFactor * mechanicalIncrement / heatCapacity;
    state.temperatureRate = ( state.temperature - oldTemperature ) / std::max( state.dt, smallNorm );
  }
  else if( options.temperatureMode == TemperatureMode::isothermal ||
           options.temperatureMode == TemperatureMode::prescribed )
  {
    state.temperatureRate = 0.0;
    state.temperature = state.prescribedTemperature;
  }

  return mechanicalIncrement;
}

void ConstitutivePointSnapshot::capture( ContinuumBase & model )
{
  m_realArrays.clear();
  m_hasStress = model.getStress().size( 0 ) > 0 && model.getStress().size( 1 ) > 0;
  if( m_hasStress )
  {
    auto const stress = model.getStress();
    auto const oldStress = model.getOldStress();
    for( localIndex i = 0; i < 6; ++i )
    {
      m_stress[i] = stress[0][0][i];
      m_oldStress[i] = oldStress[0][0][i];
    }
  }

  for( string const & name : mutableRealArrayNames() )
  {
    if( !model.hasWrapper( name ) )
    {
      continue;
    }

    dataRepository::WrapperBase const & wrapper = model.getWrapperBase( name );
    int const rank = wrapper.numArrayDims();
    if( rank < 1 || rank > 3 )
    {
      continue;
    }

    RealArraySnapshot snapshot;
    snapshot.name = name;
    snapshot.rank = rank;

    if( rank == 1 )
    {
      array1d< real64 > const & array = model.getReference< array1d< real64 > >( name );
      snapshot.dims[0] = array.size( 0 );
      snapshot.values.reserve( snapshot.dims[0] );
      for( localIndex i = 0; i < snapshot.dims[0]; ++i )
      {
        snapshot.values.push_back( array[i] );
      }
    }
    else if( rank == 2 )
    {
      array2d< real64 > const & array = model.getReference< array2d< real64 > >( name );
      snapshot.dims[0] = array.size( 0 );
      snapshot.dims[1] = array.size( 1 );
      snapshot.values.reserve( snapshot.dims[0] * snapshot.dims[1] );
      for( localIndex i = 0; i < snapshot.dims[0]; ++i )
      {
        for( localIndex j = 0; j < snapshot.dims[1]; ++j )
        {
          snapshot.values.push_back( array[i][j] );
        }
      }
    }
    else
    {
      array3d< real64 > const & array = model.getReference< array3d< real64 > >( name );
      snapshot.dims[0] = array.size( 0 );
      snapshot.dims[1] = array.size( 1 );
      snapshot.dims[2] = array.size( 2 );
      snapshot.values.reserve( snapshot.dims[0] * snapshot.dims[1] * snapshot.dims[2] );
      for( localIndex i = 0; i < snapshot.dims[0]; ++i )
      {
        for( localIndex j = 0; j < snapshot.dims[1]; ++j )
        {
          for( localIndex k = 0; k < snapshot.dims[2]; ++k )
          {
            snapshot.values.push_back( array[i][j][k] );
          }
        }
      }
    }

    m_realArrays.emplace_back( std::move( snapshot ) );
  }
}

void ConstitutivePointSnapshot::restore( ContinuumBase & model ) const
{
  if( m_hasStress )
  {
    auto stress = model.getStress();
    auto oldStress = model.getOldStress();
    for( localIndex i = 0; i < 6; ++i )
    {
      stress[0][0][i] = m_stress[i];
      oldStress[0][0][i] = m_oldStress[i];
    }
  }

  for( RealArraySnapshot const & snapshot : m_realArrays )
  {
    if( !model.hasWrapper( snapshot.name ) )
    {
      continue;
    }

    localIndex cursor = 0;
    if( snapshot.rank == 1 )
    {
      array1d< real64 > & array = model.getReference< array1d< real64 > >( snapshot.name );
      for( localIndex i = 0; i < std::min( snapshot.dims[0], array.size( 0 ) ); ++i )
      {
        array[i] = snapshot.values[cursor++];
      }
    }
    else if( snapshot.rank == 2 )
    {
      array2d< real64 > & array = model.getReference< array2d< real64 > >( snapshot.name );
      for( localIndex i = 0; i < std::min( snapshot.dims[0], array.size( 0 ) ); ++i )
      {
        for( localIndex j = 0; j < std::min( snapshot.dims[1], array.size( 1 ) ); ++j )
        {
          array[i][j] = snapshot.values[cursor++];
        }
      }
    }
    else if( snapshot.rank == 3 )
    {
      array3d< real64 > & array = model.getReference< array3d< real64 > >( snapshot.name );
      for( localIndex i = 0; i < std::min( snapshot.dims[0], array.size( 0 ) ); ++i )
      {
        for( localIndex j = 0; j < std::min( snapshot.dims[1], array.size( 1 ) ); ++j )
        {
          for( localIndex k = 0; k < std::min( snapshot.dims[2], array.size( 2 ) ); ++k )
          {
            array[i][j][k] = snapshot.values[cursor++];
          }
        }
      }
    }
  }
}

void MaterialPointDriver::initializeOnePointModel( ContinuumBase & model,
                                                   MaterialPointState & state )
{
  if( std::abs( determinant( state.deformationGradient ) ) <= smallNorm )
  {
    state.deformationGradient = identity();
  }
  if( std::abs( determinant( state.referenceMaterialFrame ) ) <= smallNorm )
  {
    state.referenceMaterialFrame = identity();
  }
  state.referenceMaterialFrame = MaterialFrameController::normalizedFrame( state.referenceMaterialFrame );
  state.oldDeformationGradient = state.deformationGradient;
  computePolarRotation( state.deformationGradient, state.rotationBeginning );
  state.rotationEnd = state.rotationBeginning;
  state.materialFrame = MaterialFrameController::update( state.referenceMaterialFrame,
                                                         state.deformationGradient,
                                                         state.rotationEnd,
                                                         state.materialFrameUpdate );
  state.jacobian = determinant( state.deformationGradient );
  state.volume = state.referenceVolume * state.jacobian;
  state.density = state.referenceDensity / std::max( state.jacobian, smallNorm );
  pushDependencies( model, state );

  if( model.getStress().size( 0 ) > 0 && model.getStress().size( 1 ) > 0 )
  {
    auto stress = model.getStress();
    auto oldStress = model.getOldStress();
    for( localIndex i = 0; i < 6; ++i )
    {
      stress[0][0][i] = state.stress[i];
      oldStress[0][0][i] = state.stress[i];
      state.oldStress[i] = state.stress[i];
    }
  }
}

void MaterialPointDriver::pushDependencies( ContinuumBase & model,
                                            MaterialPointState const & state )
{
  if( model.hasWrapper( "velocityGradient" ) )
  {
    array3d< real64 > & array = model.getReference< array3d< real64 > >( "velocityGradient" );
    for( localIndex i = 0; i < 3; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        array[0][i][j] = state.velocityGradient[i][j];
      }
    }
  }

  if( model.hasWrapper( "deformationGradient" ) )
  {
    array3d< real64 > & array = model.getReference< array3d< real64 > >( "deformationGradient" );
    for( localIndex i = 0; i < 3; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        array[0][i][j] = state.deformationGradient[i][j];
      }
    }
  }

  if( model.hasWrapper( "materialDirection" ) )
  {
    dataRepository::WrapperBase const & materialDirectionWrapper = model.getWrapperBase( "materialDirection" );
    int const rank = materialDirectionWrapper.numArrayDims();
    if( rank == 2 )
    {
      array2d< real64 > & array = model.getReference< array2d< real64 > >( "materialDirection" );
      for( localIndex i = 0; i < 3; ++i )
      {
        array[0][i] = state.materialFrame[0][i];
      }
    }
    else if( rank == 3 )
    {
      array3d< real64 > & array = model.getReference< array3d< real64 > >( "materialDirection" );
      for( localIndex i = 0; i < 3; ++i )
      {
        for( localIndex j = 0; j < 3; ++j )
        {
          array[0][i][j] = state.materialFrame[i][j];
        }
      }
    }
    else
    {
      throwError( "Unsupported materialDirection rank in material-point driver" );
    }
  }

  if( model.hasWrapper( "temperature" ) )
  {
    model.getReference< array1d< real64 > >( "temperature" )[0] = state.temperature;
  }
  if( model.hasWrapper( "temperatureRate" ) )
  {
    model.getReference< array1d< real64 > >( "temperatureRate" )[0] = state.temperatureRate;
  }
  if( model.hasWrapper( "specificInternalEnergy" ) )
  {
    model.getReference< array1d< real64 > >( "specificInternalEnergy" )[0] = state.specificInternalEnergy;
  }
  if( model.hasWrapper( "lengthScale" ) )
  {
    model.getReference< array1d< real64 > >( "lengthScale" )[0] = state.lengthScale;
  }
  if( model.hasWrapper( "strengthScale" ) )
  {
    model.getReference< array1d< real64 > >( "strengthScale" )[0] = state.strengthScale;
  }
  if( model.hasWrapper( "jacobian" ) )
  {
    model.getReference< array2d< real64 > >( "jacobian" )[0][0] = state.jacobian;
  }
  if( model.hasWrapper( "density" ) )
  {
    model.getReference< array2d< real64 > >( "density" )[0][0] = state.density;
  }
  if( model.hasWrapper( "volume" ) )
  {
    model.getReference< array1d< real64 > >( "volume" )[0] = state.volume;
  }
}

void MaterialPointDriver::pullOutputs( ContinuumBase & model,
                                       MaterialPointState & state )
{
  if( model.getStress().size( 0 ) > 0 && model.getStress().size( 1 ) > 0 )
  {
    auto const stress = model.getStress();
    auto const oldStress = model.getOldStress();
    for( localIndex i = 0; i < 6; ++i )
    {
      state.stress[i] = stress[0][0][i];
      state.oldStress[i] = oldStress[0][0][i];
    }
  }

  if( model.hasWrapper( "temperature" ) )
  {
    state.temperature = model.getReference< array1d< real64 > >( "temperature" )[0];
  }
  if( model.hasWrapper( "temperatureRate" ) )
  {
    state.temperatureRate = model.getReference< array1d< real64 > >( "temperatureRate" )[0];
  }
  if( model.hasWrapper( "specificInternalEnergy" ) )
  {
    state.specificInternalEnergy = model.getReference< array1d< real64 > >( "specificInternalEnergy" )[0];
  }
  if( model.hasWrapper( "density" ) )
  {
    state.density = model.getReference< array2d< real64 > >( "density" )[0][0];
  }
}

std::vector< localIndex > MaterialPointDriver::stressControlledComponents( StepPlan const & plan )
{
  std::vector< localIndex > result;
  for( ComponentControl const & control : plan.controls )
  {
    if( control.mode == ControlMode::stress )
    {
      result.push_back( control.component );
    }
  }
  return result;
}

void MaterialPointDriver::buildVelocityGradient( StepPlan const & plan,
                                                 std::vector< real64 > const & unknownStrainIncrements,
                                                 Matrix3 & velocityGradient )
{
  velocityGradient = Matrix3{};
  localIndex stressControlIndex = 0;
  for( ComponentControl const & control : plan.controls )
  {
    if( control.component < 0 || control.component >= 6 )
    {
      throwError( "Invalid Voigt component in material-point driver control" );
    }

    if( control.mode == ControlMode::strainIncrement )
    {
      setSymmetricVelocityComponent( velocityGradient, control.component, control.value / plan.dt );
    }
    else if( control.mode == ControlMode::strainRate )
    {
      setSymmetricVelocityComponent( velocityGradient, control.component, control.value );
    }
    else if( control.mode == ControlMode::trueStrainRate )
    {
      if( control.component > 2 )
      {
        throwError( "trueStrainRate is only defined for diagonal components" );
      }
      setSymmetricVelocityComponent( velocityGradient, control.component, control.value );
    }
    else if( control.mode == ControlMode::stress )
    {
      setSymmetricVelocityComponent( velocityGradient,
                                     control.component,
                                     unknownStrainIncrements[stressControlIndex] / plan.dt );
      ++stressControlIndex;
    }
  }
}

void MaterialPointDriver::updateKinematicsForTrial( MaterialPointState & state,
                                                    real64 const dt )
{
  state.dt = dt;
  state.oldDeformationGradient = state.deformationGradient;
  Matrix3 const LTimesF = multiply( state.velocityGradient, state.oldDeformationGradient );
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      state.deformationGradient[i][j] = state.oldDeformationGradient[i][j] + dt * LTimesF[i][j];
      state.fDot[i][j] = LTimesF[i][j];
    }
  }

  computePolarRotation( state.oldDeformationGradient, state.rotationBeginning );
  computePolarRotation( state.deformationGradient, state.rotationEnd );

  state.materialFrame = MaterialFrameController::update( state.referenceMaterialFrame,
                                                         state.deformationGradient,
                                                         state.rotationEnd,
                                                         state.materialFrameUpdate );

  state.strainIncrement[0] = state.velocityGradient[0][0] * dt;
  state.strainIncrement[1] = state.velocityGradient[1][1] * dt;
  state.strainIncrement[2] = state.velocityGradient[2][2] * dt;
  state.strainIncrement[3] = ( state.velocityGradient[1][2] + state.velocityGradient[2][1] ) * dt;
  state.strainIncrement[4] = ( state.velocityGradient[0][2] + state.velocityGradient[2][0] ) * dt;
  state.strainIncrement[5] = ( state.velocityGradient[0][1] + state.velocityGradient[1][0] ) * dt;

  for( localIndex i = 0; i < 6; ++i )
  {
    state.accumulatedStrain[i] += state.strainIncrement[i];
  }

  state.jacobian = determinant( state.deformationGradient );
  state.volume = state.referenceVolume * state.jacobian;
  state.density = state.referenceDensity / std::max( state.jacobian, smallNorm );
}

void MaterialPointDriver::callConstitutiveUpdate( ContinuumBase & model,
                                                  MaterialPointState & state,
                                                  real64 const dt,
                                                  bool const commit )
{
  pushDependencies( model, state );

  if( model.getOldStress().size( 0 ) > 0 && model.getOldStress().size( 1 ) > 0 )
  {
    auto oldStress = model.getOldStress();
    for( localIndex i = 0; i < 6; ++i )
    {
      oldStress[0][0][i] = state.oldStress[i];
    }
  }

  bool const hyperelasticUpdate = model.getCatalogName() == "HyperelasticMMS" ||
                                  model.getCatalogName() == "Hyperelastic" ||
                                  model.getCatalogName() == "Chiumenti";

  ConstitutivePassThruMPM< ContinuumBase >::execute( model, [&] ( auto & castedConstitutiveModel )
  {
    auto constitutiveWrapper = castedConstitutiveModel.createKernelUpdates();

    real64 stress[6] = {};
    if( hyperelasticUpdate )
    {
      real64 FminusI[3][3] = {};
      copyToCArray( state.deformationGradient, FminusI );
      for( localIndex i = 0; i < 3; ++i )
      {
        FminusI[i][i] -= 1.0;
      }
      constitutiveWrapper.hyperUpdate( 0, 0, FminusI, stress );
    }
    else
    {
      real64 strainIncrement[6] = {};
      real64 rotationBeginning[3][3] = {};
      real64 rotationEnd[3][3] = {};
      for( localIndex i = 0; i < 6; ++i )
      {
        strainIncrement[i] = state.strainIncrement[i];
      }
      copyToCArray( state.rotationBeginning, rotationBeginning );
      copyToCArray( state.rotationEnd, rotationEnd );

      SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                              0,
                                              0,
                                              dt,
                                              strainIncrement,
                                              rotationBeginning,
                                              rotationEnd,
                                              stress );
    }

    for( localIndex i = 0; i < 6; ++i )
    {
      state.stress[i] = stress[i];
    }
    if( commit )
    {
      constitutiveWrapper.saveConvergedState( 0, 0 );
    }
  } );

  pullOutputs( model, state );
}

std::vector< real64 > MaterialPointDriver::residualForUnknowns( ContinuumBase & model,
                                                                MaterialPointState const & stateBeginning,
                                                                ConstitutivePointSnapshot const & modelSnapshot,
                                                                StepPlan const & plan,
                                                                std::vector< real64 > const & unknownStrainIncrements,
                                                                bool const commit,
                                                                MaterialPointState * trialStateOut )
{
  modelSnapshot.restore( model );

  MaterialPointState trialState = stateBeginning;
  trialState.dt = plan.dt;
  buildVelocityGradient( plan, unknownStrainIncrements, trialState.velocityGradient );
  updateKinematicsForTrial( trialState, plan.dt );
  callConstitutiveUpdate( model, trialState, plan.dt, commit );

  std::vector< real64 > residual;
  residual.reserve( unknownStrainIncrements.size() );
  for( ComponentControl const & control : plan.controls )
  {
    if( control.mode == ControlMode::stress )
    {
      residual.push_back( trialState.stress[control.component] - control.target );
    }
  }

  if( trialStateOut != nullptr )
  {
    *trialStateOut = trialState;
  }
  return residual;
}

StepDiagnostics MaterialPointDriver::advanceOneStep( ContinuumBase & model,
                                                     MaterialPointState & state,
                                                     StepPlan const & plan )
{
  if( plan.dt <= 0.0 )
  {
    throwError( "Material-point driver requires positive dt" );
  }

  StepDiagnostics diagnostics;
  diagnostics.converged = false;

  MaterialPointState const stateBeginning = state;
  ConstitutivePointSnapshot modelSnapshot;
  modelSnapshot.capture( model );

  std::vector< localIndex > const stressComponents = stressControlledComponents( plan );
  std::vector< real64 > unknowns( stressComponents.size(), 0.0 );

  if( stressComponents.empty() )
  {
    MaterialPointState accepted;
    residualForUnknowns( model, stateBeginning, modelSnapshot, plan, unknowns, true, &accepted );
    diagnostics.energyIncrement = EnergyIntegrator::integrateStressPowerTrapezoid( accepted,
                                                                                   stateBeginning.stress,
                                                                                   accepted.stress,
                                                                                   plan.energyOptions );
    accepted.time = stateBeginning.time + plan.dt;
    state = accepted;
    pushDependencies( model, state );
    diagnostics.converged = true;
    diagnostics.iterations = 1;
    diagnostics.finalResidualNorm = 0.0;
    return diagnostics;
  }

  real64 referenceResidualNorm = 0.0;
  std::vector< real64 > residual;

  for( localIndex iteration = 0; iteration < plan.solveOptions.maxIterations; ++iteration )
  {
    residual = residualForUnknowns( model, stateBeginning, modelSnapshot, plan, unknowns, false, nullptr );
    real64 const residualNorm = voigtNorm( residual );
    if( iteration == 0 )
    {
      referenceResidualNorm = std::max( residualNorm, smallNorm );
      diagnostics.initialResidualNorm = residualNorm;
    }

    diagnostics.iterations = iteration + 1;
    diagnostics.finalResidualNorm = residualNorm;

    real64 const tolerance = plan.solveOptions.absoluteTolerance +
                             plan.solveOptions.relativeTolerance * referenceResidualNorm;
    if( residualNorm <= tolerance )
    {
      diagnostics.converged = true;
      break;
    }

    localIndex const n = static_cast< localIndex >( unknowns.size() );
    std::vector< std::vector< real64 > > jacobian( n, std::vector< real64 >( n, 0.0 ) );

    for( localIndex j = 0; j < n; ++j )
    {
      std::vector< real64 > perturbedUnknowns = unknowns;
      real64 const h = std::max( plan.solveOptions.minimumFiniteDifference,
                                 plan.solveOptions.finiteDifferenceScale * std::max( 1.0, std::abs( unknowns[j] ) ) );
      perturbedUnknowns[j] += h;
      std::vector< real64 > const perturbedResidual =
        residualForUnknowns( model, stateBeginning, modelSnapshot, plan, perturbedUnknowns, false, nullptr );

      for( localIndex i = 0; i < n; ++i )
      {
        jacobian[i][j] = ( perturbedResidual[i] - residual[i] ) / h;
      }
    }

    std::vector< real64 > rhs = residual;
    for( real64 & value : rhs )
    {
      value = -value;
    }
    std::vector< real64 > correction = solveLinearSystem( jacobian, rhs );

    real64 const correctionNorm = voigtNorm( correction );
    if( correctionNorm > plan.solveOptions.maximumCorrectionNorm )
    {
      real64 const scale = plan.solveOptions.maximumCorrectionNorm / correctionNorm;
      for( real64 & value : correction )
      {
        value *= scale;
      }
    }

    real64 bestNorm = std::numeric_limits< real64 >::max();
    std::vector< real64 > bestUnknowns = unknowns;
    bool acceptedCorrection = false;
    localIndex const maxLineSearchCuts = plan.solveOptions.useLineSearch ? 8 : 1;
    real64 alpha = 1.0;
    for( localIndex cut = 0; cut < maxLineSearchCuts; ++cut )
    {
      std::vector< real64 > candidate = unknowns;
      for( localIndex j = 0; j < n; ++j )
      {
        candidate[j] += alpha * correction[j];
      }

      std::vector< real64 > const candidateResidual =
        residualForUnknowns( model, stateBeginning, modelSnapshot, plan, candidate, false, nullptr );
      real64 const candidateNorm = voigtNorm( candidateResidual );

      if( candidateNorm < bestNorm )
      {
        bestNorm = candidateNorm;
        bestUnknowns = candidate;
      }
      if( candidateNorm < residualNorm )
      {
        acceptedCorrection = true;
        break;
      }
      alpha *= 0.5;
    }

    unknowns = bestUnknowns;
    if( !acceptedCorrection && bestNorm >= residualNorm )
    {
      diagnostics.message = "Line search did not reduce residual";
    }
  }

  MaterialPointState accepted;
  std::vector< real64 > const finalResidual =
    residualForUnknowns( model, stateBeginning, modelSnapshot, plan, unknowns, true, &accepted );
  diagnostics.finalResidualNorm = voigtNorm( finalResidual );
  real64 const tolerance = plan.solveOptions.absoluteTolerance +
                           plan.solveOptions.relativeTolerance * std::max( diagnostics.initialResidualNorm, smallNorm );
  diagnostics.converged = diagnostics.converged || diagnostics.finalResidualNorm <= tolerance;

  if( !diagnostics.converged && diagnostics.message.empty() )
  {
    std::ostringstream oss;
    oss << "Stress-control solve failed after " << diagnostics.iterations
        << " iterations; final residual norm = " << diagnostics.finalResidualNorm;
    diagnostics.message = oss.str();
  }

  diagnostics.energyIncrement = EnergyIntegrator::integrateStressPowerTrapezoid( accepted,
                                                                                 stateBeginning.stress,
                                                                                 accepted.stress,
                                                                                 plan.energyOptions );
  accepted.time = stateBeginning.time + plan.dt;
  state = accepted;
  pushDependencies( model, state );
  return diagnostics;
}

} // namespace materialPointDriver
} // namespace constitutive
} // namespace geos
