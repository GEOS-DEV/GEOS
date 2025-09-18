/*
 * ------------------------------------------------------------------------------------------------------------
 * (Copy the appropriate license block from one of the other files)
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFIELDS_HPP_

#include "mesh/MeshFields.hpp" // Required for the DECLARE_FIELD macro

namespace geos
{
namespace fields
{

//================================================================================
// Fields shared by multiple gravity solvers
//================================================================================

DECLARE_FIELD( mediumDensity,
               "mediumDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( volumeIntegral,
               "volumeIntegral",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "VolumeIntegral of the cell." );

//================================================================================
// Fields specific to the GravityFE solver
//================================================================================

DECLARE_FIELD( volumeIntegral2d,
               "volumeIntegral2d",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "VolumeIntegral for adjoint computation." );

DECLARE_FIELD( adjoint,
               "adjoint",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Adjoint field." );

//================================================================================
// Fields specific to the GravityFE_CompositionalMultiphaseFVM solver
//================================================================================

DECLARE_FIELD( fluidDensity,
               "fluidDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Fluid density of the cell" );

DECLARE_FIELD( rockDensity,
               "rockDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Rock density of the cell" );

DECLARE_FIELD( porosity,
               "porosity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Porosity of the cell" );

} // namespace fields
} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFIELDS_HPP_
