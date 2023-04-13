/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsFractures.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSFRACTURES_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSFRACTURES_HPP_

namespace geos
{

namespace poromechanicsFracturesKernels
{

/**
 * @brief A struct to perform volume, aperture and fracture traction updates
 */
struct StateUpdateKernel
{

  /**
   * @brief Launch the kernel function doing volume, aperture and fracture traction updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] contactWrapper the wrapper implementing the contact relationship
   * @param[in] dispJump the displacement jump
   * @param[in] pressure the pressure
   * @param[in] area the area
   * @param[in] volume the volume
   * @param[out] deltaVolume the change in volume
   * @param[out] aperture the aperture
   * @param[in] minimumHydraulicAperture the
   * @param[in] oldHydraulicAperture the old hydraulic aperture
   * @param[out] hydraulicAperture the effecture aperture
   * @param[in] fractureTraction the fracture traction
   */
  template< typename POLICY, typename POROUS_WRAPPER >
  static void
  launch( localIndex const size,
          POROUS_WRAPPER const & porousMaterialWrapper,
          arrayView2d< real64 const > const & dispJump,
          arrayView1d< real64 const > const & pressure,
          arrayView1d< real64 const > const & area,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 > const & deltaVolume,
          arrayView1d< real64 > const & aperture,
          arrayView1d< real64 const > const & minimumHydraulicAperture,
          arrayView1d< real64 const > const & oldHydraulicAperture,
          arrayView1d< real64 > const & hydraulicAperture,
          arrayView2d< real64 const > const & fractureTraction )
  {

    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      // update aperture to be equal to the normal displacement jump
      aperture[k] = dispJump[k][0]; // the first component of the jump is the normal one.

      hydraulicAperture[k] = minimumHydraulicAperture[k] + aperture[k];

      deltaVolume[k] = hydraulicAperture[k] * area[k] - volume[k];

      real64 const jump[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( dispJump[k] );
      real64 const traction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( fractureTraction[k] );

      porousMaterialWrapper.updateStateFromPressureApertureJumpAndTraction( k, 0, pressure[k],
                                                                            oldHydraulicAperture[k], hydraulicAperture[k],
                                                                            jump, traction );

    } );
  }

};

} // namespace poromechanicsFracturesKernels

} /* namespace geos */

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSFRACTURES_HPP_
