/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreRieszMFD.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPRERIESZMFD_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPRERIESZMFD_HPP_

#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geos
{

namespace hypre
{

/**
 * @brief Create the Riesz-map block preconditioner for the adaptive mixed MFD saddle point.
 *
 * The unknowns are the live fluxes q (faces of MFD cells), the pressures p_M of the MFD cells and
 * the pressures p_T of the TPFA cells; the condensed fluxes are recovered from their closure rows.
 * The preconditioner is the Riesz map of the product norm H(div; MFD region) x L2(MFD cells) x
 * H1(TPFA cells),
 *
 *   P = diag( M + B_M^T W_M^{-1} B_M ,  W_M ,  C_T + S_TT ),
 *
 * with B_M the face-to-cell incidence of the live faces into the MFD cells, W_M the L2 mass of the
 * MFD cells (plus accumulation) and S_TT the two-point Laplacian of the TPFA cells. The flux block
 * is solved by CG under one ADS cycle, the TPFA block by one BoomerAMG cycle. The interface faces
 * carry no weight: they are bounded by the normal-trace inequality. All blocks are read off the
 * assembled system, whose row scalings are identified and removed; the solver provides the de Rham
 * sub-complex of the live faces, the dof markers, the MFD flag and the L2 scale of each cell.
 *
 * @param params the linear solver parameters (dof markers and auxiliary discretization data)
 * @param precond the output preconditioner wrapper
 */
void createRieszMFD( LinearSolverParameters const & params,
                     HyprePrecWrapper & precond );

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPRERIESZMFD_HPP_*/
