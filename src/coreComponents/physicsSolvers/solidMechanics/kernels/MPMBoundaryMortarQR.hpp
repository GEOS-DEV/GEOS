/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MPM_BOUNDARY_MORTAR_QR_HPP
#define GEOS_MPM_BOUNDARY_MORTAR_QR_HPP

#include "MPMBoundaryMortar.hpp"

namespace geos
{
namespace mpm
{
namespace boundaryMortar
{

/** CPU QR backends; all retain the original projection and acceptance checks.
 *
 * Reference recomputes every remaining column norm at each pivot. Native uses
 * safeguarded partial-norm updates, with exact refresh before rank termination.
 * Lapack uses DGEQP3 and DORMQR from an optional, LP64 BLAS/LAPACK build.
 * There is no normal-equation solve, compliance, or implicit fallback to GS.
 */
enum class QRBackend { Reference, Native, Lapack };

bool lapackQRAvailable();
QRBackend selectQRBackend( std::string const & name );
char const * qrBackendName( QRBackend backend );
MinimumNormSolver minimumNormSolver( QRBackend backend );

QRResult minimumNormNative( std::vector< double > & bt, Index n, Index m,
                            std::vector< double > const & b, double rankTolerance );
QRResult minimumNormLapack( std::vector< double > & bt, Index n, Index m,
                            std::vector< double > const & b, double rankTolerance );

} // namespace boundaryMortar
} // namespace mpm
} // namespace geos
#endif
