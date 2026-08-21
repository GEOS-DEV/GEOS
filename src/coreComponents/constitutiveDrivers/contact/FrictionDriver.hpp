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

#ifndef GEOS_CONSTITUTIVEDRIVERS_CONTACT_FRICTIONDRIVER_HPP
#define GEOS_CONSTITUTIVEDRIVERS_CONTACT_FRICTIONDRIVER_HPP


#include "constitutiveDrivers/ConstitutiveDriver.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactSolverBase.hpp"
// #include "physicsSolvers/solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp"

namespace geos
{
namespace constitutive
{
class FrictionBase;
}

namespace LUsolver{

// ============================================================================
// Minimal dense linear algebra: LU with partial pivoting
// ============================================================================
 
// In-place LU decomposition of an n x n matrix stored row-major in A
// (A is overwritten: U in upper incl. diag, L multipliers below diag,
// unit diagonal for L implied). piv[i] = original row that ended up at row i.
// Returns false if the matrix is (numerically) singular.
bool luDecompose(std::vector<double>& A, int n, std::vector<int>& piv);

 
// Solve A x = b given the LU-factored A (from luDecompose) and its pivot.
std::vector<double> luSolve(const std::vector<double>& LU, int n,
                             const std::vector<int>& piv,
                             const std::vector<double>& b);
 
// Convenience: solve a SMALL fixed-size system (used for the 3x3 bubble
// condensation) via the same LU machinery, working on a flat copy.
template <int N>
std::array<std::array<double,N>,N> invertSmall(const std::array<std::array<double,N>,N>& M);

}

namespace mimicAssembly{

using IndexType = std::ptrdiff_t;
using CRSMat = CRSMatrix< real64, IndexType, IndexType >;

array2d<real64> buildD(double K, double G);

// ----------------------------------------------------------------------------
// Reference hex: trilinear shape functions and their parent-coordinate
// derivatives at (xi, eta, zeta) in [-1,1]^3. Standard node ordering.
// ----------------------------------------------------------------------------
constexpr std::array<std::array<double,3>,8> refCorners = {{
    {-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
    {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}
}};

void shapeAndGrad(double xi, double eta, double zeta,
                   std::array<double,8>& N,
                   std::array<std::array<double,3>,8>& dNdXi);
 
// Face bubble (Frigo Eq. 20a) for the face at xi_j = sign, on an axis-aligned
// hex. faceAxis in {0,1,2} = which parent coordinate is pinned to the face;
// faceSign = +1 or -1.
void bubbleAndGrad(double xi, double eta, double zeta,
                    int faceAxis, double faceSign,
                    double& b, array1d<double>& dbdXi);
 
// 2x2x2 Gauss rule
const double gp = 1.0 / std::sqrt(3.0);
const std::array<std::array<double,3>,8> gaussPts = {{
    {-gp,-gp,-gp},{ gp,-gp,-gp},{ gp, gp,-gp},{-gp, gp,-gp},
    {-gp,-gp, gp},{ gp,-gp, gp},{ gp, gp, gp},{-gp, gp, gp}
}};
const double gaussW = 1.0; // weight 1 each for 2-pt rule per axis
 
// Voigt B-column for a generic scalar shape function's physical gradient
array2d<real64> bColumn(double dNdx, double dNdy, double dNdz);
array2d<real64> computeJacobian(const std::array<array1d<real64>,8>& coords,
                      const std::array<std::array<double,3>,8>& dNdXi);
// ----------------------------------------------------------------------------
// Element geometry: axis-aligned unit cube, origin = (x0,y0,z0).
// Physical gradient = parent gradient * 2/L (L=1 here) -> factor 2.
// detJ = (L/2)^3 = 1/8.
// ----------------------------------------------------------------------------
struct HexElem
{
    // double x0, y0, z0;      // corner with min x,y,z
    std::array<array1d<real64>,8> coords;   // PHYSICAL corner coordinates (arbitrary,
                                  // not assumed axis-aligned or unit-length)

    std::array<int,8> node; // global node ids
    bool hasBubble = false;
    int bubbleId = -1;      // global bubble block index
    int bubbleFaceAxis = 0; // which parent axis the fault face sits on
    double bubbleFaceSign = 1.0;
};
 
// Local nodal-nodal (24x24), nodal-bubble (24x3), bubble-bubble (3x3)
struct LocalBlocks
{
    static constexpr int nnodes = 24;
    static constexpr int ncenters = 3;

    array2d<real64> Kuu{nnodes,nnodes};
    array2d<real64> Kub{nnodes,ncenters};
    array2d<real64> Kbb{ncenters,ncenters};
};
 
LocalBlocks integrateElement(const array2d<real64>& D, const HexElem& e);
struct ContactState { double tN = 0.0; bool open = true; double theta = 0.; double phi = 0.; };
void setElement(const ContactState&  contact, HexElem& elemA, HexElem& elemB, array2d<real64>& RR );
// Kcond = Kuu - Kub * Kbb^-1 * Kub^T  (3x3 inverse via our own LU, small enough
// that we go through the same generic solver rather than a closed form)
array2d<real64> condense(const LocalBlocks& L);

CRSMat assemble(const array2d<real64>& D, const ContactState& contact);

struct SolverStepResult
{
    SolverStepResult(std::vector<double>& u_,double gN_, int newtonIterations_, bool converged_):
    u(u_),gN(gN_),newtonIterations(newtonIterations_),converged(converged_){};
    std::vector<double> u;
    double gN = 0.;
    int newtonIterations = 0;
    bool converged = false; 
};

std::vector<SolverStepResult> solveStep(double imposedUx, ContactState& contact, CRSMat const & Kelastic, 
        double epsN, int kMaxNewton, std::function<ContactState(double,double,double)>& updateNormalTraction);


}




class FrictionDriver : public ConstitutiveDriver
{
public:
  FrictionDriver( const string & name,
                  Group * const parent );

  static string catalogName()
  { return "FrictionDriver"; }

  void postInputInitialization() override;

  bool execute() override;

  void getColumnNames( string_array & columnNames ) const override;

  template< typename FRICTION_TYPE >
  void
  runTest( FRICTION_TYPE & friction,
           const arrayView2d< real64, 1 > & table );

private:
  /**
   * @brief Get the friction model from the catalog
   */
  constitutive::FrictionBase & getFriction();
  constitutive::FrictionBase const & getFriction() const;

  void initializeTable();

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : ConstitutiveDriver::viewKeyStruct
  {
    constexpr static char const * frictionNameString()
    { return "friction"; }

    constexpr static char const * contactNameString()
    { return "contact"; }

    constexpr static char const * displacementFunctionString()
    { return "dispControl"; }

    // constexpr static char const * dJumpFunctionString()
    // { return "dJumpControl"; }

    constexpr static char const * stressFunctionRString()
    { return "stressControlsR"; }
    
    constexpr static char const * stressFunctionLString()
    { return "stressControlsL"; }

    constexpr static char const * thetaString()
    { return "yTiltAngle";}

    constexpr static char const * phiString()
    { return "zTiltAngle";}

    constexpr static char const * normalDispTolFac()
    { return "tolJumpN"; }

    constexpr static char const * normalTractionTolFac()
    { return "tolNormalTrac"; }

    constexpr static char const * slidingTolFac()
    { return "tolJumpT"; }

    constexpr static char const * iterPenNFac()
    { return "iterPenNFac"; }

    constexpr static char const * iterPenTFac()
    { return "iterPenTFac"; }

    //geometry
    constexpr static char const * area()
    { return "faceArea"; }

    constexpr static char const * volume()
    { return "neighborsVolume"; }

    constexpr static char const * shear()
    { return "neighborsShear"; }

    constexpr static char const * bulk()
    { return "neighborsBulk"; }

    constexpr static char const * simultaneous()
    { return "simultaneous"; }

    constexpr static char const * maxNewtonIterString()
    { return "maxNewtonIter"; }

  };

  // Time is defined in base class
//   enum columnKeys { NTRAC=1, STRAC0, STRAC1, NDJUMP, DSLIP0, DSLIP1, CC, FS, 
//                   NEWTRAC, SNEWTRAC0, SNEWTRAC1, 
//                   NJUMP, SLIP0, SLIP1,
//                   NTOL, TTOL, NTRACTOL,
//                   ITERPEN0, ITERPEN1, TLIM,
//                   ITER };
enum columnKeys { NTRAC=1, STRAC0, STRAC1, DISP, NDJUMP, DSLIP0, DSLIP1, CC, FS,
                NEWTRAC, SNEWTRAC0, SNEWTRAC1,
                NJUMP, SLIP0, SLIP1,
                NTOL, TTOL, NTRACTOL,
                ITERPEN0, ITERPEN1,
                FEMNTRAC, FEMGN, FEMNEWTONITER, FEMCONVERGED,   // NEW
                TLIM,
                ITER };

  string m_dispFunctionName; ///<
  // string m_dJumpFunctionName; ///<
  string m_stressFunctionsNamesR; ///<
  string m_stressFunctionsNamesL; ///<

  real64 m_theta{0.0}; ///< x-tilt of fault
  real64 m_phi{0.0};  ///< y-tilt of fault

  real64 m_normalDispTolFac{1.e-8};
  real64 m_normalTracTolFac{100.};
  real64 m_slidingTolFac{1e-5};

  real64 m_iterPenNFac{100};
  real64 m_iterPenTFac{1};

  //geometry for iterPen
  real64 m_area{0.0};
  array1d< real64 > m_volume{};
  array1d< real64 > m_shearModulus{};
  array1d< real64 > m_bulkModulus{};

   integer m_isSimultaneous{1};
   integer m_maxNewtonIter{20};              ///< Maximum Newton iterations

  string m_frictionName;               ///< frictionType identifier
};

}

#endif //GEOS_CONSTITUTIVEDRIVERS_CONTACT_FRICTIONDRIVER_HPP
