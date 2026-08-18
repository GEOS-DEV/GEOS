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

#include "FrictionDriver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutiveDrivers/LogLevelsInfo.hpp"
#include "constitutive/contact/FrictionBase.hpp"
#include "constitutive/contact/FrictionSelector.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

#define TEST_PBQ 1

namespace geos
{

  namespace LUsolver{

bool luDecompose(std::vector<double>& A, int n, std::vector<int>& piv)
{
    piv.resize(n);
    for (int i = 0; i < n; ++i) piv[i] = i;
 
    for (int k = 0; k < n; ++k)
    {
        // partial pivot: largest |entry| in column k, rows k..n-1
        int p = k;
        double maxVal = std::abs(A[k*n + k]);
        for (int i = k+1; i < n; ++i)
        {
            double v = std::abs(A[i*n + k]);
            if (v > maxVal) { maxVal = v; p = i; }
        }
        if (maxVal < 1e-13) return false; // singular / near-singular
 
        if (p != k)
        {
            for (int j = 0; j < n; ++j) std::swap(A[k*n+j], A[p*n+j]);
            std::swap(piv[k], piv[p]);
        }
 
        for (int i = k+1; i < n; ++i)
        {
            double factor = A[i*n + k] / A[k*n + k];
            A[i*n + k] = factor; // store multiplier (L entry)
            for (int j = k+1; j < n; ++j)
                A[i*n + j] -= factor * A[k*n + j];
        }
    }
    return true;
}

 
// Solve A x = b given the LU-factored A (from luDecompose) and its pivot.
std::vector<double> luSolve(const std::vector<double>& LU, int n,
                             const std::vector<int>& piv,
                             const std::vector<double>& b)
{
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = b[piv[i]]; // apply permutation
 
    // forward substitution, L has implicit unit diagonal
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            x[i] -= LU[i*n + j] * x[j];
 
    // backward substitution
    for (int i = n-1; i >= 0; --i)
    {
        for (int j = i+1; j < n; ++j)
            x[i] -= LU[i*n + j] * x[j];
        x[i] /= LU[i*n + i];
    }
    return x;
}
 
// Convenience: solve a SMALL fixed-size system (used for the 3x3 bubble
// condensation) via the same LU machinery, working on a flat copy.
template <int N>
std::array<std::array<double,N>,N> invertSmall(const std::array<std::array<double,N>,N>& M)
{
    std::vector<double> A(N*N);
    for (int i=0;i<N;++i) for (int j=0;j<N;++j) A[i*N+j] = M[i][j];
    std::vector<int> piv;
    bool ok = luDecompose(A, N, piv);
    if (!ok) std::cerr << "WARNING: singular bubble block!\n";
 
    std::array<std::array<double,N>,N> Inv{};
    for (int col=0; col<N; ++col)
    {
        std::vector<double> e(N, 0.0); e[col] = 1.0;
        std::vector<double> x = luSolve(A, N, piv, e);
        for (int row=0; row<N; ++row) Inv[row][col] = x[row];
    }
    return Inv;
}

};//end LUSolver

namespace mimicAssembly{

  array2d<real64> buildD(double K, double G)
  {
      array2d<real64> D{6,6};

      D[0][0] = D[1][1] = D[2][2] = K + 4/3*G;
      D[0][1] = D[1][0] = D[0][2] = D[2][0] = D[1][2] = D[2][1] = K - 2./3 * G;
      D[3][3] = D[4][4] = D[5][5] = G;


      return D;
  }

void shapeAndGrad(double xi, double eta, double zeta,
                   std::array<double,8>& N,
                   std::array<std::array<double,3>,8>& dNdXi)
{
    for (int a = 0; a < 8; ++a)
    {
        double xa = refCorners[a][0], ya = refCorners[a][1], za = refCorners[a][2];
        N[a] = 0.125 * (1+xa*xi) * (1+ya*eta) * (1+za*zeta);
        dNdXi[a][0] = 0.125 * xa * (1+ya*eta) * (1+za*zeta);
        dNdXi[a][1] = 0.125 * (1+xa*xi) * ya * (1+za*zeta);
        dNdXi[a][2] = 0.125 * (1+xa*xi) * (1+ya*eta) * za;
    }
}
 
// Face bubble (Frigo Eq. 20a) for the face at xi_j = sign, on an axis-aligned
// hex. faceAxis in {0,1,2} = which parent coordinate is pinned to the face;
// faceSign = +1 or -1.
void bubbleAndGrad(double xi, double eta, double zeta,
                    int faceAxis, double faceSign,
                    double& b, array1d<double>& dbdXi)
{
    std::array<double,3> c{xi, eta, zeta};
    double linear = 0.5 * (1.0 + faceSign * c[faceAxis]);
    double dlinear = 0.5 * faceSign; // wrt c[faceAxis]
 
    // product of (1 - c_i^2) over the two axes orthogonal to faceAxis
    std::array<int,2> other;
    { int k=0; for (int i=0;i<3;++i) if (i!=faceAxis) other[k++]=i; }
 
    double q0 = 1.0 - c[other[0]]*c[other[0]];
    double q1 = 1.0 - c[other[1]]*c[other[1]];
 
    b = linear * q0 * q1;
 
    array1d<real64> grad{};
    grad[faceAxis]  = dlinear * q0 * q1;
    grad[other[0]]  = linear * (-2.0 * c[other[0]]) * q1;
    grad[other[1]]  = linear * q0 * (-2.0 * c[other[1]]);
    dbdXi = grad;//[1x3]
}

// Voigt B-column for a generic scalar shape function's physical gradient
array2d<real64> bColumn(double dNdx, double dNdy, double dNdz)
{
    array2d<real64> Bc{6,3}; Bc.zero();
    // Bc.resize(6,3);

    Bc(0,0)=dNdx;             Bc(1,1)=dNdy;             Bc(2,2)=dNdz;
    Bc(3,1)=dNdz; Bc(3,2)=dNdy;
    Bc(4,0)=dNdz; Bc(4,2)=dNdx;
    Bc(5,0)=dNdy; Bc(5,1)=dNdx;
    return Bc;
}

LocalBlocks integrateElement(const array2d<real64>& D, const HexElem& e)
{
    LocalBlocks L;
    L.Kuu.zero(); L.Kub.zero(); L.Kbb.zero();
 
    const double scale = 2.0;      // dPhysical/dParent for unit-length edges
    const double detJ  = 1.0/8.0;
    const double w = 1.0;
 
    for (const auto& gpt : gaussPts)
    {
        double xi=gpt[0], eta=gpt[1], zeta=gpt[2];
        std::array<double,8> N;
        std::array<std::array<double,3>,8> dNdXi;
        shapeAndGrad(xi, eta, zeta, N, dNdXi);//TODO use GEOS's function as in H1_TriangleFace_Lagrange1_Gauss.hpp
        array2d<real64> mat(3,3); mat.zero();
        array2d<real64> mat36(6,3); mat36.zero();
 
        array2d<real64> Ba[8];
        // L.Kuu.resize(L.nnodes,L.nnodes);
        // L.Kbb.resize(L.ncenters,L.ncenters);
        // L.Kub.resize(L.nnodes,L.ncenters);

        for (int a=0;a<8;++a){
            Ba[a].resize(6,3);
            Ba[a] = bColumn(dNdXi[a][0]*scale, dNdXi[a][1]*scale, dNdXi[a][2]*scale);}
 
        for (int a=0;a<8;++a)
            for (int b=0;b<8;++b){
                // L.Kuu.block<3,3>(3*a,3*b) += gaussW*detJ * (Ba[a].transpose()*D*Ba[b]);
                mat.zero();
                #if TEST_PBQ
                LvArray::tensorOps::Rij_eq_PkiBklQlj<6,3>(mat,Ba[a],D,Ba[b]);
                #else
                mat36.zero();
                LvArray::tensorOps::Rij_eq_AkiBkj<3,6,6>(mat36, Ba[a],D);
                LvArray::tensorOps::Rij_eq_AikBkj<3,3,6>(mat,mat36,Ba[b]);
                #endif
                LvArray::tensorOps::scale<3,3>(mat, w*gaussW*detJ);
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        L.Kuu[3*a + i][3*b + j ]+=mat[i][j];
            }
 
        if (e.hasBubble)
        {
            double bub; array1d<real64> dbdXi;
            bubbleAndGrad(xi,eta,zeta, e.bubbleFaceAxis, e.bubbleFaceSign, bub, dbdXi);
            auto Bb = bColumn(dbdXi[0]*scale, dbdXi[1]*scale, dbdXi[2]*scale);
 
            for (int a=0;a<8;++a){
                // L.Kub.block<3,3>(3*a,0) += gaussW*detJ * (Ba[a].transpose()*D*Bb);
                mat.zero();
                #if TEST_PBQ
                LvArray::tensorOps::Rij_eq_PkiBklQlj<6,3>(mat,Ba[a],D,Bb);
                #else
                mat36.zero();
                LvArray::tensorOps::Rij_eq_AkiBkj<3,6,6>(mat36, Ba[a],D);
                LvArray::tensorOps::Rij_eq_AikBkj<3,3,6>(mat,mat36,Bb);
                #endif
                LvArray::tensorOps::scale<3,3>(mat, w*gaussW*detJ);
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        L.Kuu[3*a + i][ j ]+=mat[i][j];
            }


 
            // L.Kbb += gaussW*detJ * (Bb.transpose()*D*Bb);
            mat.zero();
            #if TEST_PBQ
            LvArray::tensorOps::Rij_eq_PkiBklPlj<6,3>(mat,Bb,D);
            #else
            mat36.zero();
            LvArray::tensorOps::Rij_eq_AkiBkj<3,6,6>(mat36, Bb,D);
            LvArray::tensorOps::Rij_eq_AikBkj<3,3,6>(mat,mat36,Bb);
            #endif
            LvArray::tensorOps::scale<3,3>(mat, w*gaussW*detJ);
            LvArray::tensorOps::add<3,3>(L.Kbb, mat);
        }
    }
    return L;
}

// Kcond = Kuu - Kub * Kbb^-1 * Kub^T  (3x3 inverse via our own LU, small enough
// that we go through the same generic solver rather than a closed form)
array2d<real64> condense(const LocalBlocks& L)
{
    // std::vector<double> Kbb3(9);
    // for (int i=0;i<3;++i) for (int j=0;j<3;++j) Kbb3[i*3+j] = L.Kbb(i,j);
    // std::vector<int> piv;
    // luDecompose(Kbb3, 3, piv);
 
    array2d<real64> KbbInv(3,3);
    LvArray::tensorOps::invert<L.ncenters>(KbbInv, L.Kbb);
    // for (int col=0; col<3; ++col)
    // {
    //     std::vector<double> e(3,0.0); e[col]=1.0;
    //     auto x = luSolve(Kbb3, 3, piv, e);
    //     for (int row=0; row<3; ++row) KbbInv(row,col) = x[row];
    // }
 
    array2d<real64> Kcond(24,24), mat(24,24), mat243(24,3);
    LvArray::tensorOps::copy<L.nnodes,L.nnodes>(Kcond,L.Kuu);
    // Eigen::Matrix<double,24,24> Kcond = L.Kuu - L.Kub * KbbInv * L.Kub.transpose();
    //TODO (test - test - test)
    #if TEST_PBQ
    LvArray::tensorOps::Rij_eq_PkiBklPlj<L.nnodes,L.ncenters>(mat,L.Kub,KbbInv);
    LvArray::tensorOps::scaledAdd<L.nnodes, L.nnodes>(Kcond,mat,-1);
    #else// main road
    LvArray::tensorOps::Rij_eq_AikBkj<L.nnodes,L.ncenters,L.ncenters>(mat243, L.Kub, KbbInv);
    LvArray::tensorOps::Rij_eq_AikBkj<L.nnodes,L.nnodes,L.ncenters>(mat,mat243,L.Kub);
    LvArray::tensorOps::scaledAdd<L.nnodes,L.nnodes>(Kcond,mat,-1);
    #endif

    // for (int r=0;r<24;++r)
    //     for (int c=0;c<24;++c)
    //     {
    //         double s=0.0;
    //         for (int k1=0;k1<3;++k1)
    //             for (int k2=0;k2<3;++k2)
    //                 s += L.Kub(r,k1)*KbbInv(k1,k2)*L.Kub(c,k2);
    //         Kcond(r,c) = L.Kuu(r,c) - s;
    //     }
    return Kcond;
}

CRSMat assemble(const array2d<real64>& D)
{
    HexElem elemA, elemB;
    elemA.x0=0; elemA.y0=0; elemA.z0=0;
    elemA.node = {0,1,2,3,4,5,6,7};
    elemA.hasBubble = true; elemA.bubbleId = 0;
    elemA.bubbleFaceAxis = 0; elemA.bubbleFaceSign = +1.0; // xi=+1 face -> x=1
 
    elemB.x0=1; elemB.y0=0; elemB.z0=0;
    elemB.node = {8,9,10,11,12,13,14,15};
    elemB.hasBubble = true; elemB.bubbleId = 1;
    elemB.bubbleFaceAxis = 0; elemB.bubbleFaceSign = -1.0; // xi=-1 face -> x=1
 
    const int nNodes = 16, nDof = nNodes*3;
 
    // ---------------- assemble & condense each element ---------------
    // Eigen::Matrix<double,24,24> KcondA = assembleCondensed(elemA);
    // Eigen::Matrix<double,24,24> KcondB = assembleCondensed(elemB);
    array2d<real64> KcondA = mimicAssembly::condense(integrateElement(D,elemA));
    array2d<real64> KcondB = mimicAssembly::condense(integrateElement(D,elemB));
 
    // ---------------- GEOS-style two-phase CRSMatrix assembly ----------------
    // Phase 1 (symbolic): reserve every (row,col) pair that will ever be
    // written -- elastic pairs from both elements, PLUS the 4x4 cross-face
    // contact coupling pairs (needed even though epsN only activates them
    // conditionally, since CRSMatrix pre-allocates its sparsity).
    CRSMat Kelastic;
    Kelastic.resize(nDof, nDof, 24); // generous initial row capacity
 
    auto elemDofs = [&](const HexElem& e)
    {
        std::array<IndexType,24> dofs;
        for (int a=0;a<8;++a) for (int c=0;c<3;++c) dofs[3*a+c] = 3*e.node[a]+c;
        return dofs;
    };
    auto dofsA = elemDofs(elemA);
    auto dofsB = elemDofs(elemB);
 
    for (auto row : dofsA) for (auto col : dofsA) Kelastic.insertNonZero(row, col, 0.0);
    for (auto row : dofsB) for (auto col : dofsB) Kelastic.insertNonZero(row, col, 0.0);
    

    // NOTE: the contact tangent (state-dependent, rebuilt every Newton
    // iteration) is deliberately NOT part of Kelastic's sparsity -- it's a
    // tiny 8x8 dense coupling added directly to the dense Jacobian buffer
    // each iteration (see solveStep). Kelastic only ever holds the constant
    // elastic operator.
    
    // Phase 2 helper (numeric): accumulate the CONSTANT elastic contribution
    // via addToRow -- exactly the GEOS kernel pattern (columns must be sorted
    // per row, which our simple 8-node connectivity already is after a sort).
    auto scatterElastic = [&](const std::array<IndexType,24>& dofs, const array2d<real64>& Kc)
    {
        for (int a=0;a<24;++a)
        {
            std::vector<IndexType> cols(dofs.begin(), dofs.end());
            std::vector<double> vals(24);
            for (int b=0;b<24;++b) vals[b] = Kc(a,b);
            // sort cols and vals together (addToRow requires sorted columns)
            std::vector<int> order(24); for (int i=0;i<24;++i) order[i]=i;
            std::sort(order.begin(), order.end(), [&](int i,int j){return cols[i]<cols[j];});
            std::vector<IndexType> colsSorted(24); std::vector<double> valsSorted(24);
            for (int i=0;i<24;++i){ colsSorted[i]=cols[order[i]]; valsSorted[i]=vals[order[i]]; }
            Kelastic.addToRow< RAJA::seq_atomic >(dofs[a], colsSorted.data(), valsSorted.data(), 24);
        }
    };
    scatterElastic(dofsA, KcondA);
    scatterElastic(dofsB, KcondB);

    // Kelastic now holds the constant, condensed elastic operator. Each
    // Newton iteration below extracts it to a dense buffer and adds the
    // state-dependent contact tangent on top -- Kelastic itself is never
    // touched again.

    return Kelastic;
};

std::vector<double> solveStep(double imposedUx, ContactState& contact, CRSMat const & Kelastic, double epsN, std::function<ContactState(double,double,double)>& updateNormalTraction )
    {
        //prep -- duplicate (pass element A and B ?)
        HexElem elemA, elemB;
        elemA.x0=0; elemA.y0=0; elemA.z0=0;
        elemA.node = {0,1,2,3,4,5,6,7};
        elemA.hasBubble = true; elemA.bubbleId = 0;
        elemA.bubbleFaceAxis = 0; elemA.bubbleFaceSign = +1.0; // xi=+1 face -> x=1
    
        elemB.x0=1; elemB.y0=0; elemB.z0=0;
        elemB.node = {8,9,10,11,12,13,14,15};
        elemB.hasBubble = true; elemB.bubbleId = 1;
        elemB.bubbleFaceAxis = 0; elemB.bubbleFaceSign = -1.0; // xi=-1 face -> x=1

        // const int nNodes = 16, nDof = nNodes*3;
        const auto nDof = Kelastic.numRows();
        // elemA reference corners with xi=+1 are local nodes {1,2,5,6} -> fault face "-"
        // elemB reference corners with xi=-1 are local nodes {0,3,4,7} -> fault face "+"
        std::array<int,4> faceLocalMinus = {1,2,5,6};
        std::array<int,4> faceLocalPlus  = {0,3,4,7};
        std::array<int,4> faceGlobalMinus, faceGlobalPlus;
        for (int i=0;i<4;++i)
        {
            faceGlobalMinus[i] = elemA.node[faceLocalMinus[i]];
            faceGlobalPlus[i]  = elemB.node[faceLocalPlus[i]];
        }
        std::array<std::array<IndexType,3>,4> faceDofMinus, faceDofPlus;
        for (int i=0;i<4;++i)
        for (int c=0;i<3;++c)
        {
            faceDofMinus[i][c] = 3*faceGlobalMinus[i]+c;
            faceDofPlus[i][c]  = 3*faceGlobalPlus[i]+c;
        }

        std::vector<double> u(nDof, 0.0);
        double tNold = contact.tN;
        const double faceWeight = 1.0/4.0;

        auto denseFromCRS = [&](const CRSMat& M)
        {
            std::vector<double> A(nDof*nDof, 0.0);
            auto view = M.toViewConst();
            for (IndexType row=0; row<nDof; ++row)
            {
                auto cols = view.getColumns(row);
                auto vals = view.getEntries(row);
                for (int k=0;k<cols.size();++k) A[row*nDof + cols[k]] += vals[k];
            }
            return A;
        };

        for (int newtonIt = 0; newtonIt < 20; ++newtonIt)
        {
            double jumpAvg[3] = {0.,0.,0.};
            for (int i=0;i<4;++i) for(int c=0; c<3; ++c) jumpAvg[c] += (u[faceDofPlus[i][c]] - u[faceDofMinus[i][c]])/4;
            double gN = LvArray::tensorOps::AiBi(jumpAvg,fault.nf);
 
            ContactState trial = updateNormalTraction(tNold, gN, epsN);
 
            // residual r = K_elastic*u + contact virtual work (tangent is
            // NOT part of this -- residual uses ONLY the elastic operator,
            // matching Frigo's r_u = (elastic term) + (traction term); the
            // tangent stiffness kc belongs to the JACOBIAN only, never to
            // the residual's K*u term)
            std::vector<double> Kdense = denseFromCRS(Kelastic);
            std::vector<double> r(nDof, 0.0);
            for (int i=0;i<nDof;++i)
            {
                double s=0.0;
                for (int j=0;j<nDof;++j) s += Kdense[i*nDof+j]*u[j];
                r[i]=s;
            }
            for (int i=0;i<4;++i)
            {
                r[faceDofPlus[i]]  += faceWeight*trial.tN*fault.nf[c];
                r[faceDofMinus[i]] -= faceWeight*trial.tN*fault.nf[c];
            }
 
            // Jacobian = Kglobal + contact tangent (if closed)
            std::vector<double> J = Kdense;
            if (!trial.open)
            {
                double kc = faceWeight*epsN/4.0;
                for (int i=0;i<4;++i)
                    for (int j=0;j<4;++j)
                for (int a=0;a<3;++a)
                    for (int b=0;b<3;++b)
                    {
                        int pi=faceDofPlus[i][a],  pj=faceDofPlus[j][b];
                        int mi=faceDofMinus[i][a], mj=faceDofMinus[j][b];
                        J[pi*nDof+pj] += kc * fault.nf[a] * fault.nf[b];
                        J[pi*nDof+mj] -= kc* fault.nf[a] * fault.nf[b];
                        J[mi*nDof+pj] -= kc* fault.nf[a] * fault.nf[b];
                        J[mi*nDof+mj] += kc* fault.nf[a] * fault.nf[b];
                    }
            }

            // BCs
            std::array<int,4> fixedLocal  = {0,3,4,7};
            std::array<int,4> drivenLocal = {1,2,5,6};
            std::vector<IndexType> fixedDofs, drivenDofsX, drivenDofsYZ;
            for (int ln : fixedLocal) for (int c=0;c<3;++c) fixedDofs.push_back(3*elemA.node[ln]+c);
            for (int ln : drivenLocal)
            {
                drivenDofsX.push_back(3*elemB.node[ln]+0);
                drivenDofsYZ.push_back(3*elemB.node[ln]+1);
                drivenDofsYZ.push_back(3*elemB.node[ln]+2);
            }
            std::vector<bool> isFixed(nDof,false);
            for (auto d: fixedDofs) isFixed[d]=true;
            for (auto d: drivenDofsYZ) isFixed[d]=true;
            for (auto d: drivenDofsX) isFixed[d]=true;
            
            // apply BCs by row/col elimination
            std::vector<double> rhs(nDof);
            for (int i=0;i<nDof;++i) rhs[i] = -r[i];
            for (int d : drivenDofsX) rhs[d] = imposedUx - u[d];
            for (int d : fixedDofs)    rhs[d] = 0.0 - u[d];
            for (int d : drivenDofsYZ) rhs[d] = 0.0 - u[d];
 
            for (int row=0; row<nDof; ++row)
            {
                if (!isFixed[row]) continue;
                for (int col=0; col<nDof; ++col) J[row*nDof+col] = (col==row) ? 1.0 : 0.0;
            }
            for (int col=0; col<nDof; ++col)
            {
                if (!isFixed[col]) continue;
                for (int row=0; row<nDof; ++row)
                    if (!isFixed[row]) J[row*nDof+col] = 0.0;
            }
 
            std::vector<int> piv;
            bool ok = LUsolver::luDecompose(J, (int)nDof, piv);
            if (!ok) { std::cerr << "singular Jacobian!\n"; break; }
            std::vector<double> du = LUsolver::luSolve(J, (int)nDof, piv, rhs);
 
            double duNorm = 0.0;
            for (int i=0;i<nDof;++i) { u[i]+=du[i]; duNorm += du[i]*du[i]; }
            contact.tN = trial.tN; contact.open = trial.open;
 
            if (std::sqrt(duNorm) < 1e-12) break;
        }
        return u;
    };

};//end mimicAssembly



using namespace dataRepository;
using namespace constitutive;

FrictionDriver::FrictionDriver( const string & name, Group * const parent )
  : ConstitutiveDriver( name, parent )
{
  registerWrapper( viewKeyStruct::frictionNameString(), &m_frictionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Friction model to test" );

  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of sample step to take in both jumps and traction increments" );

  // registerWrapper( viewKeyStruct::jumpFunctionString(), &m_jumpFunctionName ).
  //   setInputFlag( InputFlags::REQUIRED ).
  //   setDescription( "Name of the input function representing jump function along world x-axis" );

  registerWrapper( viewKeyStruct::dJumpFunctionString(), &m_dJumpFunctionName ).//should be derive from convergence history
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input function representing deltaDisplacementJump function along world x-axis" );

  registerWrapper( viewKeyStruct::stressFunctionLString(), &m_stressFunctionsNamesL ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input functions representing stresses in the left cell");
  
  registerWrapper( viewKeyStruct::stressFunctionRString(), &m_stressFunctionsNamesR ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input functions representing stresses in the right cell");

  registerWrapper( viewKeyStruct::thetaString(), &m_theta ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "y-Tilt angle in degree" );

  registerWrapper( viewKeyStruct::phiString(), &m_phi ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue(0.).
    setDescription( "z-Tilt angle in degree" );

  //first batch of parameters
  registerWrapper( viewKeyStruct::normalDispTolFac(), &m_normalDispTolFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "normal Displacement Tolerance (scale as inverse of average Young modulus)." );

  registerWrapper( viewKeyStruct::normalTractionTolFac(), &m_normalTracTolFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "normal Traction Tolerance" );

  registerWrapper( viewKeyStruct::slidingTolFac(), &m_slidingTolFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "tangential Displacement Tolerance" );

  registerWrapper( viewKeyStruct::iterPenNFac(), &m_iterPenNFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "normal Penalty Factor" );

  registerWrapper( viewKeyStruct::iterPenTFac(), &m_iterPenTFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "tangential Penatly Factor" );


  //geometry
  registerWrapper( viewKeyStruct::area(), &m_area ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Face area" );

  registerWrapper( viewKeyStruct::volume(), &m_volume ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Neighboring cells' volume" );

  registerWrapper( viewKeyStruct::bulk(), &m_bulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Neighboring cells' bulk Modulus" );

  registerWrapper( viewKeyStruct::shear(), &m_shearModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Neighboring cells' shear Modulus" );

  //algo tune
  registerWrapper( viewKeyStruct::simultaneous(), &m_isSimultaneous ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "isSimultaneous" );

  addLogLevel< logInfo::LogOutput >();
}

void FrictionDriver::postInputInitialization()
{
  ConstitutiveDriver::postInputInitialization();

  // Check that the functions exist
  FunctionManager & functionManager = FunctionManager::getInstance();
  // GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_jumpFunctionName ),
  //                GEOS_FMT( "Jump function with name '{}' not found", m_jumpFunctionName ),
  //                getWrapperDataContext( viewKeyStruct::jumpFunctionString() ) );

  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_dJumpFunctionName ),
                 GEOS_FMT( "dJump function with name '{}' not found", m_dJumpFunctionName ),
                 getWrapperDataContext( viewKeyStruct::dJumpFunctionString() ) );

  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_stressFunctionsNamesR ),
                 GEOS_FMT( "Stress functions with name '{}' not found", m_stressFunctionsNamesR ),
                 getWrapperDataContext( viewKeyStruct::stressFunctionRString() ) );
  
  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_stressFunctionsNamesL ),
                 GEOS_FMT( "Stress functions with name '{}' not found", m_stressFunctionsNamesL ),
                 getWrapperDataContext( viewKeyStruct::stressFunctionLString() ) );

  string_array columnNames;
  getColumnNames( columnNames );
  integer const numCols = static_cast< integer >(columnNames.size());

  // initialize functions
  // TableFunction & jumpFunction = functionManager.getGroup< TableFunction >( m_jumpFunctionName );
  TableFunction & dJumpFunction = functionManager.getGroup< TableFunction >( m_dJumpFunctionName );
  // TableFunction & tractionFunction = functionManager.getGroup< TableFunction >( m_tractionFunctionName );

  // jumpFunction.initializeFunction();
  dJumpFunction.initializeFunction();
  // tractionFunction.initializeFunction();

  // // TODO: Maybe we should take the maximum extent of jumpFunction and tractionFunction
  ArrayOfArraysView< real64 > coordinates = dJumpFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];

  // Allocate the data
  allocateTable( numCols, minTime, maxTime );

  // set input columns
  initializeTable();
}

bool FrictionDriver::execute()
{
  FrictionBase & baseFriction = getFriction();

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching Friction Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Friction ............... " << m_frictionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Type ................... " << baseFriction.getCatalogName() );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Steps .................. " << m_numSteps );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ................. " << m_outputFile );

  // create a dummy discretization with one quadrature point for
  // storing constitutive data
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  integer const numRows = m_table.size( 0 );
  discretization.resize( numRows ); // numRows elements
  baseFriction.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  constitutiveUpdatePassThru( baseFriction, [&]( auto & selectedFrictionModel )
  {
    using FRICTION_TYPE = TYPEOFREF( selectedFrictionModel );
    runTest< FRICTION_TYPE >( selectedFrictionModel, m_table );
  } );

  return false;
}

void FrictionDriver::getColumnNames( string_array & columnNames ) const
{
  columnNames.emplace_back( "time" );
  columnNames.emplace_back( "traction,normal" );
  columnNames.emplace_back( "traction,tangent1" );
  columnNames.emplace_back( "traction,tangent2" );
  columnNames.emplace_back( "delta displacement jump,normal" );
  columnNames.emplace_back( "delta displacement jump,tangent1" );
  columnNames.emplace_back( "delta displacement jump,tangent2" );
  columnNames.emplace_back( "encoded constaint (0:converged, 1:stick & gn>0 (opening), 2: interpenetration, 3: stick & gt>lim (disp-sliding), 4: tau>taulim (trac-sliding) )" );
  columnNames.emplace_back( "fracture state (0:stick, 1:slip , 2: new slip, 3: open)" );
  columnNames.emplace_back( "newtraction,normal" );
  columnNames.emplace_back( "newtraction,tangent1" );
  columnNames.emplace_back( "newtraction,tangent2" );
  columnNames.emplace_back( "displacement jump,normal" );
  columnNames.emplace_back( "displacement jump,tangent1" );
  columnNames.emplace_back( "displacement jump,tangent2" );
  columnNames.emplace_back( "derived disp tol, normal" );
  columnNames.emplace_back( "derived disp tol, tangent" );
  columnNames.emplace_back( "derived traction tol, normal" );
  columnNames.emplace_back( "iterative penalty, normal" );
  columnNames.emplace_back( "iterative penalty, tangent" );
  columnNames.emplace_back( "iteration" );


  if( dynamic_cast< CoulombFriction const * >(&getFriction()) != nullptr )
  {
    columnNames.emplace_back( "tau limit" );
  }
}

void FrictionDriver::initializeTable()
{
  integer const numRows = m_table.size( 0 );

  FunctionManager & functionManager = FunctionManager::getInstance();
  // TableFunction const & jumpFunction = functionManager.getGroup< TableFunction >( m_jumpFunctionName );
  TableFunction const & dJumpFunction = functionManager.getGroup< TableFunction >( m_dJumpFunctionName );

  TableFunction const & tractionFunctionR00 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesR );//should be Voigt full from stress
  TableFunction const & tractionFunctionR11 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesR );//1
  TableFunction const & tractionFunctionR22 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesR );//2
  TableFunction const & tractionFunctionR21 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesR );//3
  TableFunction const & tractionFunctionR02 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesR );//4
  TableFunction const & tractionFunctionR01 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesR );//5
  //
  TableFunction const & tractionFunctionL00 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesL );
  TableFunction const & tractionFunctionL11 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesL );
  TableFunction const & tractionFunctionL22 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesL );
  TableFunction const & tractionFunctionL21 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesL );
  TableFunction const & tractionFunctionL02 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesL );
  TableFunction const & tractionFunctionL01 = functionManager.getGroup< TableFunction >( m_stressFunctionsNamesL );

  real64 const cos_theta = cos( m_theta * M_PI/180.0 );
  real64 const sin_theta = sin( m_theta * M_PI/180.0 );

  real64 const cos_phi = cos( m_phi * M_PI/180.0 );
  real64 const sin_phi = sin( m_phi * M_PI/180.0 );

  for( integer index = 0; index < numRows; ++index )
  {
    real64 const time = m_table( index, TIME );

    real64 const leftStress[6]{ tractionFunctionR00.evaluate(&time),
                                tractionFunctionR11.evaluate(&time), 
                                tractionFunctionR22.evaluate(&time), 
                                tractionFunctionR21.evaluate(&time), 
                                tractionFunctionR02.evaluate(&time), 
                                tractionFunctionR01.evaluate(&time), 
                              };//should be full 6-Voigt

    real64 const rightStress[6]{ tractionFunctionL00.evaluate(&time),
                                tractionFunctionL11.evaluate(&time), 
                                tractionFunctionL22.evaluate(&time), 
                                tractionFunctionL21.evaluate(&time), 
                                tractionFunctionL02.evaluate(&time), 
                                tractionFunctionL01.evaluate(&time), 
                              };//should be full 6-Voigt
    
    real64 const n[3]{cos_theta,sin_theta*sin_phi,-cos_theta*sin_phi};

    //rotate and project stress
    for(const auto& sigma : {leftStress, rightStress }){
      
      constexpr real64 NCELLS = 2.0;

      const real64 sigman[3]{
        sigma[0]*n[0]+sigma[5]*n[1]+sigma[4]*n[2],
        sigma[5]*n[0]+sigma[1]*n[1]+sigma[3]*n[2],
        sigma[4]*n[0]+sigma[3]*n[1]+sigma[2]*n[2],
          };

      m_table( index, NTRAC )  += 1./NCELLS*(sigman[0]*cos_phi + sigman[1]*sin_phi*sin_theta + sigman[2]*cos_theta*sin_phi);
      m_table( index, STRAC0 ) += 1./NCELLS*(sigman[1]*cos_theta - sigman[2]*sin_theta);
      m_table( index, STRAC1 ) += 1./NCELLS*(-sigman[0]*sin_phi + sigman[1]*sin_theta*cos_phi + sigman[2]*cos_theta*cos_phi);

    }

    //Project and average Left and Right stress
    //jump = tLocal/(kn kt)
    // real64 const jump = jumpFunction.evaluate( &time );
    // m_table( index, NJUMP ) = m_table( index, NTRAC )/m_iterPenNFac;//might need compute tol --> move to evaluation part
    // m_table( index, SLIP0 ) = m_table( index, STRAC0 )/m_iterPenTFac;
    // m_table( index, SLIP1 ) = m_table( index, STRAC1 )/m_iterPenTFac;

    real64 const dJump = dJumpFunction.evaluate( &time );
    m_table( index, NDJUMP ) = dJump*cos_phi*cos_theta;
    m_table( index, DSLIP0 ) = dJump*cos_theta*sin_phi;
    m_table( index, DSLIP1 ) = dJump*sin_theta;

    m_table( index, CC )    = 0;
    m_table( index, FS )    = fields::contact::FractureState::Stick;
  }

  if( CoulombFriction const * coulombFriction = dynamic_cast< CoulombFriction const * >(&getFriction()) )
  {
    real64 const cohesion = coulombFriction->getCohesion();
    real64 const frictionCoeff = coulombFriction->getFrictionCoeff();
    for( integer index = 0; index < numRows; ++index )
    {
      real64 const normal_traction = m_table( index, NTRAC );
      m_table( index, TLIM ) = cohesion - normal_traction * frictionCoeff;
    }
  }
}

FrictionBase & FrictionDriver::getFriction()
{
  return getConstitutiveManager().getGroup< FrictionBase >( m_frictionName );
}

FrictionBase const & FrictionDriver::getFriction() const
{
  return getConstitutiveManager().getGroup< FrictionBase >( m_frictionName );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        FrictionDriver,
                        string const &, dataRepository::Group * const )

}
