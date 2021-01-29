/*
 * WaveEquation.hpp
 *
 *  Created on: Jan 12, 2021
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class AcousticWaveEquationSEM : public SolverBase
{
public:
  AcousticWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~AcousticWaveEquationSEM() override;

  AcousticWaveEquationSEM() = delete;
  AcousticWaveEquationSEM( AcousticWaveEquationSEM const & ) = delete;
  AcousticWaveEquationSEM( AcousticWaveEquationSEM && ) = default;

  AcousticWaveEquationSEM & operator=( AcousticWaveEquationSEM const & ) = delete;
  AcousticWaveEquationSEM & operator=( AcousticWaveEquationSEM && ) = delete;


  static string CatalogName() { return "AcousticSEM"; }

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void RegisterDataOnMesh( Group * const MeshBody ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 SolverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 ExplicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  /// Returns the value of a Ricker at time t0 with central Fourier frequency f0
  virtual
  real64 EvaluateRicker( real64 const & t0, real64 const & f0 );
  /// Returns the value of the second derivative of a Ricker at time t0 with central Fourier frequency f0
  virtual
  real64 EvaluateSecondDerivativeRicker( real64 const & t0, real64 const & f0 );


//  virtual void ApplyBoundaryConditions( real64 const time,
//                                        real64 const dt,
//                                        DomainPartition & domain,
//                                        DofManager const & dofManager,
//                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
//                                        arrayView1d< real64 > const & localRhs ) override;

  /**@}*/


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr auto varName = "varName";
  } waveEquationViewKeys;


  arrayView1d< string const > solidMaterialNames() const { return m_solidMaterialNames; }


protected:
  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;


  array1d< string > m_solidMaterialNames;


};

namespace extrinsicMeshData
{

EXTRINSIC_MESH_DATA_TRAIT( Pressure_nm1,
                           "pressure_nm1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_n,
                           "pressure_n",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_np1,
                           "pressure_np1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS,
                           "rhs",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "RHS" );

EXTRINSIC_MESH_DATA_TRAIT( MassVector,
                           "massVector",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Diagonal Mass Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector,
                           "dampingVector",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Diagonal Damping Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocity,
                           "mediumVelocity",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Medium velocity of the cell" );

EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector,
                           "stiffnessVector",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Stiffness vector contains R_h*Pressure_n." );

}


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_ */
