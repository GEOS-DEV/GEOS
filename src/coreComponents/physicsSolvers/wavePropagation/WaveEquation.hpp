/*
 * WaveEquation.hpp
 *
 *  Created on: Jan 12, 2021
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVEEQUATION_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVEEQUATION_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class WaveEquation : public SolverBase
{
public:
  WaveEquation( const std::string & name,
                Group * const parent );

  virtual ~WaveEquation() override;

  WaveEquation( WaveEquation const & ) = delete;
  WaveEquation( WaveEquation && ) = default;

  WaveEquation& operator=( WaveEquation const & ) = delete;
  WaveEquation& operator=( WaveEquation && ) = delete;


  static string CatalogName() { return "WaveEquation"; }

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
                           -1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_n,
                           "pressure_n",
                           array1d< real64 >,
                           -1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_np1,
                           "pressure_np1",
                           array1d< real64 >,
                           -1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n+1." );
}


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVEEQUATION_HPP_ */
