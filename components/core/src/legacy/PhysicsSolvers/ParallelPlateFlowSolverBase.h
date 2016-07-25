//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file ParallelPlateFlowSolverBase.h
 * @author settgast1
 * @date Aug 3, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVERBASE_H_
#define PARALLELPLATEFLOWSOLVERBASE_H_

#include "SolverBase.h"
#include "DataStructures/Tables/Table.h"


#ifdef SRC_EXTERNAL
#include "PhysicsSolvers/FractureFlow/LeakoffModels.h"
#endif

namespace PPFS
{
class FluidEOSBase;
}

class ParallelPlateFluidModelBase{

    public:
	ParallelPlateFluidModelBase(){/** empty **/};
    virtual ~ParallelPlateFluidModelBase(){/** empty **/};

	virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn){};
	virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                              const realT apb,const realT w, const realT qMag, const realT SHP_FCT) = 0;
	virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                  const realT qMag, const realT SHP_FCT) = 0; // one sided

};





class PowerlawFluidModel: public ParallelPlateFluidModelBase{

  public:
	PowerlawFluidModel():
		ParallelPlateFluidModelBase(), m_n(),m_M(),m_phiM(){/** empty **/};
	~PowerlawFluidModel(){/** empty **/};

	virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
	virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                                const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
	virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                    const realT qMag, const realT SHP_FCT); // one sided

	realT m_n; // fluid behavior index
	realT m_M; // consistency index
	realT m_phiM;// modified consistency index (phi accounts for flow between parallel plates)
};

/*
 * The Herschel Bulkley model describes a non-Newtonian fluid with a Bingham law-like yield stress
 * accompanied by a power-law fluid like growth in stress as a function of shear strain rate:
 *
 * \tau = \tau_y + k |du/dr|^n   for |\tau| > \tau_y
 * du/dr = 0                     for |\tau| < \tau_y
 *
 * The parallel plate model implemented here follows the approximate solution provided in
 * Wang and Gordaninejad 1999, Flow analysis of field controllable, electro- and magneto-rheological fluids
 * using Herschel-Bulkley model.
 *
 */
class HerschelBulkleyParallelPlateFluidModel: public PowerlawFluidModel{

    public:
	  HerschelBulkleyParallelPlateFluidModel():
	  	  PowerlawFluidModel(),m_tau_y(0){/** empty **/};
	  ~HerschelBulkleyParallelPlateFluidModel(){/** empty **/};

	virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
	virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                              const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
	virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                  const realT qMag, const realT SHP_FCT); // one sided

	realT m_tau_y; // fluid yield stress
	realT m_phi; // modified consistency index factor

    private:
	Table<1, realT> m_zp_kappaLookup;

	realT CalculateZp(realT kappa, realT n);
	/*
	 lookup table relating the dimensionless plug thickness (zp) to the dimensionless variable
	 kappa = K q^n where
	 K = ((2n+1)/n)^n 2^n/h^(n+1) * M/\tau_y = phi/2* 1/h^(n+1) *M/\tau_y
	 */
	// functions used in Newton's method
	realT analyticalVFunc(realT zp, realT n);
	realT dVdz(realT zp, realT n);

};


class ParallelPlateFlowSolverBase: public SolverBase
{
public:

  ParallelPlateFlowSolverBase( const std::string& name,
                               ProblemManagerT* const pm );

  virtual ~ParallelPlateFlowSolverBase();


  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn) ;
  virtual void RegisterFields( PhysicalDomainT& domain ) = 0 ;
  virtual void InitializeCommunications( PartitionBase& partition );
  virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);
  // realT CalculateMultiphaseViscosity(PhysicalDomainT& domain, lArray1d& faces);
  // void CalculateMultiphaseDensityAndBulkModulus(PhysicalDomainT& domain, localIndex kf);

  void AdvectMultiphaseFields(PhysicalDomainT& domain, SpatialPartition& partition, realT time, realT dt);
  void UpdateMultiphasePointers(PhysicalDomainT& domain);   
  void GenerateMultiphaseArray(rArrayPtrs& multiphasePtrs, localIndex kf, rArray1d& multiphase);
  void GenerateMultiphaseBCArray(rArray1d& faceVolumeFraction, realT time);

  // void BoundaryTracerValues( PhysicalDomainT& domain,
  //                            ObjectDataStructureBaseT& object,
  //                            BoundaryConditionBase* const bc, 
  //                            const lSet& aset, 
  //                            const realT time, 
  //                            const realT dt );

  void BoundaryTracerValuesToSet( PhysicalDomainT& domain, const realT time );
  
  realT BoundedAperture( const realT& aper );
  realT BoundedApertureDerivative( const realT& aper );



  virtual void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain)
  {
    (void)time;
    (void)dt;
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::UpdateEOS() not overridden");
  }


  // functions for implicit
  virtual realT Assemble ( PhysicalDomainT& domain,
                           Epetra_System& epetraSystem,
                           const realT& time,
                           const realT& dt )
  {
    (void)domain;
    (void)epetraSystem;
    (void)time;
    (void)dt;

    throw GPException("ParallelPlateFlowSolverBase::Assemble() not overridden");
    return 0;
  }


  virtual void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt)
  {
    (void)domain;
    (void)partition;
    (void)time;
    (void)dt;

    throw GPException("ParallelPlateFlowSolverBase::Assemble() not overridden");
  }

  virtual void RegisterTemporaryFields( PhysicalDomainT& domain )
  {
    (void)domain;

//    throw GPException("ParallelPlateFlowSolverBase::RegisterTemporaryFields() not overridden");
  }

  virtual void DeregisterTemporaryFields( PhysicalDomainT& domain )
  {
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::DeregisterTemporaryFields() not overridden");
  }

  virtual void FillTemporaryFields( PhysicalDomainT& domain )
  {
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::FillTemporaryFields() not overridden");
  }

  virtual void OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain )
  {
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::OverwriteFieldsWithTemporaryFields() not overridden");
  }



  // functions for explicit
  virtual void UpdateFlux( const realT time,
                           const realT dt,
                           PhysicalDomainT& domain,
                           SpatialPartition& partition)
  {
    (void)time;
    (void)dt;
    (void)domain;
    (void)partition;

    throw GPException("ParallelPlateFlowSolverBase::UpdateFlux() not overridden");
  }

  // functions for proppant
  virtual void GenerateSlurryParameters( PhysicalDomainT& domain)
  {
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::GenerateSlurryParameters() not overridden");
  }


  virtual void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt )
  {
    (void)domain;
    (void)partition;
    (void)time;
    (void)dt;

	    throw GPException("ParallelPlateFlowSolverBase::CalculateMassRate() not overridden");
  }


  // boundary conditions
  virtual void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain )
  {
    (void)dt;
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::CalculateAndApplyMassFlux() not overridden");
  }


  virtual void CalculateCarterLeakOff( const realT time, const realT dt, PhysicalDomainT& domain )
  {
    (void)time;
    (void)dt;
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::CalculateCarterLeakOff() not overridden");
  }

  virtual void CalculateMatrixFlowLeakOff( const realT time, const realT dt, PhysicalDomainT& domain )
  {
    (void)time;
    (void)dt;
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::CalculateMatrixFlowLeakOff() not overridden");
  }

  virtual void ApplyFluxBoundaryCondition( const realT time, const realT dt, const int cycleNumber, const int rank, PhysicalDomainT& domain )
  {
    (void)time;
    (void)dt;
    (void)cycleNumber;
    (void)rank;
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::ApplyFluxBoundaryCondition() not overridden");
  }

  virtual void FlowControlledBoundaryCondition( PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
          BoundaryConditionBase* bc ,
          const lSet& set,
          realT time )
  {
    (void)domain;
    (void)object;
    (void)bc;
    (void)set;
    (void)time;
    throw GPException("ParallelPlateFlowSolverBase::FlowControlledBoundaryCondition() not overridden");
  }


  virtual void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
          realT time,
          realT dt )
  {
    (void)domain;
    (void)time;
    (void)dt;
    throw GPException("ParallelPlateFlowSolverBase::GenerateParallelPlateGeometricQuantities() not overridden");
  }

  virtual void OverwriteOldGeometricQuantities( PhysicalDomainT& domain )
  {
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::OverwriteOldGeometricQuantities() not overridden");
  }



  virtual void DefineFlowSets( PhysicalDomainT& domain )
  {
    (void)domain;
    throw GPException("ParallelPlateFlowSolverBase::DefineFlowSets() not overridden");
  }


  virtual void CalculateNodalPressure ( PhysicalDomainT& domain, SpatialPartition& partition)
  {
    (void)domain;
    (void)partition;
    throw GPException("ParallelPlateFlowSolverBase::CalculateNodalPressure() not overridden");
  }

  using SolverBase::SetInitialGuess;
  virtual void SetInitialGuess( const PhysicalDomainT& domain,
                                realT* const local_solution  )
  {
    (void)domain;
    (void)local_solution;
    throw GPException("ParallelPlateFlowSolverBase::SetInitialGuess() not overridden");
  }

  virtual realT TwoFacePermeability( const Array1dT<R1Tensor>& edgeCenters,
                                     const rArray1d& edgeLengths,
                                     const Array1dT<R1Tensor>& faceCenters,
                                     const rArray1d& apertures,
                                     const rArray1d& packVfs,
                                     realT mu,
                                     localIndex eg, localIndex kf, localIndex kfb)
  {
    (void)edgeCenters;
    (void)edgeLengths;
    (void)faceCenters;
    (void)apertures;
    (void)packVfs;
    (void)mu;
    (void)eg;
    (void)kf;
    (void)kfb;

    throw GPException("ParallelPlateFlowSolverBase::TwoFacePermeability() not overridden");
    realT kappa = 0.0;
    return kappa;
  }

  // virtual realT TwoFacePermeability( const Array1dT<R1Tensor>& edgeCenters,
  //                                    const rArray1d& edgeLengths,
  //                                    const Array1dT<R1Tensor>& faceCenters,
  //                                    const rArray1d& apertures,
  //                                    const localIndex eg,
  //                                    const localIndex r,
  //                                    const localIndex s,
  //                                    const Array1dT<rArray1d>* const dwdu,
  //                                    rArray1d* const dkdu_r,
  //                                    rArray1d* const dkdu_s )
  // {
  //   (void)edgeCenters;
  //   (void)edgeLengths;
  //   (void)faceCenters;
  //   (void)apertures;
  //   (void)eg;
  //   (void)r;
  //   (void)s;
  //   (void)dwdu;
  //   (void)dkdu_r;
  //   (void)dkdu_s;

  //   throw GPException("ParallelPlateFlowSolverBase::TwoFacePermeability() not overridden");
  //   realT kappa = 0.0;
  //   return kappa;
  // }
  
  virtual realT TwoFacePermeability_PowerLaw(const Array1dT<R1Tensor>& edgeCenters,
                                             const rArray1d& edgeLengths,
                                             const Array1dT<R1Tensor>& faceCenters,
                                             const rArray1d& apertures,
                                             const Array1dT<R1Tensor>& fluidVelocity,
                                             const rArray1d& packVfs,
                                             realT mu,
                                             localIndex eg, localIndex kf, localIndex kfb )
  {
    (void)edgeCenters;
    (void)edgeLengths;
    (void)faceCenters;
    (void)apertures;
    (void)packVfs;
    (void)mu;
    (void)eg;
    (void)kf;
    (void)kfb;

    throw GPException("ParallelPlateFlowSolverBase::TwoFacePermeability_PowerLaw() not overridden");
    realT kappa = 0.0;
    return kappa;
  }

  int m_boundPhysicalAperture;
  Table<1,realT> const * const m_ApertureTable;

  // Synced fields
  // std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> m_syncedFieldsB;

  // Multiphase flow
  sArray1d m_tracerNames, m_tracerNamesMass, m_tracerNamesDMDT, m_tracerNamesVolumeFraction;
  rArrayPtrs m_tracerMassPtrs, m_tracerDMDTPtrs, m_volumeFractionPtrs;
  bool m_multiphaseFlow, m_proppantActive;
  int m_multiphaseMixingMode;
  realT m_multiphaseMixingThreshold_serialFlow, m_multiphaseMixingThreshold_singlePhase;
  sArray1d m_multiphaseSetnames;

  Array1dT< rArray1d > m_dwdu;
  Array1dT< iArray1d > m_dwdu_dof;
  rArray1d m_dwdw;

  rArray1d m_multiphaseDensity, m_multiphaseBulkModulus, m_multiphaseViscosity;

  

  std::map<localIndex,lArray1d> m_edgesToFaces;
  R1Tensor m_gravityVector;

  LeakoffBase m_leakoffModel;
  realT m_leakoffCoef = 0.0;
  int m_overLeakCompensation;
  // If =1, we allow flow cells to temporarily have a negative fluid mass.

  realT zeroApertureOffset() const
  { return m_zeroApertureOffset; }

  realT zeroApertureVolume() const
  { return m_zeroApertureVolume; }

  // Barton joint model
  realT m_bBarton;
  realT m_aBarton;
  realT m_wZeroStress;
  realT m_kContact;

  // upscaling variables
  realT m_subscaleAngle = 0;

protected:
  void MarkFaceSaturationTime( PhysicalDomainT& domain, const realT time,const realT dt );
  realT m_maxLeakOffRatio;


public:

  // permeability calculation
  realT m_SHP_FCT;
  realT m_mu;
  realT m_min_aperture;
  realT m_max_aperture;
  realT m_zeroApertureVolume;
  realT m_zeroApertureOffset;

//  realT m_min_hydro_aperture;
  // I am adding this for verification against high leak solutions.
  // For those situations, we need to use a small zero-stress barton aperture.  Otherwise, the transition near zero stress introduce some unacceptable error because the overall aperture is so small.
  // But we still want some some good permeability on the closed faces.  Otherwise if we use the small barton parameters, it will take forever to fill the closed fracture.
  // We will see how it goes.

  realT m_apertureTransition = 0.0;
  realT m_apertureSlopeReduction = 1.0;
  int m_permeabilityUpdateOption = 0;

  // fluid eos
  realT m_bulk_modulus;
  realT m_rho_o; // fluid density @ pressure=0 (for eos)
  realT m_press_o;

  PPFS::FluidEOSBase * m_fluidEOS;

  int m_updateFaceArea;

  // powerlaw fluid
  bool m_usePowerlawFluid;  // ultimately want a generic fluid model here
  //PowerlawFluidModel m_powerlawFluid;
  ParallelPlateFluidModelBase* m_fluidModelPtr;

  // This is the maximum pressure allowed in the flow solver.  The purpose of having this cap is to prevent sudden pressure spikes.
  // If the pressure value according to EOS is lower than half of this limiting value, the EOS will be directly used.  Otherwise, some special treatment will be adopted.
  // Its value should be significantly lower than the bulk modulus.
  // It should be at least twice of the minimum principal stress in the rock matrix.
  realT m_pressureCap;
  realT m_maxViscosity;

  std::string m_flowFaceSetName;
  std::string m_wellboreSolverName;

  lSet m_negativePressureSet;
  realT m_cutoffPressure = -1.0e99;
  realT m_negativeFluidPressureSlopeReduction = 1.0e-3;


  realT m_dT;
  realT m_farFieldPorePressure;
  int m_pressureDependentLeakoff;


  enum class dofVariable
  {
    mass,
    pressure
  };

  dofVariable m_dofVariable;

  realT m_contactOffsetCutoff = -1.0e99;


  void RegisterFields( FaceManagerT& faceManager, EdgeManagerT& edgeManager );







private:
  ParallelPlateFlowSolverBase();
};


inline realT ParallelPlateFlowSolverBase::BoundedAperture( const realT& aper )
{
  realT rval = aper;
  if( this->m_boundPhysicalAperture == 3 )
  {
    realT b = m_bBarton;
    realT a = m_aBarton;
    realT k = m_kContact;

    realT rootTerm = 1.0 - 4*a*k;
    realT intersectionPoint = 1.0e99;
    if( rootTerm > 0.0  )
    {
      intersectionPoint = ( 1.0 - sqrt(rootTerm) ) / ( 2.0 * b * k );
    }

    if( aper < intersectionPoint )
    {
      if( aper <= ( 1.0 - sqrt(a*k) ) / (b*k) )
      {
        rval = a/b * ( 1.0/( 1.0 + b * (-k * aper) ) );
      }
      else
      {
        rval += (-1.0 + 2.0*sqrt(a*k))/(b*k);
      }
    }
  }
  else if( this->m_boundPhysicalAperture == 2 )
  {
    realT input[1] = {aper};
    rval = m_ApertureTable->Lookup( input );
  }
  else
  {
  if( m_apertureTransition > 0.0 )
  {
    if( aper<m_apertureTransition )
    {
      rval =  m_apertureTransition - ( m_apertureTransition - aper ) * m_apertureSlopeReduction;
    }

    if( rval < this->m_min_aperture )
    {
      rval = this->m_min_aperture;
    }
  }
  else
  {
    if( aper<m_min_aperture )
    {
      rval = m_min_aperture;
    }
    else if( aper>m_max_aperture )
    {
      rval = m_max_aperture;
    }

  }
  }
  return rval;
}

inline realT ParallelPlateFlowSolverBase::BoundedApertureDerivative( const realT& aper )
{
  realT rval = 1;

  if( this->m_boundPhysicalAperture == 3 )
  {
    realT b = m_bBarton;
    realT a = m_aBarton;
    realT k = m_kContact;
    realT rootTerm = 1.0 - 4*a*k;
    realT intersectionPoint = 1.0e99;
    if( rootTerm > 0.0  )
    {
      intersectionPoint = ( 1.0 - sqrt(rootTerm) ) / ( 2.0 * b * k );
    }

    if( aper < intersectionPoint )
    {
      if( aper <= ( 1.0 - sqrt(a*k) ) / (b*k) )
      {
        rval = ( a * k ) / pow( ( 1.0 + b * (-k * aper) ), 2 );
      }
    }
  }
  else if( this->m_boundPhysicalAperture == 2 )
  {
    realT input[1] = {aper};
    rval = m_ApertureTable->Gradient( input );
  }
  else
  {
    if( m_apertureTransition > 0.0 )
    {
      if( aper<m_apertureTransition )
      {
        rval = m_apertureSlopeReduction;
      }
    }
    else
    {
      if( aper<m_min_aperture )
      {
        rval = 0 ;
      }
      else if( aper>m_max_aperture )
      {
        rval = 0 ;
      }
    }
  }

  return rval;
}
// Structure for sorting multiphase fields

struct multiphaseSort 
{ 
  realT mass;
  localIndex phase;
};

struct by_mass
{
  bool operator()(multiphaseSort const &a, multiphaseSort const &b)
  {
    return a.mass > b.mass;
  }
};




namespace PPFS
{
  // Equation of state for pressure
  inline realT P_EOS(realT rho, realT K_o, realT rho_o)
  {
    return K_o*(rho - rho_o)/rho;
  }

  inline realT Inverse_EOS(realT P, realT K_o, realT rho_o)
  {
    return  rho_o /( 1.0 - P / K_o );
  }

  // With pressure cap.
  inline realT P_EOS(realT rho, realT K_o, realT rho_o, realT pressureCap)
  {
    realT P;
    if (rho < rho_o)
    {
      P = 0.0;
    }
    else
    {
      P = K_o*(rho - rho_o)/rho;
      if ( P > 0.5 * pressureCap)
      {
        P = 0.5 * pressureCap + ( P - 0.5 * pressureCap) / (K_o - 0.5 * pressureCap) * 0.5 * pressureCap;
      }

    }
    return (P);
  }

  inline realT Inverse_EOS(realT P, realT K_o, realT rho_o, realT pressureCap)
  {
    realT rho;
    if (P <= pressureCap * 0.5)
    {
      rho = K_o * rho_o / (K_o - P);
    }
    else
    {
      rho = K_o * rho_o / (K_o - 0.5*pressureCap - (2.0*P/pressureCap - 1.0) * (K_o - 0.5*pressureCap));
    }
    return(rho);
  }


  // Equation for d Pressure/ d rho
  inline realT dPdRho_EOS(realT K_o, realT rho, realT rho_o)
  {
    return K_o*rho_o/(rho*rho);
  }

  // with pressure cap
  inline realT dPdRho_EOS(realT rho, realT K_o, realT rho_o, realT pressureCap)
  {
	 realT P = K_o*(rho - rho_o)/rho;
	 realT dPdRho(0);
     if ( P < 0.5 * pressureCap)
     {
	   dPdRho = K_o*rho_o/(rho*rho);
     } else {
		 dPdRho *= 0.5 * pressureCap/ (K_o - 0.5 * pressureCap);
	 }

	 return dPdRho;
  }

  // Equation for K: K = rho* d Pressure/ d rho
  inline realT BulkModulus_EOS(realT rho,realT K_o, realT rho_o)
  {
    return rho*dPdRho_EOS(K_o, rho, rho_o);
  }

  // Equation of state for density
  inline realT rho_EOS(realT P,realT K_o, realT rho_o)
  {
    return (1.0 + P/K_o)*rho_o;
  }

  // with pressure cap
  inline realT rho_EOS(realT P, realT K_o, realT rho_o, realT pressureCap)
  {
	    realT rho;
	    if (P <= pressureCap * 0.5)
	    {
	      rho = K_o * rho_o / (K_o - P);
	    }
	    else
	    {
	      rho = K_o * rho_o / (K_o - 0.5*pressureCap - (2.0*P/pressureCap - 1.0) * (K_o - 0.5*pressureCap));
	    }
	    return(rho);
  }
  /**
   *
   * @param la length of the path between centers of two flow elements projected on normal attached to element a
   * @param lb length of the path between centers of two flow elements projected on normal attached to element b
   * @param apa aperture of the first flow element
   * @param apb aperture of the second flow element
   * @param w width of the flow path
   *
   * The "Permeability" returned is scaled by the cell face area
   * and divided by the product of the distance between the cells and kinematic viscosity.
   *
   *   k' = k*(h*w)/(mu*L)
   *
   * so that the mass flux from cell a to cell b is given by
   *
   *   q_Mass = k' * rho * (Pa-Pb)
   *
   */
  inline realT CalculatePermeability(const realT la,
                                     const realT lb,
                                     const realT apa,
                                     const realT apb,
                                     const realT w,
                                     const realT mu,
                                     const realT SHP_FCT)
  {
    if (apa <= 0 || apb <= 0)
      return 0.0;


    const realT ka = apa*apa*apa;
    const realT kb = apb*apb*apb;

    realT permeability = ka * kb * w/(12.0 * mu *(ka*lb + kb*la)); // harmonic mean

    permeability *= SHP_FCT;
    return permeability;
  }

  /**
   * One sided permeability calculation
   *
   * @author walsh24
   *
   * @param la length of the path between centers of two flow elements projected on normal attached to element a\
   * @param h aperture of the flow element
   * @param w width of the flow path
   *
   * NB formulation differs from that in Johnson Morris 2009
   * which uses distance between adjacent elements (2*la) rather than element and boundary (la).
   *
   * The "Permeability" returned is scaled by the cell face area
   * and divided by the product of the distance between the cells and kinematic viscosity.
   *
   *   k' = k*(h*w)/(mu*L)
   *
   * so that the mass flux from cell a to cell b is given by
   *
   *   q_Mass = k' * rho * (Pa-Pb)
   *
   */
  inline realT CalculatePermeability(const realT la,
                                     const realT h,
                                     const realT w,
                                     const realT mu,
                                     const realT SHP_FCT)
  {
    if (h <= 0) return 0.0;

    realT permeability = h*h*h *w / ( 12.0 * mu * la );

    permeability *= SHP_FCT;
    return permeability;
  }

  inline void CalculatePermeabilityAndDerivative( const realT length,
                                                  const realT aper,
                                                  const realT aper0,
                                                  const realT integrationOption,
                                                  const realT width,
                                                  const realT mu,
                                                  realT& permeability,
                                                  realT& dKappa_dAper )
  {

    const realT factor = width / ( 12.0 * mu * length );
    // forward euler
    if( integrationOption == 0 )
    {
      permeability = aper0*aper0*aper0 * factor;
      dKappa_dAper = 0.0;
    }
    // backward euler
    else if( integrationOption == 1 )
    {
      permeability = aper*aper*aper * factor;
      dKappa_dAper = 3 * aper*aper* factor;
    }
    // trapazoid rule
    else if ( integrationOption == 2 )
    {
      permeability = 0.5 * ( aper0*aper0*aper0 + aper*aper*aper ) * factor;
      dKappa_dAper = 1.5 * aper*aper * factor;
    }
    // simpsons rule
    else if ( integrationOption == 3 )
    {
      permeability = 0.25 * ( aper0*aper0*aper0
                            + aper0*aper0*aper
                            + aper0*aper*aper
                            + aper*aper*aper) * factor;
      dKappa_dAper = 0.25 * ( aper0*aper0
                            + 2*aper0*aper
                            + 3*aper*aper) * factor;
    }
    else
    {
      throw GPException("ParallelPlateFlowSolverBase.h::CalculatePermeabilityAndDerivative() invalid aperPow");
    }

    return;
  }

  inline realT CalculatePRhoGravity(const R1Tensor xa,
                                     const R1Tensor xb,
                                     const realT rho_a,
                                     const realT rho_b,
                                     const R1Tensor g)
  {
    R1Tensor vec = xb;
    vec -= xa;
    realT dist = Dot(vec, g);
    realT Prho = dist * (rho_a + rho_b) * 0.5;
    Prho *= (rho_a + rho_b) * 0.5;
    return Prho;
  }

  inline realT CalculatePRhoGravity(const R1Tensor xa,
                                     const R1Tensor xb,
                                     const realT rho_a,
                                     const R1Tensor g)
  {
    R1Tensor vec = xb;
    vec -= xa;
    realT dist = Dot(vec, g);
    realT Prho = dist * rho_a * rho_a;
    return Prho;
  }

  /**
     *
     * @param la length of the path between centers of two flow elements projected on normal attached to element a
     * @param lb length of the path between centers of two flow elements projected on normal attached to element b
     * @param apa aperture of the first flow element
     * @param apb aperture of the second flow element
     * @param w width of the flow path
     * @param phiM = M' = modified consistency index     (units Pa.s^n)
     * @param qMag magnitude of fluid velocity
     * @param n fluid behavior index  (dimensionless)
     *
     * Power Law fluid: \sigma_{xy} = M(2 \dot{\epsilon}_{xy})^{n}
     * n = 1 : Newtonian fluid
     * n < 1 : Shear thinning
     * n > 1 : Shear thickening
     *
     * For a powerlaw fluid between two parallel plates:
     * dp/dx_{i} = - \phi M q_i |q|^{n-1}/h^{2n+1}
     * Where \phi = 2^{n+1}(2*n+1)^n/n^n
     * See Adachi and Detournay,
     * Self-similar solution of a plane-strain fracture driven by a power-law fluid (2002)
     *
     * The "Permeability" returned is scaled by the cell face area
     * and divided by the product of the distance between the cells and kinematic viscosity.
     *
     *   k' = k*(h*w)/(mu*L)
     *
     * so that the mass flux from cell a to cell b is given by
     *
     *   q_Mass = k' * rho * (Pa-Pb)
     *
     */
    inline realT CalculatePermeability_PowerLawFluid(const realT la,
                                                     const realT lb,
                                                     const realT apa,
                                                     const realT apb,
                                                     const realT w,
                                                     const realT qMag,
                                                     const realT phiM,
                                                     const realT n,
                                                     const realT SHP_FCT)
    {
      if (apa <= 0 || apb <= 0)
        return 0.0;


      const realT ka = pow(apa,2*n+1);
      const realT kb = pow(apb,2*n+1);

      realT permeability = ka * kb * w*pow(qMag,1-n)/(phiM *(ka*lb + kb*la)); // harmonic mean

      permeability *= SHP_FCT;
      return permeability;
    }

    // one sided permeability calculation for powerlaw fluid
    inline realT CalculatePermeability_PowerLawFluid(const realT la,
                                                     const realT h,
                                                     const realT w,
                                                     const realT qMag,
                                                     const realT phiM,
                                                     const realT n,
                                                     const realT SHP_FCT)
    {
      if (h <= 0) return 0.0;

      realT permeability = pow(h,2*n+1) * w*pow(qMag,1-n)/(phiM*la);

      permeability *= SHP_FCT;
      return permeability;
    }

    // modified consistency index is required when calculating the power law fluid permeability
    inline realT CalculateModifiedConsistencyIndex(const realT M, const realT n){
    	realT phi =2*pow(2*(2*n+1)/n,n);
    	return phi*M;
    }

    /**
        *
        * @param la length of the path between centers of two flow elements projected on normal attached to element a
        * @param lb length of the path between centers of two flow elements projected on normal attached to element b
        * @param apa aperture of the first flow element
        * @param apb aperture of the second flow element
        * @param w width of the flow path
        * @param phiM = M' = modified consistency index     (units Pa.s^n)
        * @param qMag magnitude of fluid velocity
        * @param n fluid behavior index  (dimensionless)
        * @param z_p dimensionless plug width = |2 \tau_y|/(|p_{,i}|h) where h = aperture, \tau_y is the yield stress and |p_{,i}| is the magnitude of the pressure gradient
        *
        * Herschel-Bulkley fluid:
        *    \tau = \tau_y + M(2 \dot{\epsilon}_{xy})^{n}  for \tau > \tau_y
        *    \dot{\epsilon}_{xy} = 0 for \tau < \tau_y
        *
        * See Wang and Gordaninejad 1999 "Flow analysis of Field controllable Electro and Magneto-Rheological fluids using Herschel Bulkley Model",
        * Self-similar solution of a plane-strain fracture driven by a power-law fluid (2002)
        *
        * The "Permeability" returned is scaled by the cell face area
        * and divided by the product of the distance between the cells and kinematic viscosity.
        *
        *   k' = k*(h*w)/(mu*L)
        *
        * so that the mass flux from cell a to cell b is given by
        *
        *   q_Mass = k' * rho * (Pa-Pb)
        *
        */
       inline realT CalculatePermeability_HerschelBulkleyFluid(const realT la,
                                                        const realT lb,
                                                        const realT apa,
                                                        const realT apb,
                                                        const realT w,
                                                        const realT qMag,
                                                        const realT phiM,
                                                        const realT n,
                                                        const realT zp,
                                                        const realT SHP_FCT)
       {
         if (apa <= 0 || apb <= 0 || zp >= 1)
           return 0.0;

         realT permeability =CalculatePermeability_PowerLawFluid(la, lb,apa, apb,w, qMag, phiM, n, SHP_FCT);
         realT denom = pow(1-zp , n+1)*pow(n/(n+1)*zp +1,n);

         permeability /= denom;

         return permeability;
       }

       // one sided permeability calculation for Herschel-Bulkley fluid
       inline realT CalculatePermeability_HerschelBulkleyFluid(const realT la,
                                                        const realT h,
                                                        const realT w,
                                                        const realT qMag,
                                                        const realT phiM,
                                                        const realT n,
                                                        const realT zp,
                                                        const realT SHP_FCT)
       {
         if (h <= 0|| zp >= 1) return 0.0;


         realT permeability =CalculatePermeability_PowerLawFluid(la,h,w, qMag, phiM, n, SHP_FCT);
         realT denom = pow(1-zp , n+1)*pow(n/(n+1)*zp +1,n);

         permeability /= denom;

         return permeability;
       }




       ////////////////////////
       // Equation of state  //
       ////////////////////////

       // base class for fluid equation of state
       class FluidEOSBase{
    	   public:

         FluidEOSBase( TICPP::HierarchicalDataNode* hdn)
    	   { (void)hdn; }
    	     FluidEOSBase(realT rho):m_rho_o(rho){/*empty*/}
             virtual ~FluidEOSBase(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             virtual realT pressure(realT rho) = 0;
             virtual realT density(realT P) = 0;
             virtual realT dPdRho(realT rho) = 0;
             virtual realT dPdRho(realT rho, realT rho0) { return dPdRho(rho);}
             virtual realT rho_o(){return m_rho_o;}

             virtual realT pressure_mp(rArray1d volumeFraction, realT rho){return pressure(rho);}
             virtual realT density_mp(rArray1d volumeFraction, realT P){return density(P);}
             virtual realT dPdRho_mp(rArray1d volumeFraction, realT rho){return dPdRho(rho);}
             virtual realT rho_o_mp(rArray1d volumeFraction){return rho_o();}
             virtual realT viscosity_mp(rArray1d volumeFraction){return 0.0;}    	   
             virtual realT PressureAndApertureCloseJoint_mp( rArray1d volumeFraction, realT area, realT fluidMass,
                                                    realT totalStress, realT a, realT b, realT &aperture )
             {
               // We need to calculate both pressure and aperture based on the fluid mass and total stress.
               // The pressure-aperture feedback does not work well for closed fractures as we will have severe over-shoot/under-shoot problems.
               // See Fu and Carrigan, Modeling responss of naturally fractured geothermal reservoir to low pressure stimulation, GRC 2012 for derivation of the formulations.
               throw GPException( "PressureAndApertureCloseJoint_mp not overridden for this fluid EOS type" );
               return 0;
             }
            protected:
             realT m_rho_o;
       };

       // Linear EOS
       class LinearEOS: public FluidEOSBase{
    	   public:

    	     LinearEOS(TICPP::HierarchicalDataNode* hdn);
    	     LinearEOS(realT rho, realT Ko): FluidEOSBase(rho), m_K_o(Ko){ /* empty */};
             virtual ~LinearEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             virtual realT pressure(realT rho){return m_K_o*(rho - m_rho_o)/m_rho_o;};
             virtual realT density(realT P){return (P + m_K_o) * m_rho_o / m_K_o; };
             using FluidEOSBase::dPdRho;
             virtual realT dPdRho(realT rho)
             {
               (void)rho;
               return m_K_o/m_rho_o;
             };

             static const char* FluidEOSName(){return "LinearEOS";};
    	   protected:
             realT m_K_o;
       };

       // PressureCap EOS
       // Per Randy's suggestion, the cap EOS is simplified to use a simple cutoff cap. The old one uses a curve that approaches the pressure cap at infinite density.
       class PressureCapEOS: public LinearEOS{
    	   public:

    	     PressureCapEOS(TICPP::HierarchicalDataNode* hdn);
    	     PressureCapEOS(realT rho, realT Ko, realT pCap):LinearEOS(rho,Ko), m_pressureCap(pCap){/* empty */};
             virtual ~PressureCapEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             static const char* FluidEOSName(){return "PressureCapEOS";};

             virtual realT pressure(realT rho){
            	    realT P;
            	    if (rho < m_rho_o)
            	    {
            	      P = 0.0;
            	    }
            	    else
            	    {
                    P = std::min(m_K_o*(rho - m_rho_o)/rho, m_pressureCap);
            	    }
            	    return (P);
             };

             virtual realT density(realT P){
            	    realT rho;
            	    if (P <= m_pressureCap)
            	    {
            	      rho = m_K_o * m_rho_o / (m_K_o - P);
            	    }
            	    else
            	    {
            	      rho = m_K_o * m_rho_o / (m_K_o - m_pressureCap);
            	      // Since this function is only called by boundary conditions, the possibility of losing mass should be low.
            	    }

            	    return(rho);
             };

             using FluidEOSBase::dPdRho;
             virtual realT dPdRho(realT rho){
            	 realT rval(0.0);
            	 if(rho <= m_rho_o){
            		// return 0.0; // getting some floating point exceptions - not sure if this is the best fix
            	 } else {

            	     realT P = m_K_o*(rho - m_rho_o)/(rho);
            	     if ( P < m_pressureCap)
            	     {
            	          rval = m_K_o*m_rho_o/(rho*rho);
            	     } else {
           	            rval = 0;
            	     }
            	 }
         	     return rval;
             };

             virtual realT PressureAndApertureCloseJoint_mp( rArray1d volumeFraction, realT area, realT fluidMass,
                                                    realT totalStress, realT a, realT b, realT &aperture )
             {

               realT w0, wDry, A, B, rho;
               w0 = a/b;
               wDry = w0 - a * totalStress / (1 + b * totalStress);
               if ( fluidMass < wDry * area *m_rho_o)
               {
                 aperture = wDry;
                 return 0.0;
               }
               else
               {
                 A = m_K_o * m_rho_o * area / fluidMass;
                 B = totalStress - m_K_o + A * w0;
                 aperture = w0 - (A*a + B*b + 1 - pow( (A*a+B*b+1)*(A*a+B*b+1) - 4 * A*a*B*b,0.5))/2/A/b;
                 rho = fluidMass / area / aperture;
               }
               return pressure(rho);

             }



    	   protected:
             realT m_pressureCap;
       };


       // Pressure EOS
       class PressureEOS: public LinearEOS{
         public:

           PressureEOS(TICPP::HierarchicalDataNode* hdn);
           PressureEOS(realT rho, realT Ko):LinearEOS(rho,Ko){/* empty */};
             virtual ~PressureEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             static const char* FluidEOSName(){return "PressureEOS";};

             virtual realT pressure(realT rho){
                  realT P;
                  if (rho < m_rho_o)
                  {
                    P = 0.0;
                  }
                  else
                  {
                    P = m_K_o*(rho - m_rho_o)/rho;
                  }
                  return (P);
             };

             virtual realT density(realT P){
                  realT rho;
                  if (P <= 0.99*m_K_o)
                  {
                    rho = m_K_o * m_rho_o / (m_K_o - P);
                  }
                  else
                  {
                    rho = m_K_o * m_rho_o / (0.01*m_K_o);
                    // Since this function is only called by boundary conditions, the possibility of losing mass should be low.
                  }
                  return(rho);
             };

             using FluidEOSBase::dPdRho;
             virtual realT dPdRho(realT rho){
               realT rval(0.0);
               if(rho <= m_rho_o){
                // return 0.0; // getting some floating point exceptions - not sure if this is the best fix
               } 
               else 
               {
                rval = m_K_o*m_rho_o/(rho*rho);
               }
               return rval;
             };
       };

       // Pressure EOS
       class BiLinearEOS: public LinearEOS{
       public:

         BiLinearEOS(TICPP::HierarchicalDataNode* hdn);
         BiLinearEOS(realT rho, realT Ko, realT K_neg):
           LinearEOS(rho, Ko),
           m_K_neg(K_neg)
         {}

         virtual ~BiLinearEOS(){}

         virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

         virtual realT pressure(realT rho)
         {
           realT p=0;
           if( rho>=m_rho_o )
           {
             p = m_K_o*(rho/m_rho_o - 1.0);
           }
           else
           {
             p = m_K_neg*(rho/m_rho_o - 1.0);
           }
           return p;
         };

         virtual realT density(realT P)
         {
           realT rho=0;
           if( P>=0 )
           {
             rho = ( P / m_K_o + 1.0 ) * m_rho_o;
           }
           else
           {
             rho = ( P / m_K_neg + 1.0 ) * m_rho_o;
           }
           return rho;
         };

         using FluidEOSBase::dPdRho;
         virtual realT dPdRho(realT rho)
         {
           realT dpdrho = 0;
           if( rho>=m_rho_o )
           {
             dpdrho = m_K_o/m_rho_o;
           }
           else
           {
             dpdrho = m_K_neg/m_rho_o;
           }

           return dpdrho;
         };

         virtual realT dPdRho(realT rho, realT rhoLast)
         {
           realT rval = 0.0;
           if( fabs ( rho - rhoLast )/ rho < 1.0e-6 )
           {
             rval = dPdRho(rho);
           }
           else
           {
             rval = ( pressure(rho) - pressure(rhoLast) ) / ( rho - rhoLast );
           }
           return rval;
         }


         static const char* FluidEOSName(){return "BiLinearEOS";};
       protected:
             realT m_K_neg;
       };


       class BiNonLinearEOS: public BiLinearEOS{
       public:

         BiNonLinearEOS(TICPP::HierarchicalDataNode* hdn);
         BiNonLinearEOS(realT rho, realT Ko, realT K_neg):
           BiLinearEOS(rho, Ko, K_neg)
         {}

         virtual ~BiNonLinearEOS(){}

         virtual realT pressure(realT rho)
         {
           realT p=0;
           if( rho>=m_rho_o )
           {
             p = m_K_o*( 1.0 - m_rho_o/rho );
           }
           else
           {
             p = m_K_neg*( 1.0 - m_rho_o/rho );
           }
           return p;
         };

         virtual realT density(realT P)
         {
           realT rho=0;
           if( P>=0 )
           {
             rho = ( m_K_o * m_rho_o ) / ( m_K_o - P );
           }
           else
           {
             rho = ( m_K_neg * m_rho_o ) / ( m_K_o - P );
           }
           return rho;
         };

         using FluidEOSBase::dPdRho;
         virtual realT dPdRho(realT rho)
         {
           realT dpdrho = 0;
           if( rho>=m_rho_o )
           {
             dpdrho = m_K_o * m_rho_o / ( rho*rho );
           }
           else
           {
             dpdrho = m_K_neg * m_rho_o / ( rho*rho );
           }

           return dpdrho;
         };

         virtual realT dPdRho(realT rho, realT rhoLast)
         {
           realT rval = 0.0;
           if( fabs ( rho - rhoLast )/ rho < 1.0e-6 )
           {
             rval = dPdRho(rho);
           }
           else
           {
             rval = ( pressure(rho) - pressure(rhoLast) ) / ( rho - rhoLast );
           }
           return rval;
         }


         static const char* FluidEOSName(){return "BiNonLinearEOS";};
       };


       class MultiphasePressureEOS: public LinearEOS{
         public:

           MultiphasePressureEOS(TICPP::HierarchicalDataNode* hdn);
           MultiphasePressureEOS(rArray1d multiphaseDensity, rArray1d multiphaseBulkModulus, rArray1d multiphaseViscosity): LinearEOS(multiphaseDensity[0], multiphaseBulkModulus[0])
           {
              m_NFluids = multiphaseDensity.size();

              m_multiphaseDensity.resize(m_NFluids);
              m_multiphaseDensity = multiphaseDensity;

              m_multiphaseBulkModulus.resize(m_NFluids);
              m_multiphaseBulkModulus = multiphaseBulkModulus;
              
              m_multiphaseViscosity.resize(m_NFluids);
              m_multiphaseViscosity = multiphaseViscosity;
           };

           virtual ~MultiphasePressureEOS(){};
           virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

           static const char* FluidEOSName(){return "MultiphasePressureEOS";};

        
          virtual realT viscosity_mp(rArray1d volumeFraction)
          {
            realT tmpViscosity = 0.0;
            for (localIndex ii=0; ii<m_NFluids; ii++)
            {
              tmpViscosity += volumeFraction[ii]*m_multiphaseViscosity[ii];
            }

            return tmpViscosity;
          };

          virtual realT rho_o_mp(rArray1d volumeFraction)
          {
            realT tmpRho = 0.0;
            for (localIndex ii=0; ii<m_NFluids; ii++)
            {
              tmpRho += volumeFraction[ii]*m_multiphaseDensity[ii];
            }

            return tmpRho;
          };


          virtual realT pressure_mp(rArray1d volumeFraction, realT rho)
          {
            realT P;
            realT tmpRho = 0.0;
            realT tmpK = 0.0;

            for (localIndex ii=0; ii<m_NFluids; ii++)
            {
              tmpRho += volumeFraction[ii]*m_multiphaseDensity[ii];
              tmpK += volumeFraction[ii]*m_multiphaseBulkModulus[ii];
            }
     
            if (rho < tmpRho)
            {
              P = 0.0;
            }
            else
            {
              P = tmpK*(rho - tmpRho)/rho;
            }
            return (P);
          };

          virtual realT PressureAndApertureCloseJoint_mp( rArray1d volumeFraction, realT area, realT fluidMass,
                                                 realT totalStress, realT a, realT b, realT &aperture )
          {
            realT P;
            realT tmpRho = 0.0;
            realT tmpK = 0.0;

            for (localIndex ii=0; ii<m_NFluids; ii++)
            {
              tmpRho += volumeFraction[ii]*m_multiphaseDensity[ii];
              tmpK += volumeFraction[ii]*m_multiphaseBulkModulus[ii];
            }

            realT w0, wDry, A, B, rho;
            w0 = a/b;
            wDry = w0 - a * totalStress / (1 + b * totalStress);
            if ( fluidMass < wDry * area *tmpRho)
            {
              aperture = wDry;
              return 0.0;
            }
            else
            {
              A = tmpK * tmpRho * area / fluidMass;
              B = totalStress - tmpK + A * w0;
              aperture = w0 - (A*a + B*b + 1 - pow( (A*a+B*b+1)*(A*a+B*b+1) - 4 * A*a*B*b,0.5))/2/A/b;
              rho = fluidMass / area / aperture;
            }
            P = tmpK*(rho - tmpRho)/rho;
            return P;

          }


          virtual realT density_mp(rArray1d volumeFraction, realT P)
          {
            realT tmpRho = 0.0;
            realT tmpK = 0.0;

            for (localIndex ii=0; ii<m_NFluids; ii++)
            {
              tmpRho += volumeFraction[ii]*m_multiphaseDensity[ii];
              tmpK += volumeFraction[ii]*m_multiphaseBulkModulus[ii];
            }

            realT rho;
            if (P <= 0.99*tmpK)
            {
              rho = tmpK * tmpRho / (tmpK - P);
            }
            else
            {
              rho = tmpK * tmpRho / (0.01*tmpK);
              // Since this function is only called by boundary conditions, the possibility of losing mass should be low.
            }
            return(rho);
          };


          virtual realT dPdRho_mp(rArray1d volumeFraction, realT rho)
          {
            realT tmpRho = 0.0;
            realT tmpK = 0.0;
            for (localIndex ii=0; ii<m_NFluids; ii++)
            {
              tmpRho += volumeFraction[ii]*m_multiphaseDensity[ii];
              tmpK += volumeFraction[ii]*m_multiphaseBulkModulus[ii];
            }

            realT rval(0.0);
            if(rho <= tmpRho){
             // return 0.0; // getting some floating point exceptions - not sure if this is the best fix
            } 
            else 
            {
              rval = tmpK*tmpRho/(rho*rho);
            }
            return rval;
          };

         localIndex m_NFluids;

         protected:
             rArray1d m_multiphaseDensity, m_multiphaseBulkModulus, m_multiphaseViscosity;
       };


       // Adiabatic EOS
       class AdiabaticEOS: public FluidEOSBase{
    	   public:

    	     AdiabaticEOS(TICPP::HierarchicalDataNode* hdn);
             virtual ~AdiabaticEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             virtual realT pressure(realT rho){
            	 realT rv = (rho > 0)? m_P_o*std::pow(rho/m_rho_o,m_gamma): 0.0;
            	 rv -= m_P_ref;
            	 if(rv < 0.0) rv = 0.0;
                 return rv;
             };
             virtual realT density(realT P){return (P+m_P_ref > 0)? m_rho_o*std::pow( (P+m_P_ref)/m_P_o, 1/m_gamma ): 0.0; };

             using FluidEOSBase::dPdRho;
             virtual realT dPdRho(realT rho){
            	 realT dPdRho_return = 0;
            	 //if(rho > 0) dPdRho = m_gamma*m_P_o * std::pow(rho/m_rho_o,m_gamma)/rho;
            	 realT p = pressure(rho);
            	 if(p > 0) dPdRho_return = m_gamma*m_P_o * std::pow(rho/m_rho_o,m_gamma)/rho;
            	 return dPdRho_return;
             };

             static const char* FluidEOSName(){return "AdiabaticEOS";};
    	   protected:
             realT m_P_o;
             realT m_gamma;
             realT m_P_ref;
       };


//////////////////////////

// Fluid EOS Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newFluidEOS"
//   * A base class to derive the functions to generate FluidEOS pointers: "FluidEOSInitializer"
//   * A String-to-FluidEOS-Intializer map hidden behind the getFluidEOSCatalogue function
//   * A template to create FluidEOS initializers: "FluidEOSRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_FLUID_EOS"
//
// Most initial conditions will only need to use one or two of the parts:
//   * To register a new FluidEOS in the factory: REGISTER_FluidEOS( FluidEOSClassName )
//   * To load a FluidEOS pointer from the factory:       FluidEOSBase* aFluidEOSPtr = newFluidEOS(FluidEOSString, args );

/// The FluidEOS Factory.
FluidEOSBase* newFluidEOS(const std::string& EOSName, TICPP::HierarchicalDataNode* hdn);

/// Base class to generate new FluidEOS pointers
class FluidEOSInitializer{
public:
	virtual FluidEOSBase* initializeFluidEOS( TICPP::HierarchicalDataNode* hdn) = 0;
	virtual ~FluidEOSInitializer() = 0;
};
inline FluidEOSInitializer::~FluidEOSInitializer() { }

/// Interface to the FluidEOS name -> FluidEOS initializer map
std::map<std::string, FluidEOSInitializer*> & getFluidEOSCatalogue();

/// Return a list of supported FluidEOS names
void getFluidEOSNames( std::vector<std::string>& nameList);

/// Template for creating classes derived from FluidEOSInitializer
template< class FluidEOSType >
class FluidEOSRegistrator : public FluidEOSInitializer{

public:
	FluidEOSRegistrator(void){
		std::string FluidEOSName = std::string(FluidEOSType::FluidEOSName() );
		getFluidEOSCatalogue() [FluidEOSName] = this;
	};

	FluidEOSBase* initializeFluidEOS( TICPP::HierarchicalDataNode* hdn ){
		return new FluidEOSType( hdn );
	};
};



} // end namespace

/// Compiler directive to simplify autoregistration
#define REGISTER_FLUID_EOS( ClassName ) namespace{ PPFS::FluidEOSRegistrator<PPFS::ClassName> reg_##ClassName; }


#endif /* PARALLELPLATEFLOWSOLVERBASE_H_ */
