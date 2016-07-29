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
 * @file ParallelPlateFlowSolverBase.cpp
 * @author settgast1
 * @date Aug 3, 2011
 */

#include "Common/Common.h"
#include "PhysicsSolverStrings.h"

#include "ParallelPlateFlowSolverBase.h"
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"


using namespace PS_STR;

ParallelPlateFlowSolverBase::ParallelPlateFlowSolverBase( const std::string& name,
                                                          ProblemManagerT* const pm ):
SolverBase(name,pm),
m_ApertureTable(nullptr),
m_syncedFieldsB(),
m_tracerNames(),
m_tracerNamesMass(),
m_tracerNamesDMDT(),
m_tracerNamesVolumeFraction(),
m_tracerMassPtrs(),
m_tracerDMDTPtrs(),
m_volumeFractionPtrs(),
m_edgesToFaces(),
m_gravityVector(0.0)
{
}

ParallelPlateFlowSolverBase::~ParallelPlateFlowSolverBase()
{
  delete m_fluidEOS;
}


void ParallelPlateFlowSolverBase::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  SolverBase::ReadXML(hdn);

  // Permeability
  m_mu = hdn->GetAttributeOrDefault("mu","1.0e-3 N.s/m^2");
  m_SHP_FCT = hdn->GetAttributeOrDefault<realT>("shapefactor",1.0);
  m_min_aperture = hdn->GetAttributeOrDefault(MinimumApertureStr,"1 um");
  m_max_aperture = hdn->GetAttributeOrDefault(MaximumApertureStr,"4 mm");
  //m_min_hydro_aperture = hdn->GetAttributeOrDefault("minHydroAperture", 0.0);
  // This parameter seems to cause local stability issues.
  m_zeroApertureVolume = hdn->GetAttributeOrDefault("ZeroApertureVolume","0 mm");
  m_zeroApertureOffset = hdn->GetAttributeOrDefault("ZeroApertureOffset","0 mm");


  m_boundPhysicalAperture = hdn->GetAttributeOrDefault("boundPhysicalAperture",0);
  m_apertureTransition = hdn->GetAttributeOrDefault("apertureTransition",0.0);
  m_apertureSlopeReduction = hdn->GetAttributeOrDefault("apertureSlopeReduction", 1.0);

  m_permeabilityUpdateOption = hdn->GetAttributeOrDefault<int>("permeabilityUpdateOption",0);
  R1Tensor zeroVector;
  zeroVector *= 0.0;
  m_gravityVector = hdn->GetAttributeOrDefault<R1Tensor>("gravityVector", zeroVector);
  m_negativeFluidPressureSlopeReduction = hdn->GetAttributeOrDefault<realT>("negPressSlopeReduction", 1.0);
  m_cutoffPressure = hdn->GetAttributeOrDefault<realT>("cutoffPressure", -1.0e99);

  m_contactOffsetCutoff = hdn->GetAttributeOrDefault<realT>("contactOffsetCutoff", -1.0e99);
  // Fluid EOS
  /*
  m_bulk_modulus = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
  m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
  m_press_o = hdn->GetAttributeOrDefault("press_o","0.0 Pa");
  m_pressureCap = hdn->GetAttributeOrDefault("pressurecap","1.0e8 Pa");
  if (m_pressureCap > m_bulk_modulus)
  {
    throw GPException("The pressure cap is higher than the bulk modulus!");
  }
  */

  m_maxLeakOffRatio = hdn->GetAttributeOrDefault("maxLeakOffRatio","0.3");


#ifdef SRC_EXTERNAL

  TICPP::HierarchicalDataNode* leakoffNode = hdn->GetChild("LeakoffModel");

  if(leakoffNode)
  {
    m_leakoffModel.ReadXML(leakoffNode);
    m_leakoffCoef = m_leakoffModel.GetLeakoffCoefficient();

  } else {
    std::cout << "No leakoff model specified." <<std::endl;
    m_leakoffCoef = 0.0;
  }

  std::string tempString = hdn->GetAttributeString("leakoffCoefficient");
  tempString = hdn->GetAttributeString("CartersLeakoffCoefficient");
  if(!tempString.empty())
    throw GPException("Error! The old way of specifying leak off coefficient in the solver is not supported anymore.  Please use a leak off subblock.");

#endif

  std::string ApertureTableName = hdn->GetAttributeStringOrDefault("apertureTable","");
  if( !ApertureTableName.empty() )
  {
    Table<1,realT>* table = TableManager::Instance().GetTable<1,realT>(ApertureTableName);
    //m_zeroAperture = m_ApertureTable->Lookup(realT(0.0),TableInterpolation::linear);
    std::vector<realT> xvals = table->AxisValues(0);
    std::vector<realT> yvals = table->Values();

    if( xvals.back() != 0.0 )
    {
      throw GPException("Invalid aperture table. Last coord must be zero!!");
    }
    if( xvals.size() < 2 )
    {
      throw GPException("Invalid aperture table. Must have more than two points specified!!");
    }
    int n=xvals.size()-1;
    realT slope = (yvals[n]-yvals[n-1]) / (xvals[n]-xvals[n-1]) ;

    if( slope >= 1.0 )
    {
      throw GPException("Invalid aperture table. slope of last two points >= 1 is invalid!!");
    }


    this->m_apertureTransition = (yvals[n] - slope * xvals[n] ) / ( 1.0 - slope );

    for( int i=0 ; i<10 ; ++i )
    {
      xvals.push_back(m_apertureTransition*pow(10.0,i));
      yvals.push_back(m_apertureTransition*pow(10.0,i));
    }

    table->SetAxisValues(0,xvals);
    table->SetValues(yvals);

    const_cast<Table<1,realT>*& >(m_ApertureTable) = table;
//    std::cout<<xvals<<std::endl;
//    std::cout<<yvals<<std::endl;
  }


  // Multiphase tracers
  m_tracerNames = hdn->GetStringVector("inputFluidNames");
  m_multiphaseFlow = (m_tracerNames.size() >= 1);
 
  TICPP::HierarchicalDataNode* fluidEOShdn = hdn->GetChild("FluidEOS");
  if(fluidEOShdn)
  {
    std::string fluidEOSname = fluidEOShdn->GetAttributeString("name");
    m_fluidEOS = PPFS::newFluidEOS(fluidEOSname,fluidEOShdn);
  } 
  else 
  {
    // build pressure cap model
    m_bulk_modulus = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
    m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
    m_press_o = hdn->GetAttributeOrDefault("press_o","0.0 Pa");
    m_pressureCap = hdn->GetAttributeOrDefault("pressurecap","1.0e8 Pa");
    bool disablePressureCap = hdn->GetAttributeOrDefault<bool>("disablePressureCap", 0);

    m_maxViscosity = hdn->GetAttributeOrDefault<realT>("maxViscosity", 1.0e5);

    if (m_pressureCap > m_bulk_modulus)
    {
//      throw GPException("The pressure cap is higher than the bulk modulus!");
    }
    
    if (m_multiphaseFlow)
    {      
      m_multiphaseDensity = hdn->GetAttributeVector<realT>("multiphaseDensity", ",");
      m_multiphaseBulkModulus = hdn->GetAttributeVector<realT>("multiphaseBulkModulus", ",");
      m_multiphaseViscosity = hdn->GetAttributeVector<realT>("multiphaseViscosity", ",");
      m_multiphaseMixingMode = hdn->GetAttributeOrDefault<int>("multiphaseMixingMode", 0);
      m_multiphaseMixingThreshold_serialFlow = hdn->GetAttributeOrDefault<realT>("multiphaseThreshold_serialFlow", 0.95);
      m_multiphaseMixingThreshold_singlePhase = hdn->GetAttributeOrDefault<realT>("multiphaseThreshold_singlePhase", 0.95);
      m_multiphaseSetnames = hdn->GetStringVector("multiphaseSetnames");

      if (m_multiphaseSetnames.size() == 0)
      {
        std::cout << "Warning: No setnames were specified for the multiphase fluid solver" << std::endl;
      }

      std::cout << "Using MultiphasePressureEOS" << std::endl;
      m_fluidEOS = new PPFS::MultiphasePressureEOS(m_multiphaseDensity, m_multiphaseBulkModulus, m_multiphaseViscosity);
      // m_fluidEOS = new PPFS::MultiphasePressureCapEOS(m_multiphaseDensity, m_multiphaseBulkModulus, m_multiphaseViscosity, m_pressureCap);
    }
    else if (disablePressureCap)
    {
      std::cout << "Using PressureEOS" << std::endl;
      m_fluidEOS = new PPFS::PressureEOS(m_rho_o, m_bulk_modulus);
    }
    else
    {
      std::cout << "Using PressureCapEOS" << std::endl;
      m_fluidEOS = new PPFS::PressureCapEOS(m_rho_o, m_bulk_modulus, m_pressureCap);
    }
    

  }


  // Viscosity model

  m_usePowerlawFluid = false;
  TICPP::HierarchicalDataNode* powerLawFluidNode = hdn->GetChild("PowerlawFluid");
  if(powerLawFluidNode){
    //m_powerlawFluid.ReadXML(powerLawFluidNode);
    m_usePowerlawFluid = true;
    m_fluidModelPtr = new PowerlawFluidModel();
    m_fluidModelPtr->ReadXML(powerLawFluidNode);
  }

  TICPP::HierarchicalDataNode* hbFluidNode = hdn->GetChild("HerschelBulkleyFluid");
  if(hbFluidNode){
    m_usePowerlawFluid = true;
    m_fluidModelPtr = new HerschelBulkleyParallelPlateFluidModel();
    m_fluidModelPtr->ReadXML(hbFluidNode);

  }
  if(!m_usePowerlawFluid){
    m_fluidModelPtr = new PowerlawFluidModel();
  }

  m_updateFaceArea = hdn->GetAttributeOrDefault("updateFaceArea","0");


  // Timestep
  m_courant = hdn->GetAttributeOrDefault<realT>("courant",0.5);

  // Faceset
  m_flowFaceSetName = hdn->GetAttributeString("flowFaceSet");
  if(m_flowFaceSetName.empty()) m_flowFaceSetName = hdn->GetAttributeString("faceset");



  // Proppant solver
  m_proppantActive = FALSE;

  // Wellbore solver
  m_wellboreSolverName = hdn->GetAttributeString("wellboreSolverName");


  // Barton Joint Parameters
  std::string temp = hdn->GetAttributeString("BartonJointParameters"); // aperture at zero effective stress; reference stress; aperture at ref stress
    if( !temp.empty() )
    {
      R1Tensor tempArray;
      tempArray.StrVal( temp );
      m_wZeroStress = tempArray[0];
      realT stressRef = tempArray[1];
      realT wRef = tempArray[2];

      if (wRef < 0.99 * m_wZeroStress)
      {
        m_bBarton = (m_wZeroStress - wRef) / stressRef / wRef;
        m_aBarton = m_wZeroStress * (m_wZeroStress - wRef) / stressRef / wRef;
      }
      else
      {
        m_bBarton = 0;
        m_aBarton = 0;
        m_min_aperture = m_wZeroStress;
      }

      m_kContact = hdn->GetAttributeOrDefault<realT>("contactK",0.0);

    }
    else
    {
      m_bBarton = 0;
      m_aBarton = 0;
    }
}


void ParallelPlateFlowSolverBase::RegisterFields( FaceManagerT& faceManager, EdgeManagerT& edgeManager )
{

  faceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  faceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);
  faceManager.AddKeyedDataField<FieldInfo::pressure>();
  faceManager.AddKeyedDataField<FieldInfo::density>();
  faceManager.AddKeyedDataField<FieldInfo::mass>();
  faceManager.AddKeylessDataField<int>( "flowFaceType", true, true );
  faceManager.AddKeylessDataField<realT>( "faceArea", true, false );
  faceManager.AddKeylessDataField<realT>("flowRateBC", false, true);
  faceManager.AddKeylessDataField<realT>("referenceDensity", true, true);
  faceManager.AddKeylessDataField<realT>("faceMu",true,true);

  faceManager.AddKeylessDataField<realT>("contactOffset",true,true);

  edgeManager.AddKeylessDataField<int>( "flowEdgeType", true, true );
  edgeManager.AddMap< UnorderedVariableOneToManyRelation >( "edgeToFlowFaces");
  edgeManager.AddKeylessDataField<realT>(PermeabilityStr,true,true);
  edgeManager.AddKeylessDataField<realT>(PressureStr,true,true);
  edgeManager.AddKeylessDataField<realT>(VolumetricFluxStr,true,true);
  edgeManager.AddKeylessDataField<realT>("length",true,true);
  edgeManager.AddKeylessDataField<R1Tensor>("center",true,true);

  edgeManager.AddKeylessDataField<realT>(ViscosityStr,true,true);
  edgeManager.AddKeylessDataField<realT>("edgeDensity",true,true);

#ifdef SRC_EXTERNAL
  if (m_leakoffModel.GetLeakoffCoefficient() > 0)
  {
    faceManager.AddKeylessDataField<realT>("initialSaturatedTime", true, false);
    faceManager.AddKeylessDataField<realT>("totalLeakedThickness", true, true);
    faceManager.AddKeylessDataField<realT>("faceToMatrixLeakOffRate", true, true);
  }
#endif

  // Multiphase fields
  if (m_multiphaseFlow)
  {
    localIndex n_fluids = m_tracerNames.size();
    m_tracerNamesMass.resize(n_fluids);
    m_tracerNamesDMDT.resize(n_fluids);
    m_tracerNamesVolumeFraction.resize(n_fluids);

    m_tracerMassPtrs.resize(n_fluids);
    m_tracerDMDTPtrs.resize(n_fluids);
    m_volumeFractionPtrs.resize(n_fluids);

    for(localIndex ii=0; ii<n_fluids; ii++)
    {
      // Mass
      char sbuffer0 [100];
      int N = sprintf(sbuffer0, "fluidMass_%s", m_tracerNames[ii].c_str());
      std::string tmp_name = sbuffer0;
      std::string tmp_name_short = tmp_name.substr(0, N);
      faceManager.AddKeylessDataField<realT>(tmp_name_short, true, true);
      m_tracerNamesMass[ii] = tmp_name_short;

      // dM/dt
      char sbuffer1 [100];
      N = sprintf(sbuffer1, "fluidDMDT_%s", m_tracerNames[ii].c_str());
      tmp_name = sbuffer1;
      tmp_name_short = tmp_name.substr(0, N);
      faceManager.AddKeylessDataField<realT>(tmp_name_short, true, true);
      m_tracerNamesDMDT[ii] = tmp_name_short;

      // Volume Fraction
      char sbuffer2 [100];
      N = sprintf(sbuffer2, "fluidVolumeFraction_%s", m_tracerNames[ii].c_str());
      tmp_name = sbuffer2;
      tmp_name_short = tmp_name.substr(0, N);
      faceManager.AddKeylessDataField<realT>(tmp_name_short, true, true);
      m_tracerNamesVolumeFraction[ii] = tmp_name_short;
    }
  }


  // Define common visualization fields
  m_commonFields.push_back("Aperture");
  m_commonFields.push_back("Pressure");
  m_commonFields.push_back("Mass");
  m_commonFields.push_back("Volume");
  m_commonFields.push_back("FaceArea");
  m_commonFields.push_back("birthTime");
  m_commonFields.push_back("faceMu");
  m_commonFields.push_back("flowFaceType");
  m_commonFields.push_back("ghostRank");

  if (m_multiphaseFlow)
  {
    for(localIndex ii=0; ii<m_tracerNamesVolumeFraction.size(); ii++)
    {
      m_commonFields.push_back(m_tracerNamesVolumeFraction[ii]);
    }
  }


}


void ParallelPlateFlowSolverBase::InitializeCommunications( PartitionBase& partition )
{
  m_syncedFields.clear();
  m_syncedFieldsB.clear();

  // Multiphase fields
  if (m_multiphaseFlow)
  {
    for (localIndex ii=0; ii<m_tracerNames.size(); ii++)
    {
      m_syncedFieldsB[PhysicalDomainT::FiniteElementFaceManager].push_back(m_tracerNamesVolumeFraction[ii]);
    }
  }
    
  partition.SetBufferSizes(m_syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
  partition.SetBufferSizes(m_syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);
}


void ParallelPlateFlowSolverBase::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  // iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  // iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");

  // flowFaceType = -1;
  // flowEdgeType = -1;
  

  // Viscosity fields
  rArray1d& viscosity = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);
  rArray1d& edgeDensity = domain.m_feEdgeManager.GetFieldData<realT>("edgeDensity");
  viscosity = m_mu;
  edgeDensity = m_rho_o;

  rArray1d& refDensity = domain.m_feFaceManager.GetFieldData<realT>("referenceDensity");
  rArray1d& faceMu = domain.m_feFaceManager.GetFieldData<realT>("faceMu");
  refDensity = m_rho_o;
  faceMu = m_mu;  

  // Multiphase tracers
  if (m_multiphaseFlow)
  {
    const rArray1d& faceFluidVolume = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();


    for(localIndex ii=0; ii<m_tracerNames.size(); ii++)
    {
      
      rArray1d& fluidTracerVolume = domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesMass[ii]);
      rArray1d& fluidTracerDVDT = domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesDMDT[ii]);
      rArray1d& fluidTracerVolumeFraction = domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesVolumeFraction[ii]);

      fluidTracerDVDT = 0.0;
      if (ii == 0)
      {
        fluidTracerVolume = faceFluidVolume;
        fluidTracerVolumeFraction = 1.0;
      }  
      else
      {
        fluidTracerVolumeFraction = 0.0;
      }
    }
  }

#ifdef SRC_EXTERNAL
  if (m_leakoffModel.GetLeakoffCoefficient() > 0)
  {
    rArray1d& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
    initialSaturatedTime = std::numeric_limits<realT>::max();
  }
#endif
}

void ParallelPlateFlowSolverBase::AdvectMultiphaseFields(PhysicalDomainT& domain, SpatialPartition& partition, realT time, realT dt)
{
 
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const rArray1d& faceFluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
  const rArray1d& faceFluidMass  = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& faceFluidDensity  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  // const rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>( PS_STR::FaceAreaStr );
  const Array1dT<R1Tensor>& edgeCenters = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& edgeLengths = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr  );
  // const rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const Array1dT<R1Tensor>& faceVelocities = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );
  const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);
  const rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);
  const rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");
  std::vector<multiphaseSort> mpsort(m_tracerNames.size());
  
  // Apply volume fractions to given setnames (or to flow BC's if not defined)
  if (m_multiphaseSetnames.size() > 0)
  {
    BoundaryTracerValuesToSet( domain, time );
  }
  
  // Set initial values:
  // Strictly speaking, we should use old density (beginning-of-step density) here, but it might not be worthwhile to register a new field.
  for(localIndex ii=0; ii<faceFluidVolume.size(); ii++)
  {
    for (localIndex jj=0; jj<m_tracerNames.size(); jj++)
    {
      (*m_tracerMassPtrs[jj])[ii] = (*m_volumeFractionPtrs[jj])[ii]*faceFluidVolume_old[ii]*faceFluidDensity[ii];
    }
  }

  // Set dV/dT values to zero
  for(localIndex ii=0; ii<m_tracerNames.size(); ii++)
  {
    (*m_tracerDMDTPtrs[ii]) = 0.0;
  }


  // Estimate fluid flow between faces
  // This has to be done in exactly the same way as how it's done in Assemble.
  // A small consistency under normal circumstances might get magnified under certain conditions and may yield great error in results.
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  localIndex source, receiver;  
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    int numFaces = itr->second.size();
         
    for (int ia = 0; ia < numFaces-1; ++ia)
    {
      localIndex kf = itr->second[ia];

        for (int ib = ia+1; ib < numFaces; ++ib)
        {
          localIndex kfb = itr->second[ib];
          
          // I am calling Pxrho a "potential"
          realT deltaPotential = (faceFluidPressure[kf] + edgeDeltaP[eg] * 0.5) * faceFluidDensity[kf]
                                -(faceFluidPressure[kfb] - edgeDeltaP[eg] * 0.5) * faceFluidDensity[kfb] ;
          if ( deltaPotential > 0)
          {
            source = kf;
            receiver = kfb;
          }
          else
          {
            source = kfb;
            receiver = kf;
            deltaPotential *= -1.0;
          }

          // Estimate flow rate between faces
          realT kappa = 0.0;
          if(m_usePowerlawFluid)
          {
            kappa = TwoFacePermeability_PowerLaw(edgeCenters, edgeLengths, faceCenters, 
                                                 apertures, faceVelocities, facePackVfs, edgeMus[eg],
                                                 eg, kf, kfb);
          }
          else
          {
            kappa = TwoFacePermeability(edgeCenters, edgeLengths, faceCenters, 
                                        apertures, facePackVfs, edgeMus[eg],
                                        eg, kf, kfb);
          }
          
          realT estimatedMassFluxRate = kappa*deltaPotential;

          // Estimate change in fluid fraction mass
          if( kappa > 0.0 )
          {
            bool seriesFlow=FALSE;  

            if (m_multiphaseMixingMode == 1)
            {
              // Test for series flow in between the faces
              
              // Sort phases by mass
              for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
              {
                mpsort[jj].mass = (*m_tracerMassPtrs[jj])[source] + (*m_tracerMassPtrs[jj])[receiver];
                mpsort[jj].phase = jj;
              }
              std::sort(mpsort.begin(), mpsort.end(), by_mass());

              localIndex primary = mpsort[0].phase;
              bool isSerial = ((*m_volumeFractionPtrs[primary])[receiver] >= m_multiphaseMixingThreshold_serialFlow);
              bool isSinglePhase = (((*m_volumeFractionPtrs[primary])[source]*faceFluidVolume[source] + (*m_volumeFractionPtrs[primary])[receiver]*faceFluidVolume[receiver]) / (faceFluidVolume[source] + faceFluidVolume[receiver]) >= m_multiphaseMixingThreshold_singlePhase);

              if (isSerial && !isSinglePhase)
              {
                seriesFlow = TRUE;
              }
            }

            if (seriesFlow)
            {
              realT seriesDMDT = 0.0;
              for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
              {
                realT maxPhaseDMDT = (*m_tracerMassPtrs[jj])[source]/dt;
                realT dmdt = std::min(estimatedMassFluxRate-seriesDMDT, maxPhaseDMDT);

                (*m_tracerDMDTPtrs[mpsort[jj].phase])[source] -= dmdt;
                (*m_tracerDMDTPtrs[mpsort[jj].phase])[receiver] += dmdt;

                seriesDMDT += dmdt;
                if ((estimatedMassFluxRate - seriesDMDT) <= 0.0)
                {
                  break;
                }
              }
            }
            else
            {
              for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
              {
                realT qm = estimatedMassFluxRate * (*m_volumeFractionPtrs[jj])[source];
                (*m_tracerDMDTPtrs[jj])[source] -= qm;
                (*m_tracerDMDTPtrs[jj])[receiver] += qm;
              }
            }

          }
        }
      }
    }
    

  // Set volume fraction values
  for (localIndex ii=0; ii<faceFluidVolume.size(); ii++)
  {
    if ((flowFaceType[ii] == 1) && (faceFluidVolume[ii] > 0.0))
    {
      // Update fluid volumes
      realT totalFluidMass = 0.0;
      for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
      {
        realT newFluidMass = (*m_tracerMassPtrs[jj])[ii] + (*m_tracerDMDTPtrs[jj])[ii]*dt;
        // newFluidMass = std::max(newFluidMass, 0.0); //Commenting this out following Randy's suggestion as this line would add mass.  Having negative fluid mass will not break the flow solver.  It can be a signal of bigger problems.

        (*m_tracerMassPtrs[jj])[ii] = newFluidMass;
        totalFluidMass += newFluidMass;
      }

      // Update fluid volume fractions, correct volume to match actual value
      if (totalFluidMass > 0.0)
      {
        for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
        {
          (*m_volumeFractionPtrs[jj])[ii] = (*m_tracerMassPtrs[jj])[ii]/totalFluidMass;
          (*m_tracerMassPtrs[jj])[ii] = (*m_volumeFractionPtrs[jj])[ii]*faceFluidMass[ii];
        }
      }
    }
  }

  // Sync fields
  partition.SynchronizeFields(m_syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);           
}

void ParallelPlateFlowSolverBase::UpdateMultiphasePointers(PhysicalDomainT& domain)
  {
    for (localIndex ii=0; ii<m_tracerNames.size(); ii++)
    {
      m_volumeFractionPtrs[ii] = &(domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesVolumeFraction[ii]));
      m_tracerMassPtrs[ii] = &(domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesMass[ii]));
      m_tracerDMDTPtrs[ii] = &(domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesDMDT[ii]));
    }
  }   

  void ParallelPlateFlowSolverBase::GenerateMultiphaseArray(rArrayPtrs& multiphasePtrs, localIndex kf, rArray1d& multiphaseVolumeFraction)
  {
    realT scale = 0.0;
    for (localIndex ii=0; ii<m_tracerNames.size(); ii++)
    {
      multiphaseVolumeFraction[ii] = (*multiphasePtrs[ii])[kf];
      scale += multiphaseVolumeFraction[ii];
    }
    multiphaseVolumeFraction /= scale;
  }

  void ParallelPlateFlowSolverBase::GenerateMultiphaseBCArray(rArray1d& multiphaseVolumeFraction, realT time)
  {
    realT scale = 0.0;
    rArray1d t(1);
    t[0] = time;

    for(localIndex ii=0; ii<m_tracerNames.size(); ii++)
    {
      const realT val = TableManager::Instance().LookupTable<1>(m_tracerNames[ii], t);
      multiphaseVolumeFraction[ii] = val;
      scale += multiphaseVolumeFraction[ii];
    }
    multiphaseVolumeFraction /= scale;
  }

  // void ParallelPlateFlowSolverBase::BoundaryTracerValues( PhysicalDomainT& domain,ObjectDataStructureBaseT& object,
  //                                                         BoundaryConditionBase* const bc, const lSet& aset, const realT time, const realT dt )
  // {
  //   const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  //   const rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  //   // Read input tables
  //   rArray1d faceVolumeFraction(m_tracerNames.size());
  //   GenerateMultiphaseBCArray(faceVolumeFraction, time);

  //   // Set values at BC locations
  //   for(localIndex i =0; i < bc->m_setNames.size(); ++i)
  //   {
  //     std::map< std::string, lSet >::iterator setMap = domain.m_feFaceManager.m_Sets.find( bc->m_setNames[i] );
      
  //     if( setMap != domain.m_feFaceManager.m_Sets.end() )
  //     {
  //       lSet& set = setMap->second;
  //       lSet::const_iterator b=set.begin();

  //       for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a)
  //       {
  //         if (flowFaceType[*a]==1)
  //         {
  //           for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
  //           {
  //             (*m_volumeFractionPtrs[jj])[*a] = faceVolumeFraction[jj];
  //             (*m_tracerMassPtrs[jj])[*a] = faceVolumeFraction[jj]*faceFluidVolume[*a];
  //           }
  //         }
  //       }
  //     }
  //   }
  // }


  void ParallelPlateFlowSolverBase::BoundaryTracerValuesToSet( PhysicalDomainT& domain, const realT time )
  {
    const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
    const rArray1d& faceFluidMass  = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();    

    // Read input tables
    rArray1d faceVolumeFraction(m_tracerNames.size());
    GenerateMultiphaseBCArray(faceVolumeFraction, time);

    // Set values at BC locations
    for(localIndex ii=0; ii<m_multiphaseSetnames.size(); ii++)
    {
      std::map< std::string, lSet >::iterator setMap = domain.m_feFaceManager.m_Sets.find( m_multiphaseSetnames[ii] );
      
      if( setMap != domain.m_feFaceManager.m_Sets.end() )
      {
        lSet& set = setMap->second;
        lSet::const_iterator b=set.begin();

        for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a)
        {
          if (flowFaceType[*a]==1)
          {
            for(localIndex jj=0; jj<m_tracerNames.size(); jj++)
            {
              (*m_volumeFractionPtrs[jj])[*a] = faceVolumeFraction[jj];
              (*m_tracerMassPtrs[jj])[*a] = faceVolumeFraction[jj]*faceFluidMass[*a];
            }
          }
        }
      }
    }
  }

  void ParallelPlateFlowSolverBase::MarkFaceSaturationTime( PhysicalDomainT& domain, const realT time,const realT dt )
  {
    rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
    const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
    rArray1d& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");

      for( localIndex kf = 0; kf < domain.m_feFaceManager.DataLengths(); ++kf )
      {
        if (flowFaceType[kf] == 1)
        {
          if (pressure[kf] > 0.0 && initialSaturatedTime[kf] == std::numeric_limits<realT>::max())
          {
            initialSaturatedTime[kf] = std::max(time - dt/2, 0.0);
          }
        }

      }
  }



// realT ParallelPlateFlowSolverBase::CalculateMultiphaseViscosity(PhysicalDomainT& domain, lArray1d& faces)
// {
//   if (m_multiphaseFlow)
//   {
//     const rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
//     localIndex numFlowFaces = itr->second.size();
//     realT viscosity = 0.0;
//     realT volume = 0.0;

//     for(localIndex ii=0; ii<m_tracerNamesVolumeFraction.size(); ii++)
//     {
//       const rArray1d& fluidTracerVolumeFraction = domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesVolumeFraction[ii]);
      
//       for( localIndex jj=0; jj<numFlowFaces; jj++)
//       {
//         tmpVolume = fluidTracerVolumeFraction[jj]*faceFluidVolume[jj];
//         volume += tmpVolume;
//         viscosity += m_multiphaseViscosity[ii]*tmpVolume;
//       }
//     }
    
//     // Weight by total volume of fluid in connected elements
//     return viscosity/volume;    
//   }
//   else 
//   {      
//     return m_mu;
//   }
// }

// void ParallelPlateFlowSolverBase::CalculateMultiphaseDensityAndBulkModulus(PhysicalDomainT& domain, localIndex kf)
// {
//   if (m_multiphaseFlow)
//   {
//     realT density = 0.0;
//     realT bulkModulus = 0.0;
//     realT concentration = 0.0;

//     for(localIndex ii=0; ii<m_tracerNamesVolumeFraction.size(); ii++)
//     {
//       const rArray1d& fluidTracerVolumeFraction = domain.m_feFaceManager.GetFieldData<realT>(m_tracerNamesVolumeFraction[ii]);
        
//       if (m_multiphaseMixingMode == 0)
//       {
//         density += m_multiphaseDensity[ii]*fluidTracerVolumeFraction[kf];
//         bulkModulus += m_multiphaseBulkModulus[ii]*fluidTracerVolumeFraction[kf];
//       }
//       else
//       {
//         if (fluidTracerVolumeFraction[kf] > concentration)
//         {
//           concentration = fluidTracerVolumeFraction[kf];
//           density = m_multiphaseDensity[ii];
//           bulkModulus = m_multiphaseBulkModulus[ii];
//         }
//       }
//     }

//     this->m_rho_o = density;
//     this->m_bulk_modulus = bulkModulus;
     
//   }
// }


/**
 *  Power law fluid model
 */

void PowerlawFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  m_M = hdn->GetAttributeOrDefault("ConsistencyIndex","1.0e-3 N.s/m^2");
  m_n = hdn->GetAttributeOrDefault("FluidBehaviorIndex","1");
  m_phiM = PPFS::CalculateModifiedConsistencyIndex(m_M,m_n);
}

realT PowerlawFluidModel::CalculatePermeability(const realT la, const realT lb,
                                                const realT apa, const realT apb,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){
    return PPFS::CalculatePermeability_PowerLawFluid(la,lb,apa,apb, w, qMag, m_phiM, m_n, SHP_FCT);
}

// one sided permeability
realT PowerlawFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){
    return PPFS::CalculatePermeability_PowerLawFluid(l,ap, w, qMag, m_phiM, m_n, SHP_FCT);
}



/**
 *  Herschel Bulkley Parallel Plate Fluid Model
 */

void HerschelBulkleyParallelPlateFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  PowerlawFluidModel::ReadXML(hdn);
  m_tau_y = hdn->GetAttributeOrDefault("YieldStress","0.0");
  m_phi = m_phiM/m_M;

  /* build lookup table */
  rArray1d xs;
  rArray1d values;
  xs.resize(16);
  values.resize(16);
  for(int i =0 ; i < 16; i++)
  {
    realT x = -3 + 0.2*i;
    xs[i] = x;
    realT kappa = pow(10,x);
    values[i] = CalculateZp(kappa,m_n);
  }
  m_zp_kappaLookup.SetAxisValues(0,xs);
  m_zp_kappaLookup.SetValues(values);

  for(int i =0 ; i < 16; i++)
  {
    std::cout <<  "k: " <<  pow(10,xs[i]) << " v: " <<  values[i] << "\n";
  }
}

// analytical solution relating dimensionless parameter V to zp (dimensionless plug thickness)
realT HerschelBulkleyParallelPlateFluidModel::analyticalVFunc(realT zp, realT n){
  realT V = pow(1-zp,n+1)*pow(n/(n+1)*zp + 1,n);
  return V;
}
// derivative of analyticalVFunc wrt zp
realT HerschelBulkleyParallelPlateFluidModel::dVdz(realT zp, realT n){

    realT deriv = -(n+1)*pow(1-zp,n) *pow(n/(n+1)*zp + 1,n) +
            pow(1-zp,n+1)*(n*n/(n+1))*pow(n/(n+1)*zp + 1,n-1);
    return deriv;
}

// calculate Zp as a function of kappa = Kq^n
realT HerschelBulkleyParallelPlateFluidModel::CalculateZp(realT Kqn, realT n){

  realT zp = 0.5;
  realT dzp = 1;
  int numIters = 0;
  realT tol = 1e-8;
  // newton's method
  while(fabs(dzp) > tol*zp && numIters < 500){
    realT f = Kqn*zp - analyticalVFunc(zp,n);
  realT dfdz = Kqn - dVdz(zp,n);
  dzp = -f/dfdz;
  zp += dzp;
  numIters += 1;
  }
  if(numIters == 500)
    throw GPException("HerschelBulkleyParallelPlateFluidModel: Zp calculation failed to converge");
  return zp;
}

realT HerschelBulkleyParallelPlateFluidModel::CalculatePermeability(const realT la, const realT lb,
                                                const realT apa, const realT apb,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){

  realT zp(0.0);
  if(m_tau_y > 0.0){
    realT kappa = pow( (2*m_n+1)/(2*m_n),m_n)*pow( 4*qMag/w,m_n)* m_M/m_tau_y;
    if(kappa > 1000.0){
      zp = 1.0/kappa;
    } else if(kappa < 1e-3){
      realT denom = pow(1 + m_n/(m_n+1),m_n);
      zp = 1 - pow(kappa/denom , 1.0/(m_n+1));
    } else {
      realT xx[1];
      xx[0] = log10(kappa);
      zp = m_zp_kappaLookup.Lookup(xx);
    }
  }
  return PPFS::CalculatePermeability_HerschelBulkleyFluid(la,lb,apa,apb, w, qMag, m_phiM, m_n, zp, SHP_FCT);
}

// one sided permeability
realT HerschelBulkleyParallelPlateFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){

  realT zp(0.0);
  if(m_tau_y > 0.0){
    realT kappa = pow( (2*m_n+1)/(2*m_n),m_n)*pow( 4*qMag/w,m_n)* m_M/m_tau_y;
    if(kappa > 1000.0){
      zp = 1.0/kappa;
    } else if(kappa < 1e-3){
      realT denom = pow(1 + m_n/(m_n+1),m_n);
      zp = 1 - pow(kappa/denom , 1.0/(m_n+1));
    } else {
      realT xx[1];
      xx[0] = log10(kappa);
      zp = m_zp_kappaLookup.Lookup(xx);
    }
  }
  return PPFS::CalculatePermeability_HerschelBulkleyFluid(l,ap, w, qMag, m_phiM, m_n,zp, SHP_FCT);
}



//// Fluid EOS
///////////////

void PPFS::FluidEOSBase::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
}


//// Linear EOS

PPFS::LinearEOS::LinearEOS(TICPP::HierarchicalDataNode* hdn):
FluidEOSBase(hdn),
m_K_o(){
    ReadXML(hdn);
}

void PPFS::LinearEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  FluidEOSBase::ReadXML(hdn);
  m_K_o = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
}

REGISTER_FLUID_EOS( LinearEOS )



//// Pressure cap

PPFS::PressureCapEOS::PressureCapEOS(TICPP::HierarchicalDataNode* hdn):
    LinearEOS(hdn),
        m_pressureCap(){
    ReadXML(hdn);

}

void PPFS::PressureCapEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  LinearEOS::ReadXML(hdn);
  m_pressureCap = hdn->GetAttributeOrDefault(BulkModulusStr,"1.0e8 Pa");

  if (m_pressureCap > m_K_o)
  {
      throw GPException("PressureCapEOS: The pressure cap is higher than the bulk modulus!");
  }
}
REGISTER_FLUID_EOS( PressureCapEOS )


////  Multiphase pressure cap eos
PPFS::PressureEOS::PressureEOS(TICPP::HierarchicalDataNode* hdn):
    LinearEOS(hdn){}
void PPFS::PressureEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  LinearEOS::ReadXML(hdn);
}
REGISTER_FLUID_EOS( PressureEOS )


////  Multiphase pressure cap eos
PPFS::BiLinearEOS::BiLinearEOS(TICPP::HierarchicalDataNode* hdn):
    LinearEOS(hdn)
{}

void PPFS::BiLinearEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  LinearEOS::ReadXML(hdn);
  m_K_neg = hdn->GetAttributeOrDefault("negBulkModulus","2.0e9 Pa");

}
REGISTER_FLUID_EOS( BiLinearEOS )



////  Multiphase pressure cap eos
PPFS::MultiphasePressureEOS::MultiphasePressureEOS(TICPP::HierarchicalDataNode* hdn):
    LinearEOS(hdn){}
void PPFS::MultiphasePressureEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  LinearEOS::ReadXML(hdn);
}
REGISTER_FLUID_EOS( MultiphasePressureEOS )


//// Adiabatic eos

PPFS::AdiabaticEOS::AdiabaticEOS(TICPP::HierarchicalDataNode* hdn):
    FluidEOSBase(hdn),
    m_P_o(),
    m_gamma(),
    m_P_ref()
{
    ReadXML(hdn);
}

void PPFS::AdiabaticEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  FluidEOSBase::ReadXML(hdn);
  m_P_o = hdn->GetAttributeOrDefault("P_o","1e5 Pa");
  m_P_ref = hdn->GetAttributeOrDefault("P_ref","0");
  m_rho_o = hdn->GetAttributeOrDefault("rho_o","1.204 kg/m^3");  // reset default rho_o to reflect gas
  m_gamma = hdn->GetAttributeOrDefault("gamma","7.0/5.0"); // diatomic gas
}

REGISTER_FLUID_EOS( AdiabaticEOS )


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////
//
// Fluid EOS factory

typedef std::map<std::string, PPFS::FluidEOSInitializer*> FluidEOSCatalogueType;

FluidEOSCatalogueType & PPFS::getFluidEOSCatalogue(){
  static FluidEOSCatalogueType theCatalogue ;
  return theCatalogue;
}

void PPFS::getFluidEOSNames( std::vector<std::string>& nameList){

  using namespace PPFS;
  for(FluidEOSCatalogueType::const_iterator it = getFluidEOSCatalogue().begin();
      it != getFluidEOSCatalogue().end(); ++it){
        nameList.push_back(it->first);
  }
}

PPFS::FluidEOSBase* PPFS::newFluidEOS(const std::string& FluidEOSName , TICPP::HierarchicalDataNode* hdn)
{
  using namespace PPFS;

  FluidEOSInitializer* FluidEOSInitializer = getFluidEOSCatalogue()[FluidEOSName];
  FluidEOSBase *theNewFluidEOS = NULL;

  if(!FluidEOSInitializer)
      throw GPException("Could not create unrecognized FluidEOS: "+ FluidEOSName);

  theNewFluidEOS = FluidEOSInitializer->initializeFluidEOS( hdn );

  return theNewFluidEOS;
}

