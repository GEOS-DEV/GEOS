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
 * @file WriteFieldToFile.cpp
 * @author walsh24
 * @date June 19, 2011
 */

#include "SolverFactory.h"
#include "WriteFieldToFile.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"


#include "DataStructures/VectorFields/ElementRegionT.h"



//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

WriteFieldToFile::WriteFieldToFile(  const std::string& name,
                                     ProblemManagerT* const pm ):
SolverBase(name,pm),
m_isFirstTime(true)
{
}

WriteFieldToFile::~WriteFieldToFile()
{
  // TODO Auto-generated destructor stub
}

/*
 * <WriteFieldToFile name="wftf"      * name of the solver
 *            fileprefix="data"       * prefix for output files
 *            object="Face"           * location of field (Default=Node)
 *            fieldtype="Vector Scalar"          * type of field to write to file 
 *            fieldname="ReferencePosition Ca"   * name of field to write to file
 *            times="1s 10s"          * output times (alternatively use dt="1s" to set output frequency)
 *            setnames="Zmax"          * set to write to file - whole field written if not given
 *   />  
 */
void WriteFieldToFile::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{
  SolverBase::ReadXML( hdn );

  m_filePrefix = hdn->GetAttributeString("fileprefix");

  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object", PhysicalDomainT::FiniteElementNodeManagerStr() );
  objectTypeStr = hdn->GetAttributeStringOrDefault("objecttype",objectTypeStr);
  /*
  if(objectTypeStr.empty())
    throw GPException("Cannot write field to file without object type");
  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object", PhysicalDomainT::FiniteElementNodeManagerStr() );
  {
    std::string oldStr = hdn->GetAttributeStringOrDefault("objecttype", "" ); // throw error if using old syntax
    if(!oldStr.empty()) {
      throw GPException("UpdateParentIndicies: Attempting to set objecttype - use 'object' instead.");
    }
  }
  */

  m_objectType = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);

  m_fieldNames = hdn->GetStringVector("fieldname"); Trim(m_fieldNames);
  sArray1d fieldTypesStrs = hdn->GetStringVector("fieldtype");
  m_fieldTypes.resize(fieldTypesStrs.size());
  for( unsigned i = 0; i < fieldTypesStrs.size(); ++i)
      m_fieldTypes[i] = fromString<FieldType>(fieldTypesStrs[i]);
  if(m_fieldTypes.size() == 0) m_fieldTypes = std::vector<FieldType>(m_fieldNames.size(),FieldInfo::realField);
  
  if( m_fieldTypes.size() != m_fieldNames.size() ) 
      throw GPException("Error WriteFieldToFile: number of field names does not match number of field types.");
  
  m_regionName = hdn->GetAttributeStringOrDefault("regionname","");
  if( m_objectType==PhysicalDomainT::FiniteElementElementRegion && m_regionName.empty() ){
	  throw GPException("Error WriteFieldToFile: regionname must be supplied with Element object field types.");
  }

  // output frequency
  m_dt = hdn->GetAttributeOrDefault("dt",0.0);

  // get a list of output times
  std::vector<realT> times = hdn->GetAttributeVector<realT>("times");

  if(times.size() > 0)
  {

    m_outputTimes.assign( times.begin(),times.end() );
    m_outputTimes.push_back( std::numeric_limits<double>::max() );
    if(m_dt > 0.0){
      throw GPException("Error WriteFieldToFile: attempting to set output period (dt) and list of output times.");
    }
  } else {
    m_outputTimes.push_back( 0.0 );
    m_outputTimes.push_back( m_dt);
  }

  {
    sArray1d tempSetName;
    tempSetName = hdn->GetStringVector("setname");
    if (!tempSetName.empty())
      throw GPException("Error WriteFieldToFile: 'setname' is no longer supported.  Use 'setnames' instead.");
  }
  m_setNames = hdn->GetStringVector("setnames");
  // if(m_setNames.empty()) m_setNames = hdn->GetStringVector("setnames");

  // if true writes a single file, otherwise writes one file per output time
  m_appendToFile = hdn->GetAttributeOrDefault<bool>("append",false);
        
  // header string
  m_headerString = "#";
  for( unsigned i = 0; i < m_fieldNames.size(); ++i)
  {
    int fSize = FieldSize( m_fieldTypes[i] );
  	for(int j = 0; j < fSize; ++j){
      m_headerString += " " + m_fieldNames[i];
      if(fSize > 1) m_headerString += "["+toString<unsigned>(j)+"]";
  	}
  }
  m_headerString += "\n";
  
}

void WriteFieldToFile::RegisterFields( PhysicalDomainT& domain )
{
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectType,m_regionName);

  for(unsigned i =0; i < m_fieldNames.size();++i)
    objectManager.AddKeylessDataField( m_fieldTypes[i],  m_fieldNames[i], true, true );  
}



/**
 * 
 * 
**/
double WriteFieldToFile::TimeStep( const realT& time,
                                 const realT& dt ,
                                 const int cycleNumber,
                                 PhysicalDomainT& domain,
                                 const sArray1d& namesOfSolverRegions ,
                                 SpatialPartition& partition,
                                 FractunatorBase* const fractunator)
{  

  realT dt_return = dt;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max()*0.9;

  if(time >= m_outputTimes.front()*(1.0-1e-15) )
  {
    if(partition.m_rank == 0)
      std::cout << "Writing data field(s) to file." << std::endl;

    std::string filename = m_filePrefix;
    if(partition.m_size > 1) filename += "_" + toString(partition.m_rank);
    if(!m_appendToFile) filename += "_" + toString(time);
    filename += ".txt";

    if(m_appendToFile){
      std::ofstream fStream;//(filename.c_str(),std::ios::out | std::ios::app );
      if( m_isFirstTime ){
      	fStream.open(filename.c_str(),std::ios::out );
      	fStream << m_headerString;
      	m_isFirstTime =  false;
      } else {
      	fStream.open(filename.c_str(),std::ios::out | std::ios::app );
      }
      fStream << "# Time " << time << "\n";
      fStream.close();
    } else {
      // output header
      std::ofstream fStream(filename.c_str(),std::ios::out);
      fStream << m_headerString;
      fStream << "# Time " << time << "\n"; // slightly redundant as time is in filename
      fStream.close();
    }
  
    ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectType,m_regionName);

    if(m_objectType == PhysicalDomainT::FiniteElementElementRegion){
    	ElementRegionT *elementRegion = dynamic_cast<ElementRegionT *>(&objectManager);
    	elementRegion->UpdateElementFieldsWithGaussPointData();
    }

    if( m_setNames.empty() )
    {
  
      objectManager.WriteAsciiFieldData(m_fieldTypes,  m_fieldNames,filename,true);
      
    }else{
      for(unsigned i = 0; i < m_setNames.size(); ++i)
      {
        lSet& set = objectManager.GetSet(m_setNames[i]);
        objectManager.WriteAsciiFieldData(m_fieldTypes,  m_fieldNames,filename,set,true);
      }
    }

    if(m_outputTimes.back() < std::numeric_limits<double>::max() ) 
    {
        m_outputTimes.push_back(m_outputTimes.back() + m_dt);
    }

    m_outputTimes.pop_front();
  }

  return dt_return;
}


/**
 *
 * SetMaxStableTimeStep
 *
**/
void WriteFieldToFile::SetMaxStableTimeStep( const realT& time,
                                             PhysicalDomainT& ,
                                             const sArray1d& namesOfSolverRegions,
                                             SpatialPartition& partition __attribute__((unused)) )
{

  if(m_dt > 0.0)
  {
    m_stabledt.m_maxdt =  m_dt;
  }
}

REGISTER_SOLVER( WriteFieldToFile )
