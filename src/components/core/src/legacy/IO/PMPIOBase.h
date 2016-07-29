/*
 * PMPIOBase.h
 *
 *  Created on: Jul 17, 2012
 *      Author: settgast1
 */

#ifndef PMPIOBASE_H_
#define PMPIOBASE_H_

#include <mpi.h>
#include "pmpio.h"
#include <string>
#include <stdio.h>

template< typename FILETYPE >
class PMPIOBase
{
public:
  PMPIOBase();
  virtual ~PMPIOBase();

  /// Initializes the silo library
  virtual void Initialize( const PMPIO_iomode_t readwrite ) = 0;

  /// finishes up the silo library usage
  void Finish()
  {
    PMPIO_Finish(m_baton);
  }

  /// Wait for the Baton when doing PMPIO
  FILETYPE* WaitForBaton(const int domainNumber, const int cycleNum);

  /// Wait for the Baton when doing PMPIO
  FILETYPE* WaitForBaton( const std::string& filename );

  /// Hand off the Baton when doing PMPIO
  void HandOffBaton()
  {
    PMPIO_HandOffBaton(m_baton, m_dbFilePtr);
  }

  FILETYPE* GetFilePtr() {return m_dbFilePtr;}

protected:

  FILETYPE* m_dbFilePtr;

  /// total number of "parallel" files to write out
  int m_numGroups;

  /// the pmpio baton. A processor needs this to write to the file.
  PMPIO_baton_t *m_baton;

  /// root of the filename that we will be reading/writing
  std::string m_fileRoot;


};



template< typename FILETYPE >
PMPIOBase<FILETYPE>::PMPIOBase()
{
  // TODO Auto-generated constructor stub

}

template< typename FILETYPE >
PMPIOBase<FILETYPE>::~PMPIOBase()
{
  // TODO Auto-generated destructor stub
}



// *********************************************************************************************************************
/**
 *
 * @param domainNumber
 * @param cycleNum
 * @param fileName
 * @param nsName
 */
template< typename FILETYPE >
FILETYPE* PMPIOBase<FILETYPE>::WaitForBaton(const int domainNumber, const int cycleNum)
{

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int groupRank = PMPIO_GroupRank(m_baton, rank);
  char fileName[200] =
  { 0 };
  char dirName[200] =
  { 0 };

  if (groupRank == 0)
    sprintf(fileName, "%s_%06d", m_fileRoot.c_str(), cycleNum);
  else
    sprintf(fileName, "%s_%06d.%03d", m_fileRoot.c_str(), cycleNum, groupRank);

  sprintf(dirName, "domain_%04d", domainNumber);

  m_dbFilePtr = (FILETYPE *) PMPIO_WaitForBaton(m_baton, fileName, dirName);

  return m_dbFilePtr;
}

template< typename FILETYPE >
FILETYPE* PMPIOBase<FILETYPE>::WaitForBaton( const std::string& filename )
{

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int groupRank = PMPIO_GroupRank(m_baton, rank);
  char fileName[200] =
  { 0 };
  char dirName[200] =
  { 0 };

  if (groupRank == 0)
    sprintf(fileName, "%s", filename.c_str() );
  else
    sprintf(fileName, "%s.%03d", filename.c_str(), groupRank);

//  sprintf(dirName, "domain_%04d", domainNumber);

  m_dbFilePtr = (FILETYPE *) PMPIO_WaitForBaton(m_baton, fileName, dirName);

  return m_dbFilePtr;
}



#endif /* PMPIOBASE_H_ */
