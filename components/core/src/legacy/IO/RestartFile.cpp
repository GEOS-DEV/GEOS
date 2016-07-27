/*
 * RestartFile.cpp
 *
 *  Created on: Jul 18, 2012
 *      Author: settgast1
 */

#include "RestartFile.h"


static void* PMPIO_ReadRestartCreate(const char *fname, const char *nsname, void *userData);
static void* PMPIO_WriteRestartCreate(const char *fname, const char *nsname, void *userData);

static void* PMPIO_ReadRestartOpen(const char *fname, const char *nsname, PMPIO_iomode_t ioMode, void *userData);
static void* PMPIO_WriteRestartOpen(const char *fname, const char *nsname, PMPIO_iomode_t ioMode, void *userData);

static void PMPIO_RestartClose(void *file, void *userData);


RestartFile::RestartFile()
{
  // TODO Auto-generated constructor stub

}

RestartFile::~RestartFile()
{
  // TODO Auto-generated destructor stub
}

void RestartFile::Initialize( const PMPIO_iomode_t readwrite )
{
  // Ensure all procs agree on numGroups, driver and file_ext

  MPI_Bcast(&m_numGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Initialize PMPIO, pass a pointer to the driver type as the user data.

  if( readwrite==PMPIO_READ )
  {
    m_baton = PMPIO_Init(m_numGroups, readwrite, MPI_COMM_WORLD, 1,
                         PMPIO_ReadRestartCreate, PMPIO_ReadRestartOpen, PMPIO_RestartClose,NULL);
  }
  else if( readwrite==PMPIO_WRITE )
  {
    m_baton = PMPIO_Init(m_numGroups, readwrite, MPI_COMM_WORLD, 1,
                         PMPIO_WriteRestartCreate, PMPIO_WriteRestartOpen, PMPIO_RestartClose,NULL);
  }

}


static void* PMPIO_ReadRestartCreate(const char *fname, const char *nsname, void *userData)
{
  /*
    iBinStream restartFile;
    restartFile.open( fname );

    return (void *) &restartFile;
    */
  throw GPException("RestartFile::PMPIO_ReadRestartCreate. this shouldn't be called");
}

static void* PMPIO_WriteRestartCreate(const char *fname, const char *nsname, void *userData)
{
    oBinStream* const restartFile = new oBinStream;
    restartFile->open( fname, true );

    return (void *) restartFile;
}




static void* PMPIO_ReadRestartOpen(const char *fname, const char *nsname, PMPIO_iomode_t ioMode, void *userData)
{
  iBinStream* const restartFile = new iBinStream;
  restartFile->open( fname );

  return (void *) restartFile;
}

static void* PMPIO_WriteRestartOpen(const char *fname, const char *nsname, PMPIO_iomode_t ioMode, void *userData)
{
  oBinStream* const restartFile = new oBinStream;
  restartFile->open( fname );

  return (void *) restartFile;
}

/*-----------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     Callbacks
 * Purpose:     Impliment the close callback
 *-----------------------------------------------------------------------------
 */
static void PMPIO_RestartClose(void *file, void *userData)
{
  BinStream* fileStream = static_cast<BinStream*>(file);
  fileStream->close();
  if( file!=NULL )
  {
    delete fileStream;
  }
}
