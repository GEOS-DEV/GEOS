/*
 * GlobalIndexManager.h
 *
 *  Created on: Feb 21, 2012
 *      Author: settgast1
 */

#ifndef GLOBALINDEXMANAGER_H_
#define GLOBALINDEXMANAGER_H_

#include "typedefs.h"


#if 0
class GlobalIndexManager
{
public:

  static GlobalIndexManager& Instance()
  {
    static GlobalIndexManager theGlobalIndexManager;
    return theGlobalIndexManager;
  }

  int OwningRank( const globalIndex& gIndex )
  {
    return ( (gIndex & m_rankBits)>>m_numIndexBits );
  }

  localIndex OwningIndex( const globalIndex& gIndex )
  {
    return ( gIndex & m_indexBits );
  }


private:

  const int m_numIndexBits;
  const globalIndex m_indexBits;
  const globalIndex m_rankBits;

  GlobalIndexManager():
    m_numIndexBits(40),
    m_indexBits(0x000000ffffffffff),
    m_rankBits( 0xffffff0000000000)
  {}


  ~GlobalIndexManager() {}
  GlobalIndexManager( const GlobalIndexManager&);
  GlobalIndexManager& operator=( const GlobalIndexManager&);
};

#elif 0
class GlobalIndexManager
{
public:

  static GlobalIndexManager& Instance()
  {
    static GlobalIndexManager theGlobalIndexManager;
    return theGlobalIndexManager;
  }

  int OwningRank( const globalIndex& gIndex )
  {
    return ( gIndex/m_firstRank );
  }

  localIndex OwningIndex( const globalIndex& gIndex )
  {
    return ( gIndex%m_firstRank );
  }

  void SetIndex


private:

  const globalIndex m_firstRank;

  GlobalIndexManager():
    m_firstRank(1000000000000)
  {}


  ~GlobalIndexManager() {}
  GlobalIndexManager( const GlobalIndexManager&);
  GlobalIndexManager& operator=( const GlobalIndexManager&);
};
#else
class GlobalIndexManager
{
public:


  static int OwningRank( const globalIndex& gIndex )
  {
    return ( gIndex/firstRank() );
  }

  static localIndex OwningIndex( const globalIndex& gIndex )
  {
    return ( gIndex%firstRank() );
  }


  static globalIndex Index( const int rank, const localIndex lid )
  {
    return( firstRank()*rank + lid );
  }

  static localIndex LargestLocalIndex()
  {
    return firstRank()-1;
  }

private:

  static globalIndex firstRank()
  {
    return 1000000000000;
//  18446744073709551615
//    return 1000;
  }


  GlobalIndexManager()
  {}

  ~GlobalIndexManager() {}
  GlobalIndexManager( const GlobalIndexManager&);
  GlobalIndexManager& operator=( const GlobalIndexManager&);
};

#endif



#endif /* GLOBALINDEXMANAGER_H_ */
