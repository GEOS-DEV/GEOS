#ifndef GET_AVAILABLE_MEMORY_HPP
#define GET_AVAILABLE_MEMORY_HPP

#include <sys/stat.h>
#if defined __MACH__
#include <sys/types.h>
#include <sys/sysctl.h>

#include <mach/host_info.h>
#include <mach/mach_host.h>
#include <mach/task_info.h>
#include <mach/task.h>
#endif

namespace geos
{
/**
 * @brief Retieves current available memory on host
 * @return the available memory in bytes.
 */
static inline size_t getAvailableMemory()
{
#if defined(__APPLE__) && defined(__MACH__)
  int mib[6];
  mib[0] = CTL_HW;
  mib[1] = HW_PAGESIZE;

  int pagesize;
  size_t length;
  length = sizeof( pagesize );
  if( sysctl( mib, 2, &pagesize, &length, NULL, 0 ) < 0 )
  {
    fprintf( stderr, "getting page size" );
  }

  mach_msg_type_number_t count = HOST_VM_INFO_COUNT;

  vm_statistics_data_t vmstat;
  if( host_statistics( mach_host_self(), HOST_VM_INFO, ( host_info_t ) &vmstat, &count ) != KERN_SUCCESS )
  {
    fprintf ( stderr, "Failed to get VM statistics." );
  }

  return vmstat.free_count * pagesize;
#else
  GEOS_ERROR_IF( percent > 100, "Error, percentage of memory should be smallerer than -100, check lifoOnHost (should be greater that -100)" );
  return (size_t)sysconf( _SC_AVPHYS_PAGES ) *(size_t) sysconf( _SC_PAGESIZE );
#endif
}
}

// remove these definitions from mach/boolean.h that can conflict with GEOS code (eg. InputFlags::FALSE)
#undef TRUE
#undef FALSE

#endif
