*****************************
Basic profiling with CALIPER
*****************************

GEOSX is equipped with `Caliper <https://github.com/LLNL/Caliper>`_ timers.
We integrate Caliper into GEOSX by marking source-code sections of interest such as compuational kernels or initialization steps.
Caliper is included in the GEOSX TPL library and is built by adding the following cmake configuration to a host-config file.

.. code-block:: sh

   option( ENABLE_CALIPER "Enables CALIPER" On )



GEOSX/Caliper Annotation Macros
=====================================

The following macros may be used to annotate GEOSX:

* ``GEOSX_MARK_BEGIN(name)`` - Marks user defined code region. 

* ``GEOSX_MARK_END(name)`` - Marks end of user defined code region.

* ``GEOSX_MARK_LOOP_BEGIN(loop, loopName)`` - Marks the start of a loop. **Will not work with a lambda which captures by value**

* ``GEOSX_MARK_LOOP_ITERATION`` - Marks a loop iteration.

*  ``GEOSX_MARK_LOOP_END(loop)`` - Marks end of a loop.

*  ``GEOSX_MARK_FUNCTION`` - Marks a function.

As an example of annotation, we consider the following example:
   
.. code-block:: sh

  void scatter() {
    GEOSX_MARK_FUNCTION; // Mark the function. Exports "function"="scatter"
    // ...
  }

  int main(int argc, const char* argv[]){

    GEOSX_MARK_FUNCTION;

    GEOSX_MARK_BEGIN("setup");
    //intialization step
    GEOSX_MARK_END("setup");

    GEOSX_MARK_LOOP_BEGIN(myloop, elemLoop);
    for(int elem = 0; i < noElements; ++i){

       GEOSX_MARK_LOOP_ITERATION(myLoop,i);
       //computation...

       scatter();
    }
    GEOSX_MARK_LOOP_END(myLoop);
    
    return 0;
  }


Configuring  Caliper
=================================
  
Configuration for CALIPER is done by exporting environment variables, the simplest
way to get started is setting the following variable

* ``CALI_CONFIG_PROFILE=runtime-report``

We refer the reader to the full Caliper tutorial `Caliper <https://github.com/LLNL/Caliper>`_ , for further details.

Post running application, the output will be of the following form (where time is given in microseconds). 
  
.. code-block:: sh

   Path          sum#time.duration 
   main                 145.000000 
     elemLoop            79.000000 
       scatter         4210.000000 
    setup                28.000000 


To print min/max/avg time per MPI rank for annotated regions, one may export the following variables

* ``CALI_SERVICES_ENABLE=aggregate,event,mpi,mpireport,timestamp``
* ``CALI_TIMER_SNAPSHOT_DURATION=true``
* ``CALI_TIMER_INCLUSIVE_DURATION=false``
* ``CALI_MPIREPORT_CONFIG="SELECT statistics(sum#time.duration) GROUP BY annotation,function,loop FORMAT tree"``

