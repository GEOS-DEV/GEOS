*****************************
Caliper Profiling
*****************************

GEOSX is equipped with `Caliper timers <https://github.com/LLNL/Caliper>`_.
We integrate Caliper into GEOSX by marking source-code sections of interest such as computational kernels or initialization steps.
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
   
.. code-block:: c++

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


Configuring Caliper
=================================
  
Caliper configuration is done by specifying a string to initialize Caliper with via the
`-t` option. A few options are listed below but we refer the reader to
`Caliper <https://github.com/LLNL/Caliper/blob/releases/v2.3.0/doc/ConfigManagerAPI.md>`_ for the full Caliper tutorial.

* ``-t runtime-report`` Will make Caliper print aggregated timing information to standard out.
* ``-t runtime-report,aggregate_across_ranks=false`` Will make Caliper write per rank timing information to standard out.
    This isn't useful when using more than one rank but it does provide more information for single rank runs.
* ``-t spot()`` Will make Caliper output a `.cali` timing file that can be viewed in the Spot web server.


Using Adiak
=================================
Adiak is a library that allows the addition of meta-data to the Caliper Spot output, it is enabled with Caliper.
This meta-data allows you to easily slice and dice the timings available in the Spot web server. To export meta-data
use the `adiak::value` function.

See `Adiak <https://github.com/LLNL/Adiak/blob/f27ba674b88c2435e5e3245acbda9fc0a57bf88f/docs/Adiak%20API.docx>`_
for the full Adiak documentation.


Using Spot
=================================
To use Spot you will need an LC account and a directory full of `.cali` files you would like to analyse.
Point your browser to `Spot <https://lc.llnl.gov/spot2>`_ and open up the directory containing the timing files.
