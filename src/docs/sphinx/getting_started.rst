###############################################################################
Getting Started with GEOSX
###############################################################################

Getting Ready
=================================
GEOSX resides in a git repository hosted at https://github.com/GEOSX/GEOSX.
It is suggested that you setup ssh keys, and use ssh for your clones as discussed 
`here <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`__.
It is suggested that a directory be created to host the various clones that one may require for an effective development workflow, so the first step is to make a directory

1. Setup working directory

.. code-block:: sh

  mkdir geosx
  cd geosx

Downloading the Code
=================================
There are currently two separate repositories that should be downloaded.
The first is the main repository, which may be cloned and initialized by the following steps: 

2. Clone the main repository

.. code-block:: sh

   git clone git@github.com:GEOSX/GEOSX.git
   
   cd GEOSX
   
   git submodule init
   
   git submodule update
   
   cd ..


3. Clone the third-party libraries

.. code-block:: sh

   git clone git@github.com:GEOSX/thirdPartyLibs.git
   cd thirdPartyLibs
   git submodule init
   git submodule update
   cd ..


Compiling the Code
=================================

GEOSX compilations are typically driven by a hostconfig file, which reside in GEOSX/host-configs.
If your platform does not have a host-config in the repository, you are encouraged to maintain one.
If you are running on an LC system, there is already a hostconfig and copy of the thirdPartyLibs installed, and you can skip step 4.

The first step in compiling GEOSX is to run cmake and generate the makefiles.
Starting with the third-party libraries, the config-build.script will run cmake for you.
Note that the 'make' step should be run serially, as the indiviudal package builds are run in parallel by default.

4. Configure and make the third party libraries

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -hc ../GEOSX/host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make hdf5
   make -j1

The next step is to compile the main code. 
Again, the config-build sets up cmake for you.

5. Configure and make the main code

.. code-block:: sh

   cd ../../GEOSX
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make -j4

   
Running the Code
=================================

GEOSX executables read in a XML input file. A simple example XML is located
`here <https://github.com/GEOSX/GEOSX/blob/develop/src/components/core/tests/PhysicsSolvers/LaplaceFEM.xml/>`__. 
To execute a serial run enter the following command from a working directory:

.. code-block:: sh

    path-to-geosx-bin/geosx -i path-to-xml/LaplaceFEM.xml



    Caliper Timers
=================================

GEOSX is equipped with Caliper timers, `<https://github.com/LLNL/Caliper>`.
We integrate Caliper into GEOSX by marking source-code sections of interest such as compuational kernels or initialization steps.
Caliper is included in the GEOSX TPL library and is built by adding the following cmake configuration to a host-config file.

.. code-block:: sh

   option( ENABLE_CALIPER "Enables CALIPER" On )


In particular, the following macros may be used to annotate GEOSX.

* ``GEOSX_MARK_BEGIN(name)`` - Marks user defined code region. 

* ``GEOSX_MARK_END(name)`` - Marks end of user defined code region.

* ``GEOSX_CXX_MARK_LOOP_BEGIN(loop, loopName)`` - Marks the start of a loop.

* ``GEOSX_CXX_MARK_LOOP_ITERATION`` - Marks a loop iteration.

*  ``GEOSX_CXX_MARK_LOOP_END(loop)`` - Marks end of a loop.

*  ``GEOSX_CXX_MARK_FUNCTION`` - Marks a function.

As an example of annotation, we consider the following example:
   
.. code-block:: sh

  void scatter() {
    GEOSX_CXX_MARK_FUNCTION; // Mark the function. Exports "function"="scatter"
    // ...
  }

  int main(int argc, const char* argv[]){

    GEOSX_CXX_MARK_FUNCTION;

    GEOSX_MARK_BEGIN("setup");
    //intialization step
    GEOSX_MARK_END("setup");

    GEOSX_CXX_MARK_LOOP_BEGIN(elemLoop, "elemLoop");
    for(int elem = 0; i < noElements; ++i){

       GEOSX_CXX_MARK_LOOP_ITERATION(elemLoop,i);
       //computation...

       scatter();
    }
    GEOSX_CXX_MARK_LOOP_END(elemLoop, "elemLoop");
    
    return 0;
  }

Configuration for CALIPER is done by exporting environment variables, the simplest
way to get started is setting the following variable

* ``CALI_CONFIG_PROFILE=flat-function-profile``

We refer the reader to the full Caliper tutorial `<https://github.com/LLNL/Caliper>`, for further details.  

Post running the application a .cali file will be generated. Viewing the output is done using the cali-query
command: ``cali-query -t *.cali``. The location of cali-query will depend on where caliper is installed. To simplify usage,
we recommend defining an alias.

* ``alias cali-query=thirdPartyLibs/build-cori-intel-release/caliper/src/caliper-build/src/tools/cali-query/cali-query``

Below is the output generated by invoking cali-query. By default output is given in nano-seconds.
  
.. code-block:: sh

   time.duration time.timestamp time.offset function     annotation time.inclusive.duration loop     iteration#elemloop 
          550     1536782188         550 
          101     1536782188         651 main         
           36     1536782188         687 main         setup                           36 
           34     1536782188         721 main                                            
           48     1536782188         769 main                                            elemloop 
           47     1536782188         816 main                                            elemloop                  0 
           17     1536782188         833 main/scatter                                 17 elemloop                  0 
           12     1536782188         845 main                                         76 elemloop                  0 
           23     1536782188         868 main                                            elemloop                    
           33     1536782188         901 main                                            elemloop                  1 
            7     1536782188         908 main/scatter                                  7 elemloop                  1 
            7     1536782188         915 main                                         47 elemloop                  1 
            8     1536782188         923 main                                        202 elemloop                    
            9     1536782188         932 main                                        382                                     
            
