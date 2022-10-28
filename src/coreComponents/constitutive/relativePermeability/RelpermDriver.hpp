//
// Created by root on 10/24/22.
//

#ifndef GEOSX_RELPERMDRIVER_HPP_
#define GEOSX_RELPERMDRIVER_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geosx {


    class RelpermDriver : public TaskBase{

    public:
        RelpermDriver( const string & name, Group * const parent);

        static string catalogName() { return "RelpermDriver"; }

        void postProcessInput() override;

        virtual bool execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
                              integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                              integer const GEOSX_UNUSED_PARAM( eventCounter ),
                              real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                              DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override;

        /**
         * @brief Run test using loading protocol in table
         * @param i Fluid constitutive model
         * @param table Table with input/output time history
         */
        template< typename RELPERM_TYPE >
        void runTest(RELPERM_TYPE & relperm, arrayView3d< real64 > const & table );

        /**
         * @brief Ouput table to file for easy plotting
         */
        void outputResults();

        /**
         * @brief Read in a baseline table from file and compare with computed one (for unit testing purposes)
         */
        void compareWithBaseline();

    private:

        /**
         * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
         */
        struct viewKeyStruct
        {
            constexpr static char const * relpermNameString() { return "relperm"; }
            constexpr static char const * wettingFunctionString() { return "wettingControl"; }
            constexpr static char const * nonwettingFunctionString() { return "nonwettingControl"; }
            constexpr static char const * numStepsString() { return "steps"; }
            constexpr static char const * outputString() { return "output"; }
            constexpr static char const * baselineString() { return "baseline"; }
        };

        integer m_numSteps;      ///< Number of load steps
        integer m_numColumns;    ///< Number of columns in data table (depends on number of fluid phases)
        integer m_numPhases;     ///< Number of fluid phases
        integer m_numComponents; ///< Number of fluid components

        string m_relpermName;               ///< relPermType identifier
        string m_wettingFunctionName;    ///< Time-dependent function controlling pressure
        string m_nonwettingFunctionName; ///< Time-dependent function controlling temperature
        string m_outputFile;              ///< Output file (optional, no output if not specified)

        array3d< real64 > m_table; ///< Table storing time-history of input/output

        Path m_baselineFile; ///< Baseline file (optional, for unit testing of solid models)

        enum columnKeys { TIME, WET, NONWET=3  }; ///< Enumeration of "input" column keys for readability

        static constexpr real64 m_baselineTol = 1e-3; ///< Comparison tolerance for baseline results
    };


}

#endif //GEOSX_RELPERMDRIVER_HPP_
