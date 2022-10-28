//
// Created by root on 10/24/22.
//

#include "RelpermDriver.hpp"

#include "common/MpiWrapper.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

namespace geosx {

    using namespace dataRepository;
    using namespace constitutive;

    RelpermDriver::RelpermDriver(const geosx::string &name, geosx::dataRepository::Group *const parent) :
            TaskBase( name, parent )
    {
        enableLogLevelInput();

        registerWrapper( viewKeyStruct::relpermNameString(), &m_relpermName).
                setInputFlag(InputFlags::REQUIRED).
                setDescription("Relperm model to test");

        registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
                setInputFlag( InputFlags::REQUIRED ).
                setDescription( "Number of load steps to take" );

        registerWrapper( viewKeyStruct::outputString(), &m_outputFile ).
                setInputFlag( InputFlags::OPTIONAL ).
                setApplyDefaultValue( "none" ).
                setDescription( "Output file" );

        registerWrapper( viewKeyStruct::baselineString(), &m_baselineFile ).
                setInputFlag( InputFlags::OPTIONAL ).
                setApplyDefaultValue( "none" ).
                setDescription( "Baseline file" );
    }


    void RelpermDriver::outputResults()
    {
        // TODO: improve file path output to grab command line -o directory
        //       for the moment, we just use the specified m_outputFile directly

        FILE * fp = fopen( m_outputFile.c_str(), "w" );

        fprintf( fp, "# column 1 = time\n" );
        fprintf( fp, "# columns %d-%d = phase vol fractions\n", 2, 1+m_numPhases );
        fprintf( fp, "# columns %d-%d = phase relperm\n", 2+m_numPhases, 1+2*m_numPhases );

        for( integer n=0; n<m_table.size( 0 ); ++n )
        {
            for( integer col=0; col<m_table.size( 2 ); ++col )
            {
                fprintf( fp, "%.4e ", m_table( n, 0, col ) );
            }
            fprintf( fp, "\n" );
        }
        fclose( fp );


    }


    bool RelpermDriver::execute(const geosx::real64 GEOSX_UNUSED_PARAM( time_n ),
                                const geosx::real64 GEOSX_UNUSED_PARAM( dt ),
                                const geosx::integer GEOSX_UNUSED_PARAM( cycleNumber ),
                                const geosx::integer GEOSX_UNUSED_PARAM( eventCounter ),
                                const geosx::real64 GEOSX_UNUSED_PARAM( eventProgress ),
                                geosx::DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
                                {
                                    // this code only makes sense in serial

                                    GEOSX_THROW_IF( MpiWrapper::commRank() > 0, "PVTDriver should only be run in serial", std::runtime_error );


                                    constitutive::ConstitutiveManager & constitutiveManager = this->getGroupByPath< constitutive::ConstitutiveManager >( "/Problem/domain/Constitutive" );
                                    constitutive::RelativePermeabilityBase& baseRelperm = constitutiveManager.getGroup< constitutive::RelativePermeabilityBase >( m_relpermName );

                                    if( getLogLevel() > 0 )
                                    {
                                        GEOSX_LOG_RANK_0( "Launching PVT Driver" );
                                        GEOSX_LOG_RANK_0( "  Relperm .................. " << m_relpermName );
                                        GEOSX_LOG_RANK_0( "  Type ................... " << baseRelperm.getCatalogName() );
                                        GEOSX_LOG_RANK_0( "  No. of Phases .......... " << m_numPhases );
                                        GEOSX_LOG_RANK_0( "  Steps .................. " << m_numSteps );
                                        GEOSX_LOG_RANK_0( "  Output ................. " << m_outputFile );
                                        GEOSX_LOG_RANK_0( "  Baseline ............... " << m_baselineFile );
                                    }

                                    // create a dummy discretization with one quadrature point for
                                    // storing constitutive data

                                    conduit::Node node;
                                    dataRepository::Group rootGroup( "root", node );
                                    dataRepository::Group discretization( "discretization", &rootGroup );

                                    discretization.resize( 1 );   // one element
                                    baseRelperm.allocateConstitutiveData( discretization, 1 );   // one quadrature point

                                    constitutiveUpdatePassThru( baseRelperm, [&] ( auto & selectedRelpermModel )
                                    {
                                        using RELPERM_TYPE = TYPEOFREF( selectedRelpermModel );
                                        runTest< RELPERM_TYPE >( selectedRelpermModel, m_table );
                                    } );

                                    // move table back to host for output
                                    m_table.move( LvArray::MemorySpace::host );

                                    if( m_outputFile != "none" )
                                    {
                                        outputResults();
                                    }

                                    //TODO
//                                    if( m_baselineFile != "none" )
//                                    {
//                                        compareWithBaseline();
//                                    }

                                    return false;


                            }


                            void RelpermDriver::postProcessInput()
                            {


                                constitutive::ConstitutiveManager & constitutiveManager = this->getGroupByPath< constitutive::ConstitutiveManager >( "/Problem/domain/Constitutive" );
                                constitutive::RelativePermeabilityBase& baseRelperm = constitutiveManager.getGroup< constitutive::RelativePermeabilityBase >( m_relpermName );

                                m_numPhases = baseRelperm.numFluidPhases();


                                using PT = constitutive::RelativePermeabilityBase::PhaseType;
                                integer const ipWater = baseRelperm.getPhaseOrder()[PT::WATER];
                                integer const ipOil   = baseRelperm.getPhaseOrder()[PT::OIL];
                                integer const ipGas   = baseRelperm.getPhaseOrder()[PT::GAS];

                                real64 minSw, minSnw;
                                if( baseRelperm.numFluidPhases() > 2 ) {
                                    minSw = baseRelperm.getPhaseMinVolumeFraction()[ipWater];
                                    minSnw = baseRelperm.getPhaseMinVolumeFraction()[ipGas];
                                }
                                else
                                {
                                    if(ipWater<0)
                                    {
                                        minSw = 0;
                                        minSnw = baseRelperm.getPhaseMinVolumeFraction()[ipGas];
                                    }
                                    else if(ipGas<0)
                                    {
                                        minSnw = 0;
                                        minSw = baseRelperm.getPhaseMinVolumeFraction()[ipGas];
                                    }
                                }
                                real64 const dSw = (1-minSw-minSnw) / m_numSteps;
                                // set input columns


                                if( m_numPhases > 2) {
                                    m_table.resize( (m_numSteps+1)*(m_numSteps+1), 1, 1+2*m_numPhases );
                                    for (integer ni = 0; ni < m_numSteps + 1; ++ni) {
                                        for (integer nj = 0; nj < m_numSteps + 1; ++nj) {

                                            integer index = ni * (m_numSteps + 1) + nj;
                                            m_table(index, 0, TIME) = minSw + index * dSw;
                                            m_table(index, 0, ipWater + 1) = minSw + nj * dSw;
                                            m_table(index, 0, ipGas + 1) = minSnw + ni * dSw;
                                            m_table(index, 0, ipOil + 1) =
                                                    1. - m_table(index, 0, ipWater + 1) - m_table(index, 0, ipOil + 1);
                                        }
                                    }
                                }
                                else
                                {
                                    m_table.resize( m_numSteps+1, 1, 1+2*m_numPhases );
                                    for (integer ni = 0; ni < m_numSteps + 1; ++ni) {
                                        integer index = ni;
                                        m_table(index, 0, TIME) = minSw + index * dSw;
                                        if(ipWater<0) {
                                            m_table(index, 0, ipGas + 1) = minSnw + ni * dSw;
                                            m_table(index, 0, ipOil+1 ) = 1. - m_table(index, 0, ipGas+1);
                                        }
                                        else if(ipGas<0){
                                            m_table(index, 0, ipWater + 1) = minSw + ni * dSw;
                                            m_table(index, 0, ipOil+1 ) = 1. - m_table(index, 0, ipWater+1);
                                        }
                                    }

                                }


                            }

                            //TODO impl
                            //    void RelpermDriver::compareWithBaseline()

    REGISTER_CATALOG_ENTRY( TaskBase,
                            RelpermDriver,
                            string const &, dataRepository::Group * const )



};