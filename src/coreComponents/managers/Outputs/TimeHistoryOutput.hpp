#include "OutputBase.hpp"

namespace geosx
{
  class TimeHistoryOutput : public OutputBase
  {
    TimeHistoryOutput( std::string const & name,
                       Group * const parent );
    virtual ~TimeHistoryOutput() override;

    static string CatalogName() { return "TimeHistoryOutput"; }
    virtual void SetupDirectoryStructure();

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override;

    /// Write one final output as the code exits
    virtual void Cleanup( real64 const time_n,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override;

    struct viewKeysStruct
    {
      static constexpr auto xxxString = "xxxString";
    } timeHistoryOutputViewKeys;
    private:
      //
  };
}


void TimeHistoryOutput::Execute( real64 const time_n,
                                 real64 const dt,
                                 integer const cycleNumber,
                                 integer const eventCounter,
                                 real64 const eventProgress,
                                 dataRepository::Group * domain )
{
  HDFFile hdf_file("time_history");
  array2d<real64> field; // retrieve from a manager depending on the input xml
  array1d<localIndex> local_indices; // retreive from probe
  array2d<localIndex> components(num_cmps); // parse from input xml

  HDFTimeHistoryTableIO veloc_hist<decltype(field)>("Velocity","vel",0,local_indices,components,"nd");
  veloc_hist.OpenTable( hdf_file );
  veloc_hist.AppendRow( time_n, veloc_field );
  veloc_hist.CloseTable( );
}