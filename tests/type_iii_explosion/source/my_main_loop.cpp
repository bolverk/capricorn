#include <cassert>
#include "hdsim.hpp"
#include "my_main_loop.hpp"
#include "hdf5_diagnostics.hpp"

void my_main_loop(HydroSimulation& sim)
{
  write_snapshot_to_hdf5(sim,"initial.h5");

  assert(sim.getSnapshot().cells.size()>10);
  const size_t index_ref = sim.getSnapshot().cells.size()-10;
  const double density_ref =
    sim.getSnapshot().cells.at(index_ref).density;
  while(sim.getSnapshot().cells.at(index_ref).density<
	1.01*density_ref){
    sim.timeAdvance();
    assert(sim.getCycle()<1e6);
  }

  write_snapshot_to_hdf5(sim,"final.h5");
}
