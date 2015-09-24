#include "sim_data.hpp"
#include "my_main_loop.hpp"
#ifdef WITH_MPI
#include "boost/mpi.hpp"
#endif // WITH_MPI

int main(void)
{
#ifdef WITH_MPI
  boost::mpi::environment env;
  boost::mpi::communicator world;
#endif // WITH_MPI
  SimData sim_data
#ifdef WITH_MPI
     (world)
#endif // WITH_MPI
     ;
  HydroSimulation& sim = sim_data.getSim();
  my_main_loop
    (sim
#ifdef WITH_MPI
     ,world
#endif
     );
  return 0;
}
