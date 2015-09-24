#ifndef MY_MAIN_LOOP
#define MY_MAIN_LOOP 1

#ifdef WITH_MPI
#include "boost/mpi.hpp"
#endif // WITH_MPI

void my_main_loop
(HydroSimulation& sim
#ifdef WITH_MPI
 ,const boost::mpi::communicator& world
#endif // WITH_MPI
 );

#endif // MY_MAIN_LOOP
