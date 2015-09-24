#ifndef CALC_INIT_COND_HPP
#define CALC_INIT_COND_HPP 1

#include "hydro_snapshot.hpp"
#ifdef WITH_MPI
#include <boost/mpi.hpp>
#endif // WITH_MPI

HydroSnapshot calc_init_cond
(
#ifdef WITH_MPI
 const boost::mpi::communicator& world
#endif // WITH_MPI
 );

#endif // CALC_INIT_COND
