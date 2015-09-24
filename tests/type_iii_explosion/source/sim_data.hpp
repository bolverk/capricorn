#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "hdsim.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "rigid_wall.hpp"
#include "pcm.hpp"
#ifdef WITH_MPI
#include <boost/mpi.hpp>
#endif // WITH_MPI

class SimData
{
public:

  SimData
  (
#ifdef WITH_MPI
   const boost::mpi::communicator& world
#endif // WITH_MPI
   );

  HydroSimulation& getSim(void);

private:
  const Spherical geom_;
  const IdealGas eos_;
  const HLL rs_;
  const RigidWall bc_;
  const PCM sr_;
  HydroSimulation sim_;
};

#endif // SIM_DATA_HPP
