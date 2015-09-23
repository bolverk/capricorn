#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "hdsim.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "rigid_wall.hpp"
#include "pcm.hpp"

class SimData
{
public:

  SimData(void);

  HydroSimulation& getSim(void);

private:
  const IdealGas eos_;
  const HLL rs_;
  const RigidWall bc_;
  const PCM sr_;
  HydroSimulation sim_;
};

#endif // SIM_DATA_HPP
