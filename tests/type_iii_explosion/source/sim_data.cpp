#include "calc_init_cond.hpp"
#include "sim_data.hpp"

SimData::SimData(void):
  eos_(5./3.),
  rs_(eos_),
  bc_(rs_),
  sr_(),
  sim_
  (calc_init_cond(),
   eos_,
   rs_,
   bc_,
   sr_) {}

HydroSimulation& SimData::getSim(void)
{
  return sim_;
}
