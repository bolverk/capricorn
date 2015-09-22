#include <cmath>
#include "hdsim.hpp"
#include "hydrodynamics.hpp"

HydroSimulation::HydroSimulation(const HydroSnapshot& hs,
				 const EquationOfState& eos,
				 const RiemannSolver& rs,
				 const BoundaryCondition& bc,
				 const SpatialReconstruction& sr):
  hs_(hs), eos_(eos), rs_(rs), bc_(bc), sr_(sr),
  masses_(calc_masses(hs)),
  momenta_(calc_momenta(hs)),
  energies_(calc_energies(hs,eos)),
  cfl_(0.3), time_(0), cycle_(0), cold_flows_() {}

void HydroSimulation::timeAdvance(double dt_candidate)
{

  cold_flows_.initialize_entropies(hs_,eos_);

  const double dt = dt_candidate<=0 ?
    calc_max_time_step(hs_,eos_) :
    fmin(dt_candidate, cfl_*calc_max_time_step(hs_,eos_));

  const vector<RiemannSolution> fluxes = 
    calc_fluxes(hs_,rs_,bc_,sr_);

  update_momenta(fluxes,dt,momenta_);

  update_energies(fluxes,dt,energies_);

  update_grid(fluxes,dt,hs_.grid);

  /*
  hs_.cells = retrieve_primitives(hs_.grid,
				  masses_,
				  momenta_,
				  energies_,
				  eos_);
  */
  hs_.cells = cold_flows_.retrieve_primitives
    (hs_.grid, masses_, momenta_, energies_, eos_);

  time_ += dt;
  ++cycle_;
}

void HydroSimulation::activateColdFlows(double threshold)
{
  cold_flows_.activate(threshold);
}

const HydroSnapshot& HydroSimulation::getSnapshot(void) const
{
  return hs_;
}

double HydroSimulation::getTime(void) const
{
  return time_;
}

int HydroSimulation::getCycle(void) const
{
  return cycle_;
}
