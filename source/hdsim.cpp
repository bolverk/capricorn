#include <cmath>
#include "hdsim.hpp"
#include "hydrodynamics.hpp"
#include <boost/mpi.hpp>

HydroSimulation::HydroSimulation
(const Geometry& geom,
 const HydroSnapshot& hs,
 const EquationOfState& eos,
 const RiemannSolver& rs,
 const BoundaryCondition& bc,
 const SpatialReconstruction& sr):
  geom_(geom),
  hs_(hs),
  eos_(eos),
  rs_(rs),
  bc_(bc),
  sr_(sr),
  masses_(calc_masses(hs,geom)),
  momenta_(calc_momenta(hs,geom)),
  energies_(calc_energies(hs,eos,geom)),
  cfl_(0.3),
  time_(0),
  cycle_(0),
  cold_flows_()
#ifdef WITH_MPI
  ,
  left_ghost_edge_(0),
  right_ghost_()
#endif // WITH_MPI
{
#ifdef WITH_MPI
  huddle();
#endif // WITH_MPI  
}

#ifdef WITH_MPI
void HydroSimulation::huddle(void)
{
  boost::mpi::environment env;
  boost::mpi::communicator world;
  vector<boost::mpi::request> left(2);
  vector<boost::mpi::request> right(2);
  if(world.rank()>0){
    left.at(0) = world.isend(world.rank()-1,0,hs_.cells.front());
    left.at(1) = world.irecv(world.rank()-1,1,left_ghost_edge_);
  }
  if(world.rank()<world.size()-1){
    right.at(0) = world.isend(world.rank()+1,1,hs_.grid.back());
    right.at(1) = world.irecv(world.rank()+1,0,right_ghost_);
  }
  if(world.rank()>0)
    boost::mpi::wait_all(left.begin(),left.end());
  if(world.rank()<world.size()-1)
    boost::mpi::wait_all(right.begin(),right.end());
}
#endif // WITH_MPI

void HydroSimulation::timeAdvance(double dt_candidate)
{

  cold_flows_.initialize_entropies(hs_,eos_,geom_);

  const double dt = dt_candidate<=0 ?
    calc_max_time_step(hs_,eos_) :
    fmin(dt_candidate, cfl_*calc_max_time_step(hs_,eos_));

  const vector<RiemannSolution> fluxes = 
    calc_fluxes(hs_,rs_,bc_,sr_);

  update_momenta(fluxes,dt,hs_.grid,geom_,momenta_);

  update_energies(fluxes,dt,hs_.grid,geom_,energies_);

  update_grid(fluxes,dt,hs_.grid);

  hs_.cells = cold_flows_.retrieve_primitives
    (hs_.grid, geom_, masses_, momenta_, energies_, eos_);

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
