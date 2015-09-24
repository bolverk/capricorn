#include <cmath>
#include "hdsim.hpp"
#include "hydrodynamics.hpp"

HydroSimulation::HydroSimulation
(const Geometry& geom,
 const HydroSnapshot& hs,
 const EquationOfState& eos,
 const RiemannSolver& rs,
 const BoundaryCondition& bc,
 const SpatialReconstruction& sr
#ifdef WITH_MPI
 ,const boost::mpi::communicator& world
#endif // WITH_MPI
 ):
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
  world_(world),
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
  vector<boost::mpi::request> left(2);
  vector<boost::mpi::request> right(2);
  if(world_.rank()>0){
    left.at(0) = world_.isend(world_.rank()-1,0,hs_.cells.front());
    left.at(1) = world_.irecv(world_.rank()-1,1,left_ghost_edge_);
  }
  if(world_.rank()<world_.size()-1){
    right.at(0) = world_.isend(world_.rank()+1,1,hs_.grid.back());
    right.at(1) = world_.irecv(world_.rank()+1,0,right_ghost_);
  }
  if(world_.rank()>0)
    boost::mpi::wait_all(left.begin(),left.end());
  if(world_.rank()<world_.size()-1)
    boost::mpi::wait_all(right.begin(),right.end());
}

const vector<RiemannSolution> HydroSimulation::calc_pcm_fluxes
(void) const
{
  const vector<Primitive>& cells = hs_.cells;
  vector<RiemannSolution> res(hs_.grid.size());
  if(world_.rank()==0){
    res.at(0) = bc_.calcFluxLeft(hs_);
    for(size_t i=1;i<res.size()-1;++i)
      res.at(i) = rs_
	(pair<Primitive,Primitive>
	 (cells.at(i-1),
	  cells.at(i)));
    res.back() = rs_
      (pair<Primitive,Primitive>
       (cells.back(),
	right_ghost_));
  }
  else{
    for(size_t i=0;i<res.size()-1;++i)
      res.at(i) = rs_
	(pair<Primitive,Primitive>
	 (cells.at(i),cells.at(i+1)));
    res.back() = world_.size() - 1 == world_.rank() ?
      bc_.calcFluxRight(hs_) :
      rs_(pair<Primitive,Primitive>
	  (cells.back(),right_ghost_));
  }
  return res;
}

namespace {
  template<class T> vector<T> prepend
  (const T& t,
   const vector<T>& v)
  {
    vector<T> res = v;
    res.insert(res.begin(),t);
    return res;
  }

  double calc_max_speed
  (const Primitive& cell,
   const EquationOfState& eos)
  {
    return fabs(cell.velocity)+
      eos.dp2c(cell.density,cell.pressure);
  }
}

double HydroSimulation::calc_time_step(void) const
{
  const vector<double> grid =
    world_.rank()==0 ?
    hs_.grid :
    prepend(left_ghost_edge_, hs_.grid);
  const vector<Primitive>& cells = hs_.cells;
  double max_its = 0;
  for(size_t i=0;i<cells.size();++i){
    const double max_speed =
      calc_max_speed(cells.at(i), eos_);
    const double width =
      grid.at(i+1) - grid.at(i);
    assert(width>0);
    assert(max_speed>=0);
    max_its = fmax
      (max_its,
       max_speed/width);
  }
  double res = 0;
  boost::mpi::all_reduce
    (world_,
     max_its,
     res,
     boost::mpi::maximum<double>());
  return 1.0/res;
}
#endif // WITH_MPI

void HydroSimulation::timeAdvance(double dt_candidate)
{

  cold_flows_.initialize_entropies(hs_,eos_,geom_);

#ifdef WITH_MPI
  const double dt = dt_candidate<=0 ?
    cfl_*calc_time_step() :
    fmin(cfl_*calc_time_step(),dt_candidate);
#else
  const double dt = dt_candidate<=0 ?
    calc_max_time_step(hs_,eos_) :
    fmin(dt_candidate, cfl_*calc_max_time_step(hs_,eos_));
#endif // WITH_MPI

  const vector<RiemannSolution> fluxes =
#ifdef WITH_MPI
    calc_pcm_fluxes();
#else
    calc_fluxes(hs_,rs_,bc_,sr_);
#endif // WITH_MPI

  update_momenta(fluxes,dt,hs_.grid,geom_,momenta_);

  update_energies(fluxes,dt,hs_.grid,geom_,energies_);

  update_grid(fluxes,dt,hs_.grid);

  hs_.cells = cold_flows_.retrieve_primitives
    (hs_.grid, geom_, masses_, momenta_, energies_, eos_);

#ifdef WITH_MPI
  huddle();
#endif // WITH_MPI

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
