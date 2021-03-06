/*! \file hdsim.hpp
  \author Almog Yalinewich
  \brief Hydrodynamic simulation
 */

#ifndef HDSIM_HPP
#define HDSIM_HPP 1

#include <vector>
#include "hydro_snapshot.hpp"
#include "spatial_reconstruction.hpp"
#include "equation_of_state.hpp"
#include "riemann_solver.hpp"
#include "boundary_condition.hpp"
#include "cold_flows.hpp"

using std::vector;

//! \brief Hydrodynamic simulation
class HydroSimulation
{
public:

  /*! \brief Class constructor
    \param hs Hydrodynamic snapshot
    \param eos Equation of state
    \param rs Riemann solver
    \param bc Boundary condition
    \param sr Spatial reconstruction scheme
   */
  HydroSimulation(const HydroSnapshot& hs,
		  const EquationOfState& eos,
		  const RiemannSolver& rs,
		  const BoundaryCondition& bc,
		  const SpatialReconstruction& sr);

  /*! \brief Advances the simulation in time
    \param dt_candidate Suggestion for next time step
   */
  void timeAdvance(double dt_candidate = -1);

  /*! \brief Activates the cold flows
    \param threshold Ratio between thermal and total energy
   */
  void activateColdFlows(double threshold);

  /*! \brief Returns the hydrodynamic snapshot
    \return Hydrodynamic snapshot
   */
  const HydroSnapshot& getSnapshot(void) const;

  /*! \brief Returns the simulation time
    \return Simulation time
   */
  double getTime(void) const;

  /*! \brief Returns the number of times time advance was called
    \return Number of times time advance was called
   */
  int getCycle(void) const;

private:
  HydroSnapshot hs_;
  const EquationOfState& eos_;
  const RiemannSolver& rs_;
  const BoundaryCondition& bc_;
  const SpatialReconstruction& sr_;
  const vector<double> masses_;
  vector<double> momenta_;
  vector<double> energies_;
  double cfl_;
  double time_;
  int cycle_;
  ColdFlows cold_flows_;
};

#endif // HDSIM_HPP
