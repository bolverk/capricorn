#ifndef HYDRODYNAMICS_HPP
#define HYDRODYNAMICS_HPP 1

#include "hydro_snapshot.hpp"
#include "equation_of_state.hpp"
#include "boundary_condition.hpp"
#include "spatial_reconstruction.hpp"

double calc_max_time_step(const HydroSnapshot& hs,
			  const EquationOfState& eos);

vector<RiemannSolution> calc_fluxes(const HydroSnapshot& hs,
				    const RiemannSolver& rs,
				    const BoundaryCondition& bc,
				    const SpatialReconstruction& sr);

vector<double> calc_masses(const HydroSnapshot& hs);

vector<double> calc_momenta(const HydroSnapshot& hs);

vector<double> calc_energies(const HydroSnapshot& hs,
			     const EquationOfState& eos);

void update_momenta(const vector<RiemannSolution>& fluxes,
		    double dt,
		    vector<double>& momenta);

void update_energies(const vector<RiemannSolution>& fluxes,
		    double dt,
		    vector<double>& momenta);

void update_grid(const vector<RiemannSolution>& fluxes,
		 double dt,
		 vector<double>& grid);

vector<Primitive> retrieve_primitives(const vector<double>& grid,
				      const vector<double>& masses,
				      const vector<double>& momenta,
				      const vector<double>& energies,
				      const EquationOfState& eos);

#endif // HYDRODYNAMICS_HPP
