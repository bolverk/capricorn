#include <cmath>
#include <limits>
#include <cassert>
#include "hydrodynamics.hpp"
#include "utilities.hpp"

double calc_max_time_step(const HydroSnapshot& hs,
			  const EquationOfState& eos)
{
  /*
  class MaxTimeStepCalculator: public Index2Member<pair<bool,double> >
  {
  public:

    MaxTimeStepCalculator(const HydroSnapshot& hs_i,
			  const EquationOfState& eos_i):
      hs_(hs_i), eos_(eos_i) {}

    size_t getLength(void) const
    {
      return hs_.cells.size();
    }

    pair<bool,double> operator()(size_t i) const
    {
      const double width = hs_.grid[i+1] - hs_.grid[i];
      const Primitive& cell = hs_.cells[i];
      const double max_speed = fabs(cell.velocity)+
	eos_.dp2c(cell.density, cell.pressure);
      return width/max_speed;
    }

  private:
    const HydroSnapshot& hs_;
    const EquationOfState& eos_;    
  } max_time_step_calculator(hs,eos);
  return min_term(serial_generate(max_time_step_calculator));*/

  vector<double> valid_time_steps;
  for(size_t i=0;i<hs.cells.size();++i){
    const Primitive cell = hs.cells[i];
    const double max_speed = fabs(cell.velocity)+
      eos.dp2c(cell.density, cell.pressure);
    if(std::numeric_limits<double>::epsilon()<max_speed){
      const double width = hs.grid[i+1] - hs.grid[i];
      assert(width>std::numeric_limits<double>::epsilon());
      valid_time_steps.push_back(width/max_speed);
    }
  }
  
  return min_term(valid_time_steps);
}

vector<RiemannSolution> calc_fluxes(const HydroSnapshot& hs,
				    const RiemannSolver& rs,
				    const BoundaryCondition& bc,
				    const SpatialReconstruction& sr)
{
  class FluxCalculator: public Index2Member<RiemannSolution>
  {
  public:

    FluxCalculator(const HydroSnapshot& hs_i,
		   const RiemannSolver& rs_i,
		   const BoundaryCondition& bc_i,
		   const SpatialReconstruction& sr_i):
      hs_(hs_i), rs_(rs_i), bc_(bc_i), 
      interpolated_(sr_i(hs_i)) {}

    size_t getLength(void) const
    {
      return hs_.grid.size();
    }

    RiemannSolution operator()(size_t i) const
    {
      if(i==0)
	return bc_.calcFluxLeft(hs_);
      else if(i==hs_.grid.size()-1)
	return bc_.calcFluxRight(hs_);
      else
	return rs_(interpolated_[i-1]);
    }

  private:
    const HydroSnapshot& hs_;
    const RiemannSolver& rs_;
    const BoundaryCondition& bc_;
    const vector<pair<Primitive, Primitive> > interpolated_;
  } flux_calculator(hs,rs,bc,sr);

  return serial_generate(flux_calculator);
}

vector<double> calc_masses(const HydroSnapshot& hs)
{
  class MassCalculator: public Index2Member<double>
  {
  public:

    MassCalculator(const HydroSnapshot& hs_i):
    hs_(hs_i) {}

    size_t getLength(void) const
    {
      return hs_.cells.size();
    }

    double operator()(size_t i) const
    {
      const double width = 
	hs_.grid[i+1] - hs_.grid[i];
      return width*hs_.cells[i].density;
    }

  private:
    const HydroSnapshot& hs_;
  } mass_calculator(hs);

  return serial_generate(mass_calculator);
}

vector<double> calc_momenta(const HydroSnapshot& hs)
{
  class MomentumCalculator: public Index2Member<double>
  {
  public:

    MomentumCalculator(const HydroSnapshot& hs_i):
    hs_(hs_i) {}

    size_t getLength(void) const
    {
      return hs_.cells.size();
    }

    double operator()(size_t i) const
    {
      const double width = 
	hs_.grid[i+1] - hs_.grid[i];
      return width*hs_.cells[i].density*hs_.cells[i].velocity;
    }

  private:
    const HydroSnapshot& hs_;
  } momentum_calculator(hs);

  return serial_generate(momentum_calculator);
}

vector<double> calc_energies(const HydroSnapshot& hs,
			     const EquationOfState& eos)
{
  class EnergyCalculator: public Index2Member<double>
  {
  public:

    EnergyCalculator(const HydroSnapshot& hs_i,
		     const EquationOfState& eos_i):
      hs_(hs_i), eos_(eos_i) {}

    size_t getLength(void) const
    {
      return hs_.cells.size();
    }

    double operator()(size_t i) const
    {
      const double width = 
	hs_.grid[i+1] - hs_.grid[i];
      const Primitive& cell = hs_.cells[i];
      const double kinetic_energy = 
	0.5*pow(cell.velocity,2);
      const double thermal_energy = 
	eos_.dp2e(cell.density,cell.pressure);
      return width*cell.density*
	(thermal_energy+kinetic_energy);
    }

  private:

    const HydroSnapshot& hs_;
    const EquationOfState& eos_;
  } energy_calculator(hs, eos);
  
  return serial_generate(energy_calculator);
}

void update_momenta(const vector<RiemannSolution>& fluxes,
		    double dt,
		    vector<double>& momenta)
{
  class MomentumDifferenceCalculator: public Index2Member<double>
  {
  public:

    MomentumDifferenceCalculator(const vector<RiemannSolution>& fluxes_i,
				 double dt_i):
      fluxes_(fluxes_i), dt_(dt_i) {}

    size_t getLength(void) const
    {
      return fluxes_.size() - 1;
    }

    double operator()(size_t i) const
    {
      return dt_*(fluxes_[i].pressure-
		  fluxes_[i+1].pressure);
    }

  private:

    const vector<RiemannSolution>& fluxes_;
    const double dt_;
  } momentum_difference_calculator(fluxes, dt);

  serial_increment(serial_generate(momentum_difference_calculator),
		   momenta);
}

void update_energies(const vector<RiemannSolution>& fluxes,
		     double dt,
		     vector<double>& energies)
{
  class EnergyDifferenceCalculator: public Index2Member<double>
  {
  public:

    EnergyDifferenceCalculator(const vector<RiemannSolution>& fluxes_i,
			       double dt_i):
      fluxes_(fluxes_i), dt_(dt_i) {}

    size_t getLength(void) const
    {
      return fluxes_.size()-1;
    }

    double operator()(size_t i) const
    {
      return dt_*(fluxes_[i].pressure*fluxes_[i].velocity-
		 fluxes_[i+1].pressure*fluxes_[i+1].velocity);
    }

  private:
    const vector<RiemannSolution>& fluxes_;
    const double dt_;
  } energy_difference_calculator(fluxes, dt);

  serial_increment(serial_generate(energy_difference_calculator),
		   energies);
}

void update_grid(const vector<RiemannSolution>& fluxes,
		 double dt,
		 vector<double>& grid)
{
  class DisplacementCalculator: public Index2Member<double>
  {
  public:

    DisplacementCalculator(const vector<RiemannSolution>& fluxes_i,
			   double dt_i):
      fluxes_(fluxes_i), dt_(dt_i) {}

    size_t getLength(void) const 
    {
      return fluxes_.size();
    }

    double operator()(size_t i) const
    {
      return fluxes_[i].velocity*dt_;
    }

  private:
    const vector<RiemannSolution>& fluxes_;
    const double dt_;
  } displacement_calculator(fluxes, dt);

  serial_increment(serial_generate(displacement_calculator),
		   grid);
}

vector<Primitive> retrieve_primitives(const vector<double>& grid,
				      const vector<double>& masses,
				      const vector<double>& momenta,
				      const vector<double>& energies,
				      const EquationOfState& eos)
{
  class PrimitiveCalculator: public Index2Member<Primitive>
  {
  public:

    PrimitiveCalculator(const vector<double>& grid_i,
			const vector<double>& masses_i,
			const vector<double>& momenta_i,
			const vector<double>& energies_i,
			const EquationOfState& eos_i):
      grid_(grid_i),
      masses_(masses_i),
      momenta_(momenta_i),
      energies_(energies_i),
      eos_(eos_i) {}

    size_t getLength(void) const
    {
      return masses_.size();
    }

    Primitive operator()(size_t i) const
    {
      const double width =
	grid_[i+1] - grid_[i];
      const double density = masses_[i]/width;
      const double velocity = momenta_[i]/masses_[i];
      const double thermal_energy = 
	energies_[i]/masses_[i] - 
	0.5*pow(velocity,2);
      const double pressure = 
	eos_.de2p(density, thermal_energy);
      return Primitive(density,
		       pressure,
		       velocity);
    }

  private:
    const vector<double>& grid_;
    const vector<double>& masses_;
    const vector<double>& momenta_;
    const vector<double>& energies_;
    const EquationOfState& eos_;
  } primitive_calculator(grid,
			 masses,
			 momenta,
			 energies,
			 eos);

  return serial_generate(primitive_calculator);
}
