#include "cold_flows.hpp"
#include "utilities.hpp"

ColdFlows::ColdFlows(void):
  active_(false),
  threshold_(0),
  entropies_() {}

void ColdFlows::activate(double threshold)
{
  active_ = true;
  threshold_ = threshold;
}

void ColdFlows::initialize_entropies(const HydroSnapshot& hs,
				     const EquationOfState& eos)
{
  if(active_){
    class EntropyCalculator: public Index2Member<double>
    {
    public:
      EntropyCalculator(const HydroSnapshot& hs_i,
			const EquationOfState& eos_i):
	hs_(hs_i), eos_(eos_i) {}
      
      size_t getLength(void) const
      {
	return hs_.cells.size();
      }

      double operator()(size_t i) const
      {
	const double width = hs_.grid[i+1] - hs_.grid[i];
	return width*hs_.cells[i].density*
	  eos_.dp2s(hs_.cells[i].density,
		    hs_.cells[i].pressure);
      }
      
    private:
      const HydroSnapshot& hs_;
      const EquationOfState& eos_;
    } entropy_calculator(hs, eos);

    entropies_ = serial_generate(entropy_calculator);
  }
}

vector<Primitive> ColdFlows::retrieve_primitives
(const vector<double>& grid,
 const vector<double>& masses,
 const vector<double>& momenta,
 vector<double>& energies,
 const EquationOfState& eos) const
{
  class PrimitiveCalculator: public Index2Member<Primitive>
  {
  public:

    PrimitiveCalculator(const vector<double>& grid_i,
			const vector<double>& masses_i,
			const vector<double>& momenta_i,
			vector<double>& energies_i,
			const EquationOfState& eos_i,
			const vector<double>& entropies_i,
			const bool active_i,
			const double threshold_i):
      grid_(grid_i),
      masses_(masses_i),
      momenta_(momenta_i),
      energies_(energies_i),
      eos_(eos_i),
      entropies__(entropies_i),
      active__(active_i),
      threshold__(threshold_i){}

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
      const double kinetic_energy = 
	0.5*pow(velocity,2);
      const double total_energy = energies_[i]/masses_[i];
      const double thermal_energy =
	total_energy - kinetic_energy;
      const bool cold_flow_flag = 
	active__ &&
	(threshold__>thermal_energy/total_energy ||
	 total_energy<0);
      const double pressure = cold_flow_flag ?
	eos_.ds2p(density, entropies__[i]/masses_[i]) :
	eos_.de2p(density, thermal_energy);
      if(cold_flow_flag)
	energies_[i] = masses_[i]*eos_.dp2e(density,
					    pressure);
      return Primitive(density, pressure, velocity);
    }

  private:
    const vector<double>& grid_;
    const vector<double>& masses_;
    const vector<double>& momenta_;
    vector<double>& energies_;
    const EquationOfState& eos_;
    const vector<double>& entropies__;
    const bool active__;
    const double threshold__;
  } primitive_calculator(grid,masses,momenta,energies,eos,
			 entropies_, active_, threshold_);

  return serial_generate(primitive_calculator);
}
