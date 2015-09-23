#include "sim_data.hpp"
#include "my_main_loop.hpp"

int main(void)
{
  SimData sim_data;
  HydroSimulation& sim = sim_data.getSim();
  my_main_loop(sim);
  return 0;
}
