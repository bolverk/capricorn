#include "free_flow.hpp"

FreeFlow::FreeFlow(void) {}

RiemannSolution FreeFlow::calcFluxLeft(const HydroSnapshot& hs) const
{
  return RiemannSolution(hs.cells.front().pressure,
			 hs.cells.front().velocity);
}

RiemannSolution FreeFlow::calcFluxRight(const HydroSnapshot& hs) const
{
  return RiemannSolution(hs.cells.back().pressure,
			 hs.cells.back().velocity);
}
