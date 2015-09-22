#include "vacuum.hpp"

Vacuum::Vacuum(const RiemannSolver& rs, double density):
  rs_(rs), density_(density) {}

RiemannSolution Vacuum::calcFluxLeft(const HydroSnapshot& hs) const
{
  return rs_(std::pair<Primitive,Primitive>(Primitive(density_,0,0),
					    hs.cells.front()));
}

RiemannSolution Vacuum::calcFluxRight(const HydroSnapshot& hs) const
{
  return rs_(std::pair<Primitive,Primitive>(hs.cells.back(),
					    Primitive(density_,0,0)));
}

