#ifndef RIGID_WALL_HPP
#define RIGID_WALL_HPP 1

#include "boundary_condition.hpp"

class RigidWall: public BoundaryCondition
{
public:

  RigidWall(const RiemannSolver& rs);

  RiemannSolution calcFluxLeft(const HydroSnapshot& hs) const;

  RiemannSolution calcFluxRight(const HydroSnapshot& hs) const;

private:

  const RiemannSolver& rs_;
};

#endif // RIGID_WALL_HPP
