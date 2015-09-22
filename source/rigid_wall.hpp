#ifndef RIGID_WALL_HPP
#define RIGID_WALL_HPP 1

#include "boundary_condition.hpp"

//! \brief Rigid wall boundary conditions
class RigidWall: public BoundaryCondition
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
   */
  RigidWall(const RiemannSolver& rs);

  RiemannSolution calcFluxLeft(const HydroSnapshot& hs) const;

  RiemannSolution calcFluxRight(const HydroSnapshot& hs) const;

private:

  const RiemannSolver& rs_;
};

#endif // RIGID_WALL_HPP
