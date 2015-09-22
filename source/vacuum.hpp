/*! \file vacuum.hpp
  \author Almog Yalinewich
  \brief Vacuum boundary conditions
 */

#ifndef VACUUM_HPP
#define VACUUM_HPP 1

#include "boundary_condition.hpp"

class Vacuum: public BoundaryCondition
{
public:

  Vacuum(const RiemannSolver& rs, double density=0);

  RiemannSolution calcFluxLeft(const HydroSnapshot& hs) const;

  RiemannSolution calcFluxRight(const HydroSnapshot& hs) const;

private:
  const RiemannSolver& rs_;
  const double density_;
};

#endif // VACUUM_HPP
