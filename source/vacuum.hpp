/*! \file vacuum.hpp
  \author Almog Yalinewich
  \brief Vacuum boundary conditions
 */

#ifndef VACUUM_HPP
#define VACUUM_HPP 1

#include "boundary_condition.hpp"

//! \brief Vacuum boundary conditions
class Vacuum: public BoundaryCondition
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
    \param density Tenuous density (in case vacuum is too extreme)
   */
  Vacuum(const RiemannSolver& rs, double density=0);

  RiemannSolution calcFluxLeft(const HydroSnapshot& hs) const;

  RiemannSolution calcFluxRight(const HydroSnapshot& hs) const;

private:
  const RiemannSolver& rs_;
  const double density_;
};

#endif // VACUUM_HPP
