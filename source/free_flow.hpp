/*! \file free_flow.hpp
  \author Almog Yalinewich
  \brief Free flow boundary conditions
 */

#ifndef FREE_FLOW_HPP
#define FREE_FLOW_HPP 1

#include "boundary_condition.hpp"

//! \brief Free flow boundary conditions
class FreeFlow: public BoundaryCondition
{
public:

  //! \brief Class constructor
  FreeFlow(void);

  RiemannSolution calcFluxLeft(const HydroSnapshot& hs) const;

  RiemannSolution calcFluxRight(const HydroSnapshot& hs) const;

private:
};

#endif // FREE_FLOW_HPP
