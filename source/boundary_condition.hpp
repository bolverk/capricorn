/*! \file boundary_condition.hpp
  \author Almog Yalinewich
  \brief Boundary conditions
 */

#ifndef BOUNDARY_CONDITION_HPP
#define BOUNDARY_CONDITION_HPP 1

#include "riemann_solver.hpp"
#include "hydro_snapshot.hpp"

//! \brief Abstract class for boundary conditions
class BoundaryCondition
{
public:

  /*! \brief Calculate the fluxes at the left most edge
    \param hs Hydrodynamic snapshot
    \return Fluxes
   */
  virtual RiemannSolution calcFluxLeft(const HydroSnapshot& hs) const = 0;

  /*! \brief Calculates the fluxes at the right most edge
    \param hs Hydrodynamic snapshot
    \return Fluxes
   */
  virtual RiemannSolution calcFluxRight(const HydroSnapshot& hs) const = 0;

  //! \brief Class destructor
  virtual ~BoundaryCondition(void);
};

#endif // BOUNDARY_CONDITION_HPP
