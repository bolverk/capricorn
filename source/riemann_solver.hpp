/*! \file riemann_solver.hpp
  \author Almog Yalinewich
  \brief Riemann solver
 */

#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP 1

#include <utility>
#include "primitive.hpp"

using std::pair;

//! \brief Container for values of the Riemann solution
class RiemannSolution
{
public:

  //! \brief Null constructor (sets everything to zero)
  RiemannSolution(void);

  /*! \brief Class constructor
    \param pressure_i Pressure
    \param velocity_i Velocity
   */
  RiemannSolution(double pressure_i,
		  double velocity_i);

  //! \brief Pressure
  double pressure;

  //! \brief Velocity
  double velocity;
};

//! \brief Abstract class for Riemann solver
class RiemannSolver
{
public:

  /*! \brief Solves the Riemann problem
    \param data Hydrodynamic conditions on both sides of the interface
    \return Riemann solution
   */
  virtual RiemannSolution operator()
  (const pair<Primitive,Primitive>& data) const = 0;

  //! \brief Class destructor
  virtual ~RiemannSolver(void);
};

#endif
