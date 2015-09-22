/*! \file hll.hpp
  \author Almog Yalinewich
  \brief HLL riemann solver
 */

#ifndef HLL_HPP
#define HLL_HPP 1

#include "riemann_solver.hpp"
#include "equation_of_state.hpp"

//! \brief HLL riemann solver
class HLL: public RiemannSolver
{
public:

  /*! \brief Class constructor
    \param eos Equation of state
   */
  HLL(const EquationOfState& eos);

  RiemannSolution operator()
  (const pair<Primitive,Primitive>& data) const;

private:
  const EquationOfState& eos_;
};

#endif // HLL_HPP
