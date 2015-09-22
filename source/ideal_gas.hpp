/*! \file ideal_gas.hpp
  \author Almog Yalinewich
  \brief Ideal gas equation of state
 */

#ifndef IDEAL_GAS_HPP
#define IDEAL_GAS_HPP 1

#include "equation_of_state.hpp"

//! \brief Ideal gas equation of state
class IdealGas: public EquationOfState
{
public:

  /*! \brief Class constructor
    \param adiabatic_index Adiabatic index
   */
  IdealGas(double adiabatic_index);

  double dp2c(double density, double pressure) const;

  double dp2e(double density, double pressure) const;

  double de2p(double density, double energy) const;

  double dp2s(double density, double pressure) const;

  double ds2p(double density, double entropy) const;

  /*! \brief Returns the adiabatic index
    \return Adiabatic index
   */
  double getAdiabaticIndex(void) const;

private:

  double g_;
};

#endif // IDEAL_GAS_HPP
