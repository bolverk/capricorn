/*! \file equation_of_state.hpp
  \author Almog Yalinewich
  \brief Abstract class for equation of state
 */

#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP 1

//! \brief Abstract class for equation of state
class EquationOfState
{
public:

  /*! \brief Calculates the speed of sound
    \param density Density
    \param pressure Pressure
    \return Speed of sound
   */
  virtual double dp2c(double density, double pressure) const = 0;

  /*! \brief Calculates specific energy (energy per mass)
    \param density Density
    \param pressure Pressure
    \return Specific energy (energy per mass)
   */
  virtual double dp2e(double density, double pressure) const = 0;

  /*! \brief Calculates the pressure
    \param density Density
    \param energy Specific energy (energy per mass)
    \return Pressure
   */
  virtual double de2p(double density, double energy) const = 0;

  /*! \brief Calculates the entropy
    \param density Density
    \param pressure Pressure
    \return Entropy
   */
  virtual double dp2s(double density, double pressure) const = 0;

  /*! \brief Calculates the pressure
    \param density Density
    \param entropy Entropy
    \return Pressure
   */
  virtual double ds2p(double density, double entropy) const = 0;

  //! \brief Class destructor
  virtual ~EquationOfState(void);
};

#endif // EQUATION_OF_STATE_HPP
