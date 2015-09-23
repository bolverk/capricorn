/*! \file primitive.hpp
  \author Almog Yalinewich
  \brief Primitive hydrodynamic variables
 */

#ifndef PRIMITIVE_HPP
#define PRIMITIVE_HPP 1

//! \brief Primitive hydrodynamic variables
class Primitive
{
public:

  //! \brief Null constructor (initialises everything to zero)
  Primitive(void);

  /*! \brief Class constructor
    \param density_i Density
    \param pressure_i Pressure
    \param velocity_i Velocity
   */
  Primitive(double density_i,
	    double pressure_i,
	    double velocity_i);

  //! \brief Density
  double density;

  //! \brief Pressure
  double pressure;

  //! \brief Velocity
  double velocity;

#ifdef WITH_MPI
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & density;
    ar & pressure;
    ar & velocity;
  }
#endif // WITH_MPI
};

#endif // PRIMITIVE_HPP
