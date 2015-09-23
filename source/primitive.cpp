#include "primitive.hpp"

Primitive::Primitive(void):
  density(0), pressure(0), velocity(0) {}

Primitive::Primitive(double density_i,
		     double pressure_i,
		     double velocity_i):
  density(density_i),
  pressure(pressure_i),
  velocity(velocity_i) {}

#ifdef WITH_MPI
template<class Archive>
void Primitive::serialize
(Archive& ar,
 const unsigned int /*version*/)
{
  ar & density;
  ar & pressure;
  ar & velocity;
}
#endif // WITH_MPI
