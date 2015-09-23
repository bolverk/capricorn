#include <cmath>
#include "geometry.hpp"

Geometry::~Geometry(void) {}

Planar::Planar(void) {}

double Planar::calcVolume(double radius) const
{
  return radius;
}

double Planar::calcArea(double /*radius*/) const
{
  return 1;
}

Cylindrical::Cylindrical(void) {}

double Cylindrical::calcVolume(double radius) const
{
  return M_PI*pow(radius,2);
}

double Cylindrical::calcArea(double radius) const
{
  return 2.0*M_PI*radius;
}

Spherical::Spherical(void) {}

double Spherical::calcVolume(double radius) const
{
  return 4.0*M_PI*pow(radius,3.0)/3.0;
}

double Spherical::calcArea(double radius) const
{
  return 4.0*M_PI*pow(radius,2.0);
}
  
