#include <cmath>
#include "ideal_gas.hpp"

IdealGas::IdealGas(double adiabatic_index):
  g_(adiabatic_index) {}

double IdealGas::dp2c(double density, double pressure) const
{
  return sqrt(g_*pressure/density);
}

double IdealGas::dp2e(double density, double pressure) const
{
  return (pressure/density)/(g_-1);
}

double IdealGas::de2p(double density, double energy) const
{
  return (g_-1)*energy*density;
}

double IdealGas::dp2s(double density, double pressure) const
{
  return pressure/pow(density,g_);
}

double IdealGas::ds2p(double density, double entropy) const
{
  return entropy*pow(density,g_);
}

double IdealGas::getAdiabaticIndex(void) const
{
  return g_;
}
