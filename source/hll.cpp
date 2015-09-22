#include <cmath>
#include <limits>
#include "hll.hpp"

HLL::HLL(const EquationOfState& eos):
  eos_(eos) {}

namespace {
  pair<double, double> calc_wave_speeds(const pair<Primitive,Primitive>& data,
					const EquationOfState& eos)
  {
    const double vl = data.first.velocity;
    const double cl = eos.dp2c(data.first.density,
			       data.first.pressure);
    const double vr = data.second.velocity;
    const double cr = eos.dp2c(data.second.density,
			       data.second.pressure);
    return pair<double, double>(fmin(vl-cl,vr-cr),
				fmax(vl+cl,vr+cr));
  }

  RiemannSolution calc_interface_values
  (const pair<double,double>& wave_speeds,
   const pair<Primitive,Primitive>& data)
  {
    const double dl = data.first.density;
    const double pl = data.first.pressure;
    const double vl = data.first.velocity;
    const double dr = data.second.density;
    const double pr = data.second.pressure;
    const double vr = data.second.velocity;
    const double sl = wave_speeds.first;
    const double sr = wave_speeds.second;
    const double v = (-pl+pr+dl*vl*(sl-vl)+dr*vr*(vr-sr))/
      (dl*(sl-vl)+dr*(vr-sr));
    const double p = (dl*(sl-vl)*(pr+dr*(sr-vr)*(vl-vr))+dr*pl*(vr-sr))/
      (dl*(sl-vl)+dr*(vr-sr));
    return RiemannSolution(p,v);
  }
}

RiemannSolution HLL::operator()(const pair<Primitive,Primitive>& data) const
{
  if(std::numeric_limits<double>::epsilon()>data.first.pressure &&
     std::numeric_limits<double>::epsilon()>fabs(data.first.velocity) &&
     std::numeric_limits<double>::epsilon()>data.second.pressure &&
     std::numeric_limits<double>::epsilon()>fabs(data.second.velocity))
    return RiemannSolution(0,0);
  return calc_interface_values(calc_wave_speeds(data,eos_),data);
}
