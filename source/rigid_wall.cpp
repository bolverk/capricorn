#include "rigid_wall.hpp"

namespace {
  Primitive reflect(const Primitive& p)
  {
    return Primitive
      (p.density,
       p.pressure,
       -p.velocity);		     
  }
}

RigidWall::RigidWall(const RiemannSolver& rs):
  rs_(rs) {}

RiemannSolution RigidWall::calcFluxLeft
(const HydroSnapshot& hs) const
{
  return rs_
    (pair<Primitive,Primitive>
     (reflect(hs.cells.front()),
      hs.cells.front()));
}

RiemannSolution RigidWall::calcFluxRight
(const HydroSnapshot& hs) const
{
  return rs_
    (pair<Primitive,Primitive>
     (hs.cells.back(),
      reflect(hs.cells.back())));
}
