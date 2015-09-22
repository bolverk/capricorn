#include "riemann_solver.hpp"

RiemannSolution::RiemannSolution(void):
  pressure(0), velocity(0) {}

RiemannSolution::RiemannSolution(double pressure_i,
				 double velocity_i):
  pressure(pressure_i), velocity(velocity_i) {}

RiemannSolver::~RiemannSolver(void) {}
