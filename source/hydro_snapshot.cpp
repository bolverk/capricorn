#include "hydro_snapshot.hpp"

HydroSnapshot::HydroSnapshot(const vector<double>& grid_i,
			     const vector<Primitive>& cells_i):
  grid(grid_i), cells(cells_i) {}

HydroSnapshot::HydroSnapshot(const HydroSnapshot& source):
  grid(source.grid), cells(source.cells) {}
