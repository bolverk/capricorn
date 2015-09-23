#include <cmath>
#include "calc_init_cond.hpp"

using std::size_t;

HydroSnapshot calc_init_cond(void)
{
  const double dr2r = 0.1;
  const double r_min = 1e-2;
  vector<double> cell_edges(200,0);
  cell_edges.at(0) = r_min;
  for(size_t i=1;i<cell_edges.size();++i)
    cell_edges.at(i) = cell_edges.at(i-1)*(1+dr2r);
  vector<Primitive> cells(cell_edges.size()-1);
  for(size_t i=0;i<cells.size();++i){
    const double r = 0.5*(cell_edges.at(i)+cell_edges.at(i+1));
    cells.at(i).density = pow(r,-3);
    cells.at(i).pressure = r<10*r_min ? 1e2 : 0;
    cells.at(i).velocity = 0;
  }
  return HydroSnapshot(cell_edges, cells);
}
