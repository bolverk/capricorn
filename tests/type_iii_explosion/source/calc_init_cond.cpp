#include <cmath>
#include "calc_init_cond.hpp"
#ifdef WITH_MPI
#include <boost/mpi.hpp>
#endif // WITH_MPI

using std::size_t;

#ifdef WITH_MPI
namespace{
  vector<size_t> cumsum(const vector<size_t> v)
  {
    vector<size_t> res(1,0);
    for(size_t i=0;i<v.size();++i)
      res.push_back(v.at(i)+res.back());
    return res;
  }
}
#endif // WITH_MPI

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
#ifdef WITH_MPI
  boost::mpi::environment env;
  boost::mpi::communicator world;
  vector<size_t> cell_dist_table
    (static_cast<size_t>(world.size()));
  assert(cell_dist_table.size()>0);
  for(size_t i=0;cells.size();++i)
    ++cell_dist_table.at(i%cell_dist_table.size());
  const vector<size_t> boundaries =
    cumsum(cell_dist_table);
  const size_t left_boundary = boundaries.at(static_cast<size_t>(world.rank()));
  const size_t right_boundary = boundaries.at(static_cast<size_t>(world.rank()+1));
  vector<Primitive> local_cells;
  for(size_t i=left_boundary;i<right_boundary;++i)
    local_cells.push_back(cells.at(i));
  vector<double> local_grid;
  if(world.rank()==0)
    local_grid.push_back(cell_edges.at(0));
  for(size_t i=left_boundary;i<right_boundary;++i)
    local_grid.push_back(cell_edges.at(i+1));
  return HydroSnapshot(local_grid, local_cells);
#else
  return HydroSnapshot(cell_edges, cells);
#endif //WITH_MPI
}
