#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include "hdsim.hpp"
#include "my_main_loop.hpp"
#include "hdf5_diagnostics.hpp"

using std::ofstream;
using std::endl;

namespace {
  void write_number
  (const double number,
   const string& fname)
  {
    ofstream f(fname.c_str());
    f << number << endl;
    f.close();
  }

  /*
  size_t find_positive_index
  (const vector<double>& edges,
   double r)
  {
    for(size_t i=0;i<edges.size();++i){
      if(edges.at(i)>r)
	return i;
    }
    throw "position was not found";
  }
  */

  /*
  class SelfSimilarSnapshots
  {
  public:

    SelfSimilarSnapshots(const double r_min):
      r_min_(r_min), counter_(0) {}

    void operator()(const HydroSimulation& sim) const
    {
      const HydroSnapshot hs = sim.getSnapshot();
      const double r = r_min_*pow(2.0,counter_);
      if(r>hs.grid.back())
	return;
      const size_t pindex = find_positive_index
	(hs.grid, r);
      if(!(hs.cells.at(pindex).velocity > 1e-2))
	return;
      write_snapshot_to_hdf5(sim,"snapshot_"+boost::lexical_cast<string>(counter_)+".h5");
      ++counter_;
    }

  private:
    const double r_min_;
    mutable int counter_;
    };*/

  bool shock_reached_end
  (const HydroSnapshot& hs
#ifdef WITH_MPI
   ,const boost::mpi::communicator& world
#endif // WITH_MPI
   )
  {
#ifdef WITH_MPI
    bool res;
    if(world.rank()==world.size()-1){
      assert(hs.cells.size()>10);
      res = hs.cells.at(hs.cells.size()-10).velocity>1e-2;
    }
    boost::mpi::broadcast(world,res,world.size()-1);
    return res;
#else
    return hs.cells.at(hs.cells.size()-1).velocity>1e-2;
#endif // WITH_MPI
  }
}

void my_main_loop
(HydroSimulation& sim
#ifdef WITH_MPI
 ,const boost::mpi::communicator& world
#endif // WITH_MPI
 )
{
#ifdef WITH_MPI
  write_snapshot_to_hdf5
    (sim,"initial_"+
     boost::lexical_cast<string>(world.rank())+
     ".h5");
#else
  write_snapshot_to_hdf5(sim,"initial.h5");
#endif // WITH_MPI

  assert(sim.getSnapshot().cells.size()>10);
  //  const SelfSimilarSnapshots diag(1);
  while(!shock_reached_end
	(sim.getSnapshot()
#ifdef WITH_MPI
	 ,world
#endif // WITH_MPI
	 )){
    sim.timeAdvance();
    assert(sim.getCycle()<1e6);
#ifdef WITH_MPI
    if(world.rank()==0)
#endif // WITH_MPI
      write_number(sim.getTime(),"time.txt");
    //    diag(sim);
  }

#ifdef WITH_MPI
  write_snapshot_to_hdf5
    (sim,"final_"+
     boost::lexical_cast<string>(world.rank())+
     ".h5");
#else
  write_snapshot_to_hdf5(sim,"final.h5");
#endif // WITH_MPI
}
