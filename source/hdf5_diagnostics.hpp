#ifndef HDF5_DIAGNOSTICS_HPP
#define HDF5_DIAGNOSTICS_HPP 1

#include <string>
#include "hdsim.hpp"

using std::string;

void write_snapshot_to_hdf5(const HydroSimulation& sim,
			    const string& fname);

#endif // HDF5_DIAGNOSTICS_HPP
