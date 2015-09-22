#include <H5Cpp.h>
#include "hdf5_diagnostics.hpp"
#include "utilities.hpp"

using H5::PredType;
using H5::DataSet;
using H5::FloatType;
using H5::DataSpace;
using H5::H5File;

namespace {
  void write_std_vector_to_hdf5
  (H5File& file,
   vector<double> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = num_list.size();
    DataSpace dataspace(1, dimsf);
    
    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace);

    /* To do: Replace raw pointers with smart pointers
     */
    double *data = new double[static_cast<int>(num_list.size())];
    for(size_t i=0;i<num_list.size();++i)
      data[i] = num_list[i];
    dataset.write(data, PredType::NATIVE_DOUBLE);
    delete[] data;
  }
}

void write_snapshot_to_hdf5
(const HydroSimulation& sim,
 const string& fname)
{
  H5File file(H5std_string(fname),
	      H5F_ACC_TRUNC);

  write_std_vector_to_hdf5(file,
			   vector<double>(1,sim.getTime()),
			   "time");
  write_std_vector_to_hdf5(file,
			   sim.getSnapshot().grid,
			   "grid");
  class PropertyExtractor: public Index2Member<double>
  {
  public:
    
    PropertyExtractor(const vector<Primitive>& cells,
		      const double Primitive::*mem_ptr):
      cells_(cells), mem_ptr_(mem_ptr) {}

    size_t getLength(void) const
    {
      return cells_.size();
    }

    double operator()(size_t i) const
    {
      return cells_[i].*mem_ptr_;
    }

  private:
    const vector<Primitive>& cells_;
    const double Primitive::*mem_ptr_;
  };
 
  write_std_vector_to_hdf5(file,
			   serial_generate
			   (PropertyExtractor(sim.getSnapshot().cells,
					      &Primitive::density)),
			   "density");
  write_std_vector_to_hdf5(file,
			   serial_generate
			   (PropertyExtractor(sim.getSnapshot().cells,
					      &Primitive::pressure)),
			   "pressure");
  write_std_vector_to_hdf5(file,
			   serial_generate
			   (PropertyExtractor(sim.getSnapshot().cells,
					      &Primitive::velocity)),
			   "velocity");
}
