#include "simple_io.hpp"
#include <iostream>
#include <fstream>

using std::ofstream;
using std::endl;

void write_single_number_to_file(double number,
				 const string& fname)
{
  ofstream f(fname.c_str());
  f.precision(14);
  f << number << endl;
  f.close();
}
