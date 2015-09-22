/*! \file hydro_snapshot.hpp
  \author Almog Yalinewich
  \brief Hydrodynamic snapshot
 */

#ifndef HYDRO_SNAPSHOT_HPP
#define HYDRO_SNAPSHOT_HPP 1

#include <vector>
#include "primitive.hpp"

using std::vector;

//! \brief Hydrodynamic snapshot
class HydroSnapshot
{
public:

  /*! \brief Class constructor
    \param grid_i Position of grid points (cell edges)
    \param cells_i Primitive variables
   */
  HydroSnapshot(const vector<double>& grid_i,
		const vector<Primitive>& cells_i);

  /*! \brief Copy constructor
    \param source Source
   */
  HydroSnapshot(const HydroSnapshot& source);

  //! \brief Position of grid points (cell edges)
  vector<double> grid;

  //! \brief Primitive variables
  vector<Primitive> cells;
};

#endif // HYDRO_SNAPSHOT_HPP
