/*! \file spatial_reconstruction.hpp
  \author Almog Yalinewich
  \brief Abstract class for spatial reconstruction
 */

#ifndef SPATIAL_RECONSTRUCTION_HPP
#define SPATIAL_RECONSTRUCTION_HPP 1

#include <vector>
#include "hydro_snapshot.hpp"

using std::vector;
using std::pair;

//! \brief Abstract class for spatial reconstruction
class SpatialReconstruction
{
public:

  /*! \brief Interpolate values
    \param hs Hydrodynamic snapshot
    \return List of interpolated values
   */
  virtual vector<pair<Primitive, Primitive> > operator()
  (const HydroSnapshot& hs) const = 0;

  //! \brief Null constructor
  virtual ~SpatialReconstruction(void);
};

#endif // SPATIAL_RECONSTRUCTION_HPP
