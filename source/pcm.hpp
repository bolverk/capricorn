/*! \file pcm.hpp
  \author Almog Yalinewich
  \brief Piecewise constant interpolation
 */

#ifndef PCM_HPP
#define PCM_HPP 1

#include "spatial_reconstruction.hpp"

//! \brief Piecewise constant interpolation
class PCM: public SpatialReconstruction
{
public:

  PCM(void);

  vector<pair<Primitive,Primitive> > operator()
  (const HydroSnapshot& hs) const;

};

#endif // PCM_HPP
