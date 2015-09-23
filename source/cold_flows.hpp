/*! \file cold_flows.hpp
  \author Almog Yalinewich
  \brief Deals with cold and fast flows
 */

#ifndef COLD_FLOWS_HPP
#define COLD_FLOWS_HPP 1

#include <vector>
#include "hydro_snapshot.hpp"
#include "equation_of_state.hpp"
#include "geometry.hpp"

using std::vector;

//! \brief Help dealing with cold flows
class ColdFlows
{
public:

  //! \brief Class constructor
  ColdFlows(void);

  /*! \brief Activates cold flows correction
    \param threshold Ratio between thermal and total energy below which correction will apply
   */
  void activate(double threshold);

  /*! \brief Initializes the entropies
    \param hs Hydrodynamic snapshot
    \param eos Equation of state
   */
  void initialize_entropies
  (const HydroSnapshot& hs,
   const EquationOfState& eos,
   const Geometry& geom);

  /*! \brief Retrieves primitives
    \param grid Computational grid
    \param masses Masses
    \param momenta Linear momenta
    \param energies Energies
    \param eos Equation of state
   */
  vector<Primitive> retrieve_primitives
  (const vector<double>& grid,
   const Geometry& geom,
   const vector<double>& masses,
   const vector<double>& momenta,
   vector<double>& energies,
   const EquationOfState& eos) const;

private:
  bool active_;
  double threshold_;
  vector<double> entropies_;
};

#endif // COLD_FLOWS_HPP
