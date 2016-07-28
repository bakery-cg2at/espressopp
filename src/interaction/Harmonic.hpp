/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _INTERACTION_HARMONIC_HPP
#define _INTERACTION_HARMONIC_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "FixedPairListAdressInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the Harmonic potential.
    */
    class Harmonic : public PotentialTemplate< Harmonic > {
    private:
      real K;
      real r0;

    public:
      static void registerPython();

      Harmonic(): K(0.0), r0(0.0) {
        setShift(0.0);
        setCutoff(infinity);
      }

      Harmonic(real _K, real _r0, real _cutoff, real _shift) : K(_K), r0(_r0) {
        setShift(_shift);
        setCutoff(_cutoff);
      }

      Harmonic(real _K, real _r0,  real _cutoff) : K(_K), r0(_r0) {
        autoShift = false;
        setCutoff(_cutoff);
        setAutoShift();
      }

      // Setter and getter
      void setK(real _K) {
        K = _K;
        updateAutoShift();
      }
      
      real getK() const { return K; }

      void setR0(real _r0) { 
        r0 = _r0;
        updateAutoShift();
      }
      
      real getR0() const { return r0; }

      real _computeEnergySqrRaw(real distSqr) const {
        real energy = K * pow((sqrt(distSqr) - r0), 2);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        real r = sqrt(distSqr);
        real ffactor = -2.0 * K * (r - r0) / r;
        force = dist * ffactor;
        return true;
      }

      boost::python::list getParams() {
        python::list params;

        params.append(python::make_tuple("class", "Harmonic"));
        params.append(python::make_tuple("K", K));
        params.append(python::make_tuple("r0", r0));
        params.append(python::make_tuple("cutoff", getCutoff()));
        params.append(python::make_tuple("shift", getShift()));

        return params;
      }
    };

    // provide pickle support
    struct Harmonic_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(Harmonic const& pot)
      {
        real K;
        real r0;
        real rc;
        real sh;
        K = pot.getK();
        r0 = pot.getR0();
        rc = pot.getCutoff();
        sh = pot.getShift();
        return boost::python::make_tuple(K, r0, rc, sh);
      }
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

  }
}

#endif
