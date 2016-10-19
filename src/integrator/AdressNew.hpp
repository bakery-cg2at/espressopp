/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
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
#ifndef _AdressNew_HPP
#define _AdressNew_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"
#include "FixedVSList.hpp"

#include "boost/signals2.hpp"


namespace espressopp {

namespace integrator {

class AdressNew: public Extension {

 public:
  real dhy;
  real pidhy2;
  real dex;
  real dex2;
  real dexdhy;
  real dexdhy2;

  AdressNew(shared_ptr<System> _system, shared_ptr<FixedVSList> _vslist);

  ~AdressNew();

  void setHyEx(Real3D center, real ex_, real hy_, bool is_sphere);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  virtual void connect();
  virtual void disconnect();

 private:
  boost::signals2::connection _runInit, _aftIntP;

  shared_ptr<FixedVSList> vs_list;
  Real3D adrCenter;
  bool isSphere;

  void updateWeights();

  virtual real weight(real);
  virtual real weightderivative(real);

  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}

}

#endif
