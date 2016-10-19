/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  
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

#include "python.hpp"
#include "mpi.hpp"
#include "AdressNew.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "FixedVSList.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"
#include <iomanip>
#include <sstream>

namespace espressopp {

namespace integrator {

using namespace espressopp::iterator;

LOG4ESPP_LOGGER(AdressNew::theLogger, "AdressNew");

AdressNew::AdressNew(shared_ptr<System> _system, shared_ptr<FixedVSList> _vslist)
    : Extension(_system), vs_list(_vslist) {
  LOG4ESPP_INFO(theLogger, "construct AdressNew");
  type = Extension::Adress;
}

void AdressNew::setHyEx(Real3D center, real ex_, real hy_, bool is_sphere) {
  dex = ex_;
  dhy = hy_;
  pidhy2 = M_PI / (dhy * 2.0);
  dex2 = dex * dex;
  dexdhy = dex + hy_;
  dexdhy2 = dexdhy * dexdhy;

  adrCenter = center;
  isSphere = is_sphere;
  if (isSphere)
    throw std::runtime_error("is sphere not supported yet!");
}


AdressNew::~AdressNew() {
  LOG4ESPP_INFO(theLogger, "~AdressNew");
  disconnect();
}

void AdressNew::disconnect() {
  _runInit.disconnect();
  _aftIntP.disconnect();
}

void AdressNew::connect() {
  // connection to after runInit()
  _runInit = integrator->runInit.connect(
      boost::bind(&AdressNew::updateWeights, this), boost::signals2::at_front);
  _aftIntP = integrator->aftIntP.connect(
      boost::bind(&AdressNew::updateWeights, this), boost::signals2::at_front);
}

void AdressNew::updateWeights() {
  System &system = getSystemRef();

  // Update weights of the particles.
  FixedVSList::GlobalTuples vs = vs_list->globalTuples;
  FixedVSList::GlobalTuples::iterator it = vs.begin();

  Real3D d1(0.0, 0.0, 0.0);
  real d1sq;
  real w;

  for (; it != vs.end(); ++it) {
    Particle *vp = system.storage->lookupRealParticle(it->first);
    if (vp) {
      system.bc->getMinimumImageVector(d1, vp->position(), adrCenter);
      d1sq = d1[0]*d1[0];
      w = weight(d1sq);
      vp->lambda() = w;

      // Update weights for all underlying particles.
      for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end(); ++itp) {
        Particle *at = system.storage->lookupLocalParticle(*itp);
        if (at) {
          at->lambda() = w;
        }
      }
    }
  }
}

// AdressNew Weighting function
real AdressNew::weight(real distanceSqr) {
  if (dex2 > distanceSqr) return 1.0;
  else if (dexdhy2 < distanceSqr) return 0.0;
  else {
    real argument = sqrt(distanceSqr) - dex;
    return 1.0 - (30.0 / (pow(dhy, 5.0))) * (1.0 / 5.0 * pow(argument, 5.0) - dhy / 2.0 * pow(argument, 4.0)
        + 1.0 / 3.0 * pow(argument, 3.0) * dhy * dhy);
    //return pow(cos(pidhy2 * argument),2.0); // for cosine squared weighting function
  }
}
real AdressNew::weightderivative(real distance) {
  real argument = distance - dex;
  return -(30.0 / (pow(dhy, 5.0)))
      * (pow(argument, 4.0) - 2.0 * dhy * pow(argument, 3.0) + argument * argument * dhy * dhy);
  //return -pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument); // for cosine squared weighting function
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void AdressNew::registerPython() {
  using namespace espressopp::python;

  class_<AdressNew, shared_ptr<AdressNew>, bases<Extension> >
      ("integrator_AdressNew", init<shared_ptr<System>, shared_ptr<FixedVSList> >())
      .def("connect", &AdressNew::connect)
      .def("setHyEx", &AdressNew::setHyEx)
      .def("disconnect", &AdressNew::disconnect);
}

}
}
