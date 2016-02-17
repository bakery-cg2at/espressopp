/*
  Copyright (C) 2015-2016
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

#include "DynamicResolution.hpp"
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include "python.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "FixedTupleListAdress.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
namespace integrator {
using namespace espressopp::iterator;  //NOLINT

DynamicResolution::DynamicResolution(
  shared_ptr<System> _system,
  shared_ptr<VerletListAdress> _verletList,
  shared_ptr<FixedTupleListAdress> _fixedtupleList,
  real _rate)
      : Adress(_system, _verletList, _fixedtupleList, false), rate_(_rate) {
    LOG4ESPP_INFO(theLogger, "construct DynamicResolution");

    type = Extension::Adress;

    // In this case the whole simulation box is treat as a single
    // resolution domain. Therefore there is no need to set the spatial dimensions.
    dhy = getSystem()->bc->getBoxL()[0];
    verletList->setHyEx(dhy, 0.0, true);
    verletList->setAdrCenter(dhy/2.0, 0.0, 0.0);
    verletList->setAdrRegionType(false);
    resolution_ = 0.0;
}


DynamicResolution::~DynamicResolution() {
  LOG4ESPP_INFO(theLogger, "~DynamicResolution");
  disconnect();
}

void DynamicResolution::connect() {
  Adress::connect();
  _aftIntV = integrator->aftIntV.connect(
      boost::bind(&DynamicResolution::ChangeResolution, this),
      boost::signals2::at_back);
}

void DynamicResolution::disconnect() {
  Adress::disconnect();
  _aftIntV.disconnect();
}

void DynamicResolution::set_active(bool active) {
  active_ = active;
}

void DynamicResolution::ChangeResolution() {
  if (!active_)
    return;

  // Increase the resolution lineary with the time.
  resolution_ += rate_;
  if (resolution_ > 1.0) {
    resolution_ = 1.0;
    active_ = false;
  } else if (resolution_ < 0.0) {
    resolution_ = 0.0;
    active_ = false;
  }
  updateWeights();
}

void DynamicResolution::SetPosVel() {
  System& system = getSystemRef();
  // Set the positions and velocity of CG particles & update weights.
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;
    FixedTupleListAdress::iterator it3;
    it3 = fixedtupleList->find(&vp);

    if (it3 != fixedtupleList->end()) {
      std::vector<Particle*> atList;
      atList = it3->second;

      // Compute center of mass
      Real3D cmp(0.0, 0.0, 0.0);  // center of mass position
      Real3D cmv(0.0, 0.0, 0.0);  // center of mass velocity

      for (std::vector<Particle*>::iterator it2 = atList.begin(); it2 != atList.end(); ++it2) {
        Particle &at = **it2;
        cmp += at.mass() * at.position();
        cmv += at.mass() * at.velocity();
        at.lambda() = resolution_;
        at.lambdaDeriv() = 0.0;
      }

      cmp /= vp.mass();
      cmv /= vp.mass();

      // Fold position.
      system.bc->foldPosition(cmp);

      // update (overwrite) the position and velocity of the VP
      vp.position() = cmp;
      vp.velocity() = cmv;

      vp.lambda() = resolution_;
      vp.lambdaDeriv() = 0.0;
    } else {
      std::stringstream msg;
      msg << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
      msg << " (" << vp.position() << ")";
      throw std::runtime_error(msg.str());
    }
  }
}


void DynamicResolution::updateWeights() {
  System& system = getSystemRef();
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;
    FixedTupleListAdress::iterator it3;
    it3 = fixedtupleList->find(&vp);

    vp.lambda() = resolution_;
    vp.lambdaDeriv() = 0.0;

    if (it3 != fixedtupleList->end()) {
      std::vector<Particle*> atList;
      atList = it3->second;

      // Propagate lambda/lambdaDeriv downstream to underlying atoms
      for (std::vector<Particle*>::iterator it2 = atList.begin(); it2 != atList.end(); ++it2) {
        Particle &at = **it2;
        at.lambda() = resolution_;
        at.lambdaDeriv() = 0.0;
      }
    } else {
      std::stringstream msg;
      msg << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
      msg << " (" << vp.position() << ")";
      throw std::runtime_error(msg.str());
    }
  }
}


/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void DynamicResolution::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<DynamicResolution, shared_ptr<DynamicResolution>, bases<Extension> >
    ("integrator_DynamicResolution", init<shared_ptr<System>, shared_ptr<VerletListAdress>,
                                          shared_ptr<FixedTupleListAdress>, real >())
    .add_property("active", &DynamicResolution::active, &DynamicResolution::set_active)
    .add_property("resolution", &DynamicResolution::resolution, &DynamicResolution::set_resolution)
    .add_property("rate", &DynamicResolution::rate, &DynamicResolution::set_rate)
    .def("connect", &DynamicResolution::connect)
    .def("disconnect", &DynamicResolution::disconnect)
    .def("SetPosVel", &DynamicResolution::SetPosVel);
}
}  // end namespace integrator
}  // end namespace espressopp
