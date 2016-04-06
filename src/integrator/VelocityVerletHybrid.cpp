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

//#include <iomanip>
#include "python.hpp"
#include "VelocityVerletHybrid.hpp"
#include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"
#include "bc/BC.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER(name)
#endif

namespace espressopp {
using namespace std;
namespace integrator {
using namespace interaction;
using namespace iterator;
using namespace esutil;

LOG4ESPP_LOGGER(VelocityVerletHybrid::theLogger, "VelocityVerletHybrid");

VelocityVerletHybrid::VelocityVerletHybrid(shared_ptr<System> system, shared_ptr<FixedVSList> vs_list)
    : MDIntegrator(system), vs_list_(vs_list) {
  LOG4ESPP_INFO(theLogger, "construct VelocityVerletHybrid");
  resortFlag = true;
  maxDist = 0.0;
  timeIntegrate.reset();
  resetTimers();
  System &sys = getSystemRef();
}

VelocityVerletHybrid::~VelocityVerletHybrid() {
  LOG4ESPP_INFO(theLogger, "free VelocityVerletHybrid");
}

real VelocityVerletHybrid::updateVS() {
  LOG4ESPP_INFO(theLogger, "update position and velocity of VS");
  FixedVSList::GlobalTuples vs = vs_list_->globalTuples;
  FixedVSList::GlobalTuples::iterator it = vs.begin();

  System &system = getSystemRef();
  storage::Storage &storage = *system.storage;
  const bc::BC& bc = *getSystemRef().bc;

  real maxSqDist = 0.0;

  for (; it != vs.end(); ++it) {
    Particle *vp = storage.lookupRealParticle(it->first);
    if (vp) {
      Real3D cmp(0.0, 0.0, 0.0);
      Real3D cmv(0.0, 0.0, 0.0);
      for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end(); ++itp) {
        Particle *at = storage.lookupLocalParticle(*itp);
        if (at) {
          LOG4ESPP_DEBUG(theLogger, "vp-" << vp->id() << " at:" << at->id() << " " << at->position());
          cmp += at->mass() * at->position();
          cmv += at->mass() * at->velocity();
        } else {
          std::cout << " AT particle (" << *itp << ") of VP " << vp->id() << "-"
                    << vp->ghost() << " not found in tuples ";
          std::cout << " (" << vp->position() << ")" << std::endl;
          exit(1);
        }
      }
      cmp /= vp->mass();
      cmv /= vp->mass();
      LOG4ESPP_DEBUG(theLogger, "vp-" << vp->id() << " cmp=" << cmp);

      Real3D before_cmp = cmp;
      bc.getMinimumDistance(cmp);

      Real3D d_p = vp->position() - cmp;
      maxSqDist = std::max(maxSqDist, d_p.sqr());
      LOG4ESPP_DEBUG(theLogger, "vp " << vp->id() << " old=" << vp->position() << " new=" << cmp
                                << " before_cmp=" << before_cmp);

      vp->position() = cmp;
      vp->velocity() = cmv;
    }
  }

  real maxAllSqDist;
  mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

  LOG4ESPP_INFO(theLogger, " particles in updateVS" <<
      ", max move local = " << sqrt(maxSqDist) <<
      ", global = " << sqrt(maxAllSqDist));

  return sqrt(maxAllSqDist);
}

void VelocityVerletHybrid::updateVS_vel() {
  LOG4ESPP_INFO(theLogger, "update velocity of VS");

  FixedVSList::GlobalTuples vs = vs_list_->globalTuples;
  FixedVSList::GlobalTuples::iterator it = vs.begin();

  System &system = getSystemRef();
  storage::Storage &storage = *system.storage;

  for (; it != vs.end(); ++it) {
    Particle *vp = storage.lookupRealParticle(it->first);
    if (vp) {
      Real3D cmv(0.0, 0.0, 0.0);
      for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end(); ++itp) {
        Particle *at = storage.lookupLocalParticle(*itp);
        if (at) {
          cmv += at->mass() * at->velocity();
        } else {
          std::cout << " AT particle (" << *itp << ") of VP " << vp->id() << "-"
              << vp->ghost() << " not found in tuples ";
          std::cout << " (" << vp->position() << ")" << std::endl;
          exit(1);
        }
      }
      // Updates velocity.
      cmv /= vp->mass();
      vp->velocity() = cmv;
    }
  }
}

void VelocityVerletHybrid::run(int nsteps) {
  int nResorts = 0;
  real time;
  timeIntegrate.reset();
  System &system = getSystemRef();
  storage::Storage &storage = *system.storage;
  skinHalf = 0.5 * system.getSkin();

  time = timeIntegrate.getElapsedTime();
  // signal
  runInit();
  timeRunInitS += timeIntegrate.getElapsedTime() - time;

  // Before start make sure that particles are on the right processor
  if (resortFlag) {
    LOG4ESPP_INFO(theLogger, "resort particles");
    storage.decompose();
    maxDist = 0.0;
    resortFlag = false;
  }

  bool recalcForces = true;  // TODO: more intelligent

  if (recalcForces) {
    LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

    time = timeIntegrate.getElapsedTime();
    // signal
    recalc1();
    timeRecalc1S += timeIntegrate.getElapsedTime() - time;

    updateForces();

    time = timeIntegrate.getElapsedTime();
    // signal
    recalc2();
    timeRecalc2S += timeIntegrate.getElapsedTime() - time;
  }

  LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");

  for (int i = 0; i < nsteps; i++) {
    LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

    time = timeIntegrate.getElapsedTime();
    // signal
    befIntP();
    timeBefIntPS += timeIntegrate.getElapsedTime() - time;

    LOG4ESPP_INFO(theLogger, "updating positions and velocities")
    maxDist += integrate1();
    timeInt1 += timeIntegrate.getElapsedTime() - time;

    time = timeIntegrate.getElapsedTime();
    // signal
    aftIntP();
    timeAftIntPS += timeIntegrate.getElapsedTime() - time;

    LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

    if (maxDist > skinHalf)
      resortFlag = true;

    if (resortFlag) {
      time = timeIntegrate.getElapsedTime();
      LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
      storage.decompose();
      maxDist = 0.0;
      resortFlag = false;
      nResorts++;
      timeResort += timeIntegrate.getElapsedTime() - time;
    }

    LOG4ESPP_INFO(theLogger, "updating forces")
    updateForces();

    timeIntegrate.startMeasure();
    // signal
    befIntV();
    timeBefIntVS += timeIntegrate.stopMeasure();

    time = timeIntegrate.getElapsedTime();
    integrate2();
    timeInt2 += timeIntegrate.getElapsedTime() - time;

    timeIntegrate.startMeasure();
    // signal
    aftIntV();
    timeAftIntVS += timeIntegrate.stopMeasure();
  }

  LOG4ESPP_INFO(theLogger, "finished run");
}

void VelocityVerletHybrid::resetTimers() {

}


void VelocityVerletHybrid::printTimers() {

}

void VelocityVerletHybrid::loadTimers(std::vector<real> &return_vector) {

}

static boost::python::object wrapGetTimers(class VelocityVerletHybrid *obj) {

}


real VelocityVerletHybrid::integrate1() {
  System &system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  // loop over all particles of the local cells
  int count = 0;
  real maxSqDist = 0.0; // maximal square distance a particle moves
  for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    if (cit->vp())  // propagate only real particles, skip virtual sites.
      continue;
    real sqDist = 0.0;
    LOG4ESPP_INFO(theLogger, "updating first half step of velocities and full step of positions")
    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() <<
        ", pos = " << cit->position() <<
        ", v = " << cit->velocity() <<
        ", f = " << cit->force());

    real dtfm = 0.5 * dt / cit->mass();

    // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
    cit->velocity() += dtfm * cit->force();

    // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
    Real3D deltaP = cit->velocity();

    deltaP *= dt;
    cit->position() += deltaP;
    sqDist += deltaP * deltaP;

    count++;

    maxSqDist = std::max(maxSqDist, sqDist);
  }

  // signal
  inIntP(maxSqDist);

  real maxAllSqDist;
  mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

  LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
      ", max move local = " << sqrt(maxSqDist) <<
      ", global = " << sqrt(maxAllSqDist));

  return sqrt(maxAllSqDist);
}

void VelocityVerletHybrid::integrate2() {
  LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
  System &system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  // loop over all particles of the local cells
  real half_dt = 0.5 * dt;
  for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    if (cit->vp())
      continue;
    real dtfm = half_dt / cit->mass();
    /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
    cit->velocity() += dtfm * cit->force();
  }

  // Update velocity of VS based on the AT velocities.
  updateVS_vel();

  step++;
}

void VelocityVerletHybrid::calcForces() {
  LOG4ESPP_INFO(theLogger, "calculate forces");

  initForces();

  timeIntegrate.startMeasure();
  // signal
  aftInitF();
  timeAftInitFS += timeIntegrate.stopMeasure();

  System &sys = getSystemRef();
  const InteractionList &srIL = sys.shortRangeInteractions;
  real time;
  for (size_t i = 0; i < srIL.size(); i++) {
    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
    srIL[i]->addForces();
  }
}

void VelocityVerletHybrid::updateForces() {
  LOG4ESPP_INFO(theLogger, "update ghosts, calculate forces and collect ghost forces")
  real time;
  storage::Storage &storage = *getSystemRef().storage;

  // Make sure that positions and velocity of VS sites is correct with respect to the atoms.
  real maxDist = updateVS();

  if (maxDist > skinHalf) {
    time = timeIntegrate.getElapsedTime();
    storage.decompose();
    timeComm1 += timeIntegrate.getElapsedTime() - time;
  } else {
    time = timeIntegrate.getElapsedTime();
    storage.updateGhosts();
    timeComm1 += timeIntegrate.getElapsedTime() - time;
  }

  time = timeIntegrate.getElapsedTime();
  calcForces();
  timeForce += timeIntegrate.getElapsedTime() - time;

  time = timeIntegrate.getElapsedTime();
  storage.collectGhostForces();
  timeComm2 += timeIntegrate.getElapsedTime() - time;

  // Distribute forces from VS to AT.
  distributeVSforces();

  timeIntegrate.startMeasure();
  // signal
  aftCalcF();
  timeAftCalcFS += timeIntegrate.stopMeasure();
}

void VelocityVerletHybrid::initForces() {
  // forces are initialized for real + ghost particles

  System &system = getSystemRef();
  CellList localCells = system.storage->getLocalCells();

  LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

  for (CellListIterator cit(localCells); !cit.isDone(); ++cit) {
    cit->force() = 0.0;
    cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
  }
}

void VelocityVerletHybrid::distributeVSforces() {
  // Zeros forces on ghost particles.
  System &system = getSystemRef();
  storage::Storage &storage = *system.storage;
  CellList ghostCells = storage.getGhostCells();
  LOG4ESPP_INFO(theLogger, "zeros forces on ghost particles");

  for (CellListIterator cit(ghostCells); !cit.isDone(); ++cit) {
    cit->force() = 0.0;
    cit->drift() = 0.0;
  }

  // Distribute forces to AT particles from VS.
  FixedVSList::GlobalTuples vs = vs_list_->globalTuples;
  FixedVSList::GlobalTuples::iterator it = vs.begin();

  for (; it != vs.end(); ++it) {
    Particle *vp = storage.lookupRealParticle(it->first);
    if (vp) {
      for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end(); ++itp) {
        Particle *at = storage.lookupLocalParticle(*itp);
        if (at) {
          at->force() += (at->mass() / vp->mass()) * vp->force();
        } else {
          std::cout << " AT particle (" << *itp << ") of VP " << vp->id() << "-"
              << vp->ghost() << " not found in tuples ";
          std::cout << " (" << vp->position() << ")" << std::endl;
          exit(1);
        }
      }
    }
  }
  real time = timeIntegrate.getElapsedTime();
  storage.collectGhostForces();
  timeComm2 += timeIntegrate.getElapsedTime() - time;
}

void VelocityVerletHybrid::printForces(bool withGhosts) {

}

void VelocityVerletHybrid::printPositions(bool withGhosts) {

}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void VelocityVerletHybrid::registerPython() {

  using namespace espressopp::python;

  // Note: use noncopyable and no_init for abstract classes
  class_<VelocityVerletHybrid, bases<MDIntegrator>, boost::noncopyable>
      ("integrator_VelocityVerletHybrid", init<shared_ptr<System>, shared_ptr<FixedVSList> >())
      .def("getTimers", &wrapGetTimers)
      .def("resetTimers", &VelocityVerletHybrid::resetTimers);
}
}
}
