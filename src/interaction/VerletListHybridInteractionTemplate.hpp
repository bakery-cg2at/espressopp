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

// ESPP_CLASS 
#ifndef _INTERACTION_VERLETLISTHYBRIDINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTHYBRIDINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"

#include "storage/Storage.hpp"

namespace espressopp {
namespace interaction {
template<typename _Potential>
class VerletListHybridInteractionTemplate: public Interaction {

 protected:
  typedef _Potential Potential;

 public:
  VerletListHybridInteractionTemplate(shared_ptr<VerletList> _verletList, bool _cg_potential)
      : verletList(_verletList), cgPotential(_cg_potential) {
    potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
    ntypes = 0;
    scaleFactor_ = 1.0;
    LOG4ESPP_DEBUG(_Potential::theLogger, "initialize hybrid potential cg=" << cgPotential);
  }

  virtual ~VerletListHybridInteractionTemplate() { };

  void setVerletList(shared_ptr<VerletList> _verletList) {
    verletList = _verletList;
  }

  shared_ptr<VerletList> getVerletList() {
    return verletList;
  }

  void setPotential(int type1, int type2, const Potential &potential) {
    // typeX+1 because i<ntypes
    ntypes = std::max(ntypes, std::max(type1 + 1, type2 + 1));
    potentialArray.at(type1, type2) = potential;
    LOG4ESPP_INFO(_Potential::theLogger, "added potential for type1=" << type1 << " type2=" << type2);
    if (type1 != type2) { // add potential in the other direction
      potentialArray.at(type2, type1) = potential;
      LOG4ESPP_INFO(_Potential::theLogger,
                    "automatically added the same potential for type1=" << type2 << " type2=" << type1);
    }
  }

  // this is used in the innermost force-loop
  Potential &getPotential(int type1, int type2) {
    return potentialArray.at(type1, type2);
  }

  // this is mainly used to access the potential from Python (e.g. to change parameters of the potential)
  shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
    return make_shared<Potential>(potentialArray.at(type1, type2));
  }

  void setScaleFactor(real s) {
    scaleFactor_ = s;
    if (scaleFactor_ >= 1.0)
      scaleFactor_ = 1.0;
    else if (scaleFactor_ <= 0.0)
      scaleFactor_ = 0.0;
  }

  real scaleFactor() { return scaleFactor_; }

  virtual void addForces();
  virtual real computeEnergy();
  virtual real computeEnergyAA();
  virtual real computeEnergyCG();
  virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
  virtual real computeVirial();
  virtual void computeVirialTensor(Tensor &w);
  virtual void computeVirialTensor(Tensor &w, real z);
  virtual void computeVirialTensor(Tensor *w, int n);
  virtual real getMaxCutoff();
  virtual int bondType() { return Nonbonded; }
  virtual python::list getInteractionMatrix();

 protected:
  int ntypes;
  shared_ptr<VerletList> verletList;
  esutil::Array2D<Potential, esutil::enlarge> potentialArray;

  bool cgPotential;
  real scaleFactor_;
  // not needed esutil::Array2D<shared_ptr<Potential>, esutil::enlarge> potentialArrayPtr;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template<typename _Potential>
inline void
VerletListHybridInteractionTemplate<_Potential>::
addForces() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");

  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();

    real p1lambda = p1.lambda();
    real p2lambda = p2.lambda();
    if (p1lambda < 0.0 || p2lambda < 0.0)
      continue;

    real w12 = p1lambda * p2lambda;
    real forcescale12 = w12;
    if (cgPotential) {
      forcescale12 = (1 - w12);
    }

    forcescale12 *= scaleFactor_;

    if (forcescale12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      Real3D force(0.0);
      if (potential._computeForce(force, p1, p2)) {
        p1.force() += forcescale12 * force;
        p2.force() -= forcescale12 * force;
        LOG4ESPP_TRACE(_Potential::theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " force=" << force);
      }
    }
  }
}

template<typename _Potential>
inline real
VerletListHybridInteractionTemplate<_Potential>::
computeEnergy() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up potential energies");

  real es = 0.0;
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();

    real p1lambda = p1.lambda();
    real p2lambda = p2.lambda();
    if (p1lambda < 0.0 || p2lambda < 0.0)
      continue;

    real w12 = p1lambda * p2lambda;
    real forcescale12 = w12;
    if (cgPotential) {
      forcescale12 = (1 - w12);
    }

    forcescale12 *= scaleFactor_;

    if (forcescale12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      es += forcescale12*potential._computeEnergy(p1, p2);
      LOG4ESPP_TRACE(_Potential::theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " potential energy=" << es);
    }
  }

  // reduce over all CPUs
  real esum = 0.0;
  boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, es, esum, std::plus<real>());
  return esum;
}

template<typename _Potential>
inline real
VerletListHybridInteractionTemplate<_Potential>::
computeEnergyAA() {
  return computeEnergy();
}

template<typename _Potential>
inline real
VerletListHybridInteractionTemplate<_Potential>::
computeEnergyCG() {
  return computeEnergy();
}

template<typename _Potential>
inline void
VerletListHybridInteractionTemplate<_Potential>::
computeVirialX(std::vector<real> &p_xx_total, int bins) {
  LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialX() is not yet implemented.");
}

template<typename _Potential>
inline real
VerletListHybridInteractionTemplate<_Potential>::
computeVirial() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial");

  real w = 0.0;
  for (PairList::Iterator it(verletList->getPairs());
       it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();

    real p1lambda = p1.lambda();
    real p2lambda = p2.lambda();
    if (p1lambda < 0.0 || p2lambda < 0.0)
      continue;

    real w12 = p1lambda * p2lambda;
    real forcescale12 = w12;
    if (cgPotential) {
      forcescale12 = (1 - w12);
    }

    forcescale12 *= scaleFactor_;

    if (forcescale12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      // shared_ptr<Potential> potential = getPotential(type1, type2);

      Real3D force(0.0, 0.0, 0.0);
      if (potential._computeForce(force, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
        Real3D r21 = p1.position() - p2.position();
        w = w + r21 * forcescale12*force;
      }
    }
  }

  // reduce over all CPUs
  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return wsum;
}

template<typename _Potential>
inline void
VerletListHybridInteractionTemplate<_Potential>::
computeVirialTensor(Tensor &w) {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial tensor");

  Tensor wlocal(0.0);
  for (PairList::Iterator it(verletList->getPairs());
       it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();

    real p1lambda = p1.lambda();
    real p2lambda = p2.lambda();
    if (p1lambda < 0.0 || p2lambda < 0.0)
      continue;

    real w12 = p1lambda * p2lambda;
    real forcescale12 = w12;
    if (cgPotential) {
      forcescale12 = (1 - w12);
    }

    forcescale12 *= scaleFactor_;

    if (forcescale12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      // shared_ptr<Potential> potential = getPotential(type1, type2);

      Real3D force(0.0, 0.0, 0.0);
      if (potential._computeForce(force, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
        Real3D r21 = p1.position() - p2.position();
        wlocal += Tensor(r21, forcescale12*force);
      }
    }
  }

  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}

// local pressure tensor for layer, plane is defined by z coordinate
template<typename _Potential>
inline void
VerletListHybridInteractionTemplate<_Potential>::
computeVirialTensor(Tensor &w, real z) {
  LOG4ESPP_ERROR(_Potential::theLogger, "computeVirialTensor not implemented for VerletListHybridInteraction");
}

// it will calculate the pressure in 'n' layers along Z axis
// the first layer has coordinate 0.0 the last - (Lz - Lz/n)
template<typename _Potential>
inline void
VerletListHybridInteractionTemplate<_Potential>::
computeVirialTensor(Tensor *w, int n) {
  LOG4ESPP_ERROR(_Potential::theLogger, "computeVirialTensor not implemented for VerletListHybridInteraction");
}

template<typename _Potential>
inline real
VerletListHybridInteractionTemplate<_Potential>::
getMaxCutoff() {
  real cutoff = 0.0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
      cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
    }
  }
  return cutoff;
}

template<typename _Potential>
inline boost::python::list
VerletListHybridInteractionTemplate<_Potential>::getInteractionMatrix() {
  python::list params;

  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
      params.append(python::make_tuple(i, j, getPotential(i, j).getParams()));
    }
  }

  return params;
}

}
}
#endif
