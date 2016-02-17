/*
  Copyright (C) 2015-2016
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

#include "python.hpp"
#include "VerletListAdress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace espressopp::iterator;

  LOG4ESPP_LOGGER(VerletListAdress::theLogger, "VerletList");

  /*-------------------------------------------------------------*/

    VerletListAdress::VerletListAdress(shared_ptr<System> system, real cut, real _adrCut,
                                       bool rebuildVL, real _dEx, real _dHy)
    : SystemAccess(system), adrCut(_adrCut) {
      LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << cut);

      if (!system->storage) {
         throw std::runtime_error("system has no storage");
      }
      skin = system->getSkin();
      cutverlet = cut + skin;
      cutsq = cutverlet * cutverlet;
      builds = 0;

      // AdResS stuff
      adrCutverlet = adrCut + skin;
      adrcutsq = adrCutverlet*adrCutverlet;
      adrCenterSet = false;
      setHyEx(_dHy, _dEx, rebuildVL);
      fixedtupleList = system->storage->getFixedTuples();

      // make a connection to System to invoke rebuild on resort
      connectionResort = system->storage->onParticlesChanged.connect(
          boost::bind(&VerletListAdress::rebuild, this));
    }

    /*-------------------------------------------------------------*/
    
    void VerletListAdress::setHyEx(real _hy, real _ex, bool rebuildVL) {
      // AdResS stuff
      dEx = _ex;
      dHy = _hy;
      real adressSize = dEx + dHy + skin; // adress region size
      if (dEx + dHy == 0)
        adressSize = 0;
      adrsq = adressSize * adressSize;

      if (rebuildVL)
        rebuild(); // not called if exclusions are provided
    }


    void VerletListAdress::rebuild()
    {
      vlPairs.clear();
      adrZone.clear(); // particles in adress zone
      cgZone.clear(); // particles in CG zone
      adrPairs.clear(); // pairs in adress zone
      const bc::BC& bc = *getSystemRef().bc;

      // get local cells
      CellList localcells = getSystem()->storage->getLocalCells();

      // if adrCenter is not set, the center of adress zone moves along with some particles
      if (!adrCenterSet) { // maybe now working
          bool found_adr_position;
          LOG4ESPP_TRACE(theLogger, "found moving AdResS zone")
          for (CellListIterator it(localcells); it.isValid(); ++it) {
                found_adr_position = false;
                // loop over positions
                for (std::vector<Real3D*>::iterator it2 = adrPositions.begin(); 
                        it2 != adrPositions.end() && !found_adr_position; ++it2) {
                    Real3D dist;
                    real distsq;
                    bc.getMinimumImageVectorBox(dist, it->getPos(), **it2);

                    if (getAdrRegionType()){ // spherical adress region
                       distsq=dist.sqr();
                    }
                    else {  // slab-type adress region
                       distsq=dist[0]*dist[0];
                    }

                    //std::cout << "distance " << sqrt(distsq) << "\n";
                    if (distsq <= adrsq) {
                        adrZone.insert(&(*it));
                        // do not need to loop further
                        found_adr_position = true;
                    }
                }
                // if not near enough to any adrPositions, put in cgZone
                if (!found_adr_position) {
                    cgZone.insert(&(*it)); 
                }
          }
      }
      // center of adress zone is fixed
      else {
          for (CellListIterator it(localcells); it.isValid(); ++it) {
              Real3D dist;
              real distsq;
              bc.getMinimumImageVectorBox(dist, it->getPos(), adrCenter);
              if (getAdrRegionType()){ // spherical adress region
                distsq=dist.sqr();
              }
              else {  // slab-type adress region
                distsq=dist[0]*dist[0];
              }
              if (distsq <= adrsq) {
                  adrZone.insert(&(*it));
              }
              else {
                  cgZone.insert(&(*it));
              }
          }
      }

      // add particles to adress pairs and VL
      CellList cl = getSystem()->storage->getRealCells();
      int count=0;
      
      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " and " << it->second->id());
        checkPair(*it->first, *it->second);
        count+=1;
      }
      LOG4ESPP_INFO(theLogger, "rebuilt VerletList, cutsq = " << cutsq
                   << " local size = " << vlPairs.size()
                   << " Cells: " << cl.size()
                   << " count: " << count);
      builds++;
    }
    
    /// Checkes if pair of CG particles pt1 and pt2 is valid.
    /// If one of them is in HYB/AT zone then read underlying 
    /// atomistic particles and update adrPairs list taking into
    /// account exclusion lists for atomistic particles.
    void VerletListAdress::checkPair(Particle& pt1, Particle& pt2)
    {
      Real3D d;
      real distsq;

      // Adr zone stores atomistc pairs only for CG pairs that are in AT/HY zone.
      if (adrZone.count(&pt1) == 1 || adrZone.count(&pt2) == 1) {
        FixedTupleListAdress::iterator it1 = fixedtupleList->find(&pt1);
        FixedTupleListAdress::iterator it2 = fixedtupleList->find(&pt2);
        if (it1 != fixedtupleList->end() && it2 != fixedtupleList->end()) {
          for (std::vector<Particle*>::iterator itv1 = it1->second.begin();
               itv1 != it1->second.end(); ++itv1) {
            Particle &p3 = **itv1;
            for (std::vector<Particle*>::iterator itv2 = it2->second.begin();
                 itv2 != it2->second.end(); ++itv2) {
              Particle &p4 = **itv2;
              // Excluded list for atomistic pairs, the ids are unique so it can be the same list.
              if (exList.count(std::make_pair(p3.id(), p4.id())) == 1) continue;
              if (exList.count(std::make_pair(p4.id(), p3.id())) == 1) continue;
              d = p3.position() - p4.position();
              distsq = d.sqr();
              LOG4ESPP_TRACE(theLogger, "AT_1:" << p3.id()
                  << " of " << pt1.id() << " @ " << p3.position()
                  << " AT_2:" << p4.id() << " of " << pt2.id() << " @ "
                  << p4.position() << " -> distsq=" << distsq);
              if (distsq <= adrcutsq)
                adrPairs.add(p3, p4);
            }
          }
        }
      }

      // CG particles, use standard exclude list.
      if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
      if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;
      d = pt1.position() - pt2.position();
      distsq = d.sqr();
      if (distsq <= cutsq)
        vlPairs.add(pt1, pt2);

      LOG4ESPP_TRACE(theLogger, "CG p1: " << pt1.id()
                     << " @ " << pt1.position()
             << " - CG p2: " << pt2.id() << " @ " << pt2.position()
             << " -> distsq = " << distsq);
    }

    /*-------------------------------------------------------------*/

    int VerletListAdress::totalSize() const
    {
      System& system = getSystemRef();
      int size = vlPairs.size();
      int allsize;

      mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
      return allsize;
    }

    int VerletListAdress::totalAdrSize() const
    {
      System& system = getSystemRef();
      int size = adrPairs.size();
      int allsize;

      mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
      return allsize;
    }

    int VerletListAdress::localSize() const {
      return vlPairs.size();
    }

    int VerletListAdress::localAdrSize() const {
      return adrPairs.size();
    }

    python::tuple VerletListAdress::getPair(int i) {
      if (i > vlPairs.size())
        return python::make_tuple();
      else
        return python::make_tuple(vlPairs[i].first->id(), vlPairs[i].second->id());
    }

    python::tuple VerletListAdress::getAdrPair(int i) {
      if (i > adrPairs.size())
        return python::make_tuple();
      else
        return python::make_tuple(adrPairs[i].first->id(), adrPairs[i].second->id());
    }

    bool VerletListAdress::exclude(longint pid1, longint pid2) {
        exList.insert(std::make_pair(pid1, pid2));
        return true;
    }

    void VerletListAdress::addAdrParticle(longint pid) {
          std::cout<<"Warning! Moving adres region only works with VerletListAdressInteractionTemplate.hpp"<<std::endl;
          std::cout<<"VerletListHadressInteractionTemplate.hpp would need to be modified too"<<std::endl;
          adrList.insert(pid);
    }

    void VerletListAdress::setAdrCenter(real x, real y, real z){
        adrCenter = Real3D(x, y, z);
        adrCenterSet = true;
        adrPositions.push_back(&adrCenter);
    }

    void VerletListAdress::setAdrRegionType(bool _sphereAdr){
        sphereAdr = _sphereAdr;
        if (sphereAdr) {
          std::cout<<"Warning! Spherical adres region only works with VerletListAdressInteractionTemplate.hpp"<<std::endl;
          std::cout<<"VerletListHadressInteractionTemplate.hpp would need to be modified too"<<std::endl;
        }
    }

    bool VerletListAdress::getAdrRegionType(){
        return sphereAdr;
    }

    /* not used anymore
    // types above this number are considered atomistic
    void VerletListAdress::setAtType(size_t type) {
        atType = type;
    }*/

    /*-------------------------------------------------------------*/

    VerletListAdress::~VerletListAdress()
    {
      LOG4ESPP_INFO(theLogger, "~VerletList");

      if (!connectionResort.connected()) {
        connectionResort.disconnect();
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VerletListAdress::registerPython() {
      using namespace espressopp::python;

      bool (VerletListAdress::*pyExclude)(longint pid1, longint pid2)
            = &VerletListAdress::exclude;

      void (VerletListAdress::*pyAddAdrParticle)(longint pid)
            = &VerletListAdress::addAdrParticle;

      void (VerletListAdress::*pySetAdrCenter)(real x, real y, real z)
                  = &VerletListAdress::setAdrCenter;

      void (VerletListAdress::*pySetAdrRegionType)(bool _sphereAdr)
            = &VerletListAdress::setAdrRegionType;

      class_<VerletListAdress, shared_ptr<VerletList> >
        ("VerletListAdress", init< shared_ptr<System>, real, real, bool, real, real>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletListAdress::getBuilds, &VerletListAdress::setBuilds)
        .def("totalSize", &VerletListAdress::totalSize)
        .def("totalAdrSize", &VerletListAdress::totalAdrSize)
        .def("localSize", &VerletListAdress::localSize)
        .def("localAdrSize", &VerletListAdress::localAdrSize)
        .def("exclude", pyExclude)
        .def("getPair", &VerletListAdress::getPair)
        .def("getAdrPair", &VerletListAdress::getAdrPair)
        .def("addAdrParticle", pyAddAdrParticle)
        .def("setAdrCenter", pySetAdrCenter)
        .def("setAdrRegionType", pySetAdrRegionType)
        .def("rebuild", &VerletListAdress::rebuild)
        //.def("setAtType", pySetAtType)
        ;
    }

}
