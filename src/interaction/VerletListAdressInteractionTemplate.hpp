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
#ifndef _INTERACTION_VERLETLISTADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTADRESSINTERACTIONTEMPLATE_HPP

//#include <typeinfo>

#include "System.hpp"
#include "bc/BC.hpp"

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"
#include "esutil/Array2D.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _PotentialAT, typename _PotentialCG >
    class VerletListAdressInteractionTemplate: public Interaction {
    
    protected:
      typedef _PotentialAT PotentialAT;
      typedef _PotentialCG PotentialCG;
    
    public:
      VerletListAdressInteractionTemplate
      (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList)
                : verletList(_verletList), fixedtupleList(_fixedtupleList) {

          potentialArrayAT = esutil::Array2D<PotentialAT, esutil::enlarge>(0, 0, PotentialAT());
          potentialArrayCG = esutil::Array2D<PotentialCG, esutil::enlarge>(0, 0, PotentialCG());

          // AdResS stuff
          pidhy2 = M_PI/(verletList->getHy() * 2);
          dex = verletList->getEx();
          dex2 = dex * dex;
          dhy = verletList->getHy();
          dexdhy = dex + dhy;
          dexdhy2 = dexdhy * dexdhy;
          
          ntypes = 0;
      }

      void
      setVerletList(shared_ptr < VerletListAdress > _verletList) {
        verletList = _verletList;
      }

      shared_ptr<VerletListAdress> getVerletList() {
        return verletList;
      }

      void
      setFixedTupleList(shared_ptr<FixedTupleListAdress> _fixedtupleList) {
          fixedtupleList = _fixedtupleList;
      }

      void
      setPotentialAT(int type1, int type2, const PotentialAT &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        
          potentialArrayAT.at(type1, type2) = potential;
          if (type1 != type2) { // add potential in the other direction
             potentialArrayAT.at(type2, type1) = potential;
          }
      }

      void
      setPotentialCG(int type1, int type2, const PotentialCG &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        
         potentialArrayCG.at(type1, type2) = potential;
         if (type1 != type2) { // add potential in the other direction
             potentialArrayCG.at(type2, type1) = potential;
         }
      }

      PotentialAT &getPotentialAT(int type1, int type2) {
        return potentialArrayAT.at(type1, type2);
      }

      PotentialCG &getPotentialCG(int type1, int type2) {
        return potentialArrayCG.at(type1, type2);
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();      
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Nonbonded; }
      virtual python::list getInteractionMatrix();

    protected:
      int ntypes;
      shared_ptr<VerletListAdress> verletList;
      shared_ptr<FixedTupleListAdress> fixedtupleList;
      esutil::Array2D<PotentialAT, esutil::enlarge> potentialArrayAT;
      esutil::Array2D<PotentialCG, esutil::enlarge> potentialArrayCG;

      // AdResS stuff
      real pidhy2; // pi / (dhy * 2)
      real dexdhy; // dex + dhy
      real dexdhy2; // dexdhy^2
      real dex;
      real dhy;
      real dex2; // dex^2

    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      // CG pairs, both in CG and HY/AT zone.
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        real w12 = p1.lambda() * p2.lambda();
        
        if (w12 < ALMOST_ONE) {
          const PotentialCG &potentialCG = getPotentialCG(type1, type2);
          // CG forces
          Real3D force(0.0, 0.0, 0.0);
          if(potentialCG._computeForce(force, p1, p2)) {
            force *= (1-w12);
            p1.force() += force;
            p2.force() -= force;
          }
        }
      }  // end calculating forces in CG region

      // Compute forces of atomistic particles in HY/AT zone.
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {

         // These are the two AT interacting sites.
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;

         // read weights
         real w12 = p1.lambda() * p2.lambda();
         
         if (w12 > ALMOST_ZERO) {
           const PotentialAT &potentialAT = getPotentialAT(p1.type(), p2.type());
           Real3D force(0.0, 0.0, 0.0);
           if (potentialAT._computeForce(force, p1, p2)) {
             force *= w12;
             p1.force() += force;
             p2.force() -= force;
           }
         }
      }

      // Loop over CG particles and overwrite AT forces and velocity.
      // This makes the AT particles move along with CG particles.
      
      // Note (Karsten): This is a different approach compared to H-AdResS. Here, we overwrite intra-atomistic
      // rotations and vibrations in the CG zone. This leads to failures in the kinetic energy. However, in Force-AdResS there is no energy conservation anyway.
      // Here we calculate CG forces/velocities and distribute them to AT particles. In contrast, in H-AdResS, we calculate AT forces from intra-molecular
      // interactions and inter-molecular center-of-mass interactions and just update the positions of the center-of-mass CG particles.
      std::set<Particle*> cgZone = verletList->getCGZone();    
      for (std::set<Particle*>::iterator it=cgZone.begin();
                    it != cgZone.end(); ++it) {

            Particle &vp = **it;

            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);

            if (it3 != fixedtupleList->end()) {

                std::vector<Particle*> atList1;
                atList1 = it3->second;

                Real3D vpfm = vp.force() / vp.getMass();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                        itv != atList1.end(); ++itv) {
                    Particle &at = **itv;
                    at.velocity() = vp.velocity(); // Overwrite velocity - Note (Karsten): See comment above.
                    at.force() += at.mass() * vpfm;
                }
            }
            else { // this should not happen
                std::cout << " VP particle not found in tuples: " << vp.id() << "-" << vp.ghost();
                exit(1);
                return;
            }
      }
    }
      
    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergy() {
 
      LOG4ESPP_INFO(theLogger, "compute energy of the Verlet list pairs");
      
      real e = 0.0;        
      for (PairList::Iterator it(verletList->getPairs()); 
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          int type1 = p1.type();
          int type2 = p2.type();
          real w12 = p1.lambda() * p2.lambda();
          if (w12 < ALMOST_ONE) {
            const PotentialCG &potential = getPotentialCG(type1, type2);
            e += (1.0 - w12) * potential._computeEnergy(p1, p2);
          }
      }
      
      for (PairList::Iterator it(verletList->getAdrPairs()); 
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second; 
          real w12 = p1.lambda() * p2.lambda();
          int type1 = p1.type();
          int type2 = p2.type();
          if (w12 > ALMOST_ZERO) {
            const PotentialAT &potentialAT = getPotentialAT(p1.type(), p2.type());
            e += w12*potentialAT._computeEnergy(p1, p2);
          }
      }
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      
    }
    
    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyAA() {
      real e = 0.0;        

      for (PairList::Iterator it(verletList->getAdrPairs()); 
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;                
          real w12 = p1.lambda() * p2.lambda();
          if (w12 > ALMOST_ZERO) {
            // AT energies
            const PotentialAT &potentialAT = getPotentialAT(p1.type(), p2.type());
            e += w12*potentialAT._computeEnergy(p1, p2);
          }
      }
       
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      

    }
    
    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyCG() {
      real e = 0.0;        
      for (PairList::Iterator it(verletList->getPairs()); 
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          real w12 = p1.lambda() * p2.lambda();
          if (w12 < ALMOST_ONE) {
            int type1 = p1.type();
            int type2 = p2.type();
            const PotentialCG &potential = getPotentialCG(type1, type2);
            e += (1.0-w12)*potential._computeEnergy(p1, p2);
          }
      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      
    }

    template < typename _PotentialAT, typename _PotentialCG >
    inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      //std::cout << "Warning! At the moment computeVirialX in VerletListAdressInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    
      std::set<Particle*> cgZone = verletList->getCGZone();
      for (std::set<Particle*>::iterator it=cgZone.begin();
          it != cgZone.end(); ++it) {

      Particle &vp = **it;
      vp.lambda() = 0.0;
      //weights.insert(std::make_pair(&vp, 0.0));
      }
      
      std::set<Particle*> adrZone = verletList->getAdrZone();
      for (std::set<Particle*>::iterator it=adrZone.begin();
              it != adrZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);

          if (it3 != fixedtupleList->end()) {

              std::vector<Particle*> atList;
              atList = it3->second;

              // compute center of mass
              Real3D cmp(0.0, 0.0, 0.0); // center of mass position
              Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
              //real M = vp.getMass(); // sum of mass of AT particles
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle &at = **it2;
                  //Real3D d1 = at.position() - vp.position();
                  //Real3D d1;
                  //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                  //cmp += at.mass() * d1;

                  cmp += at.mass() * at.position();
                  cmv += at.mass() * at.velocity();
              }
              cmp /= vp.getMass();
              cmv /= vp.getMass();
              //cmp += vp.position(); // cmp is a relative position
              //std::cout << " cmp M: "  << M << "\n\n";
              //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
              vp.velocity() = cmv;

              // calculate distance to nearest adress particle or center
              std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
              Real3D pa = **it2; // position of adress particle
              Real3D d1 = vp.position() - pa;
              //real d1 = vp.position()[0] - pa[0];
              real min1sq = d1.sqr(); // set min1sq before loop  // d1*d1;
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                   pa = **it2;
                   d1 = vp.position() - pa;
                   //d1 = vp.position()[0] - pa[0];
                   real distsq1 = d1.sqr(); // d1*d1;
                   //std::cout << pa << " " << sqrt(distsq1) << "\n";
                   if (distsq1 < min1sq) min1sq = distsq1;
              }

              //real min1 = sqrt(min1sq);
              //std::cout << vp.id() << " min: " << min1 << "\n";
              //std::cout << vp.id() << " dex: " << dex << "\n";
              //std::cout << vp.id() << " dex+dhy: " << dexdhy << "\n";

              // calculate weight and write it in the map
              real w;
              if (dex2 > min1sq) w = 1;
              else if (dexdhy2 < min1sq) w = 0;
              else {
                   //w = cos(pidhy2 * (sqrt(min1sq) - dex));
                   //w *= w;
                   real argument = sqrt(min1sq) - dex;
                   w = 1.0-(30.0/(pow(dhy, 5.0)))*(1.0/5.0*pow(argument, 5.0)-dhy/2.0*pow(argument, 4.0)+1.0/3.0*pow(argument, 3.0)*dhy*dhy);
              }

              vp.lambda() = w;
              //weights.insert(std::make_pair(&vp, w));

              //if (w1 == 1 || w2 == 1) std::cout << p1.id() << " ";
              //std::cout << vp.id() << " weight: " << w << "\n";
          }
          else { // this should not happen
              std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
              std::cout << " (" << vp.position() << ")\n";
              exit(1);
              return;
          }
      }
      
      int i = 0;
      int bin1 = 0;
      int bin2 = 0;
      
      System& system = verletList->getSystemRef();
      Real3D Li = system.bc->getBoxL();
      real Delta_x = Li[0] / (real)bins;
      real Volume = Li[1] * Li[2] * Delta_x;
      size_t size = bins;
      std::vector <real> p_xx_local(size);      

      for (i = 0; i < bins; ++i)
        {
          p_xx_local[i] = 0.0;
        }
      
      // Calculate CG particles.
      for (PairList::Iterator it(verletList->getPairs());                
           it.isValid(); ++it) {                                         
        Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.type();                                           
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          real vir_temp = 0.5 * dist[0] * force[0];
          
          if (p1.position()[0] > Li[0])
          {
              real p1_wrap = p1.position()[0] - Li[0];
              bin1 = floor (p1_wrap / Delta_x);    
          }         
          else if (p1.position()[0] < 0.0)
          {
              real p1_wrap = p1.position()[0] + Li[0];
              bin1 = floor (p1_wrap / Delta_x);    
          }
          else
          {
              bin1 = floor (p1.position()[0] / Delta_x);          
          }
          
          if (p2.position()[0] > Li[0])
          {
              real p2_wrap = p2.position()[0] - Li[0];
              bin2 = floor (p2_wrap / Delta_x);     
          }         
          else if (p2.position()[0] < 0.0)
          {
              real p2_wrap = p2.position()[0] + Li[0];
              bin2 = floor (p2_wrap / Delta_x);     
          }
          else
          {
              bin2 = floor (p2.position()[0] / Delta_x);          
          }             
         
          if (bin1 >= p_xx_local.size() || bin2 >= p_xx_local.size()){
              std::cout << "p_xx_local.size() " << p_xx_local.size() << "\n";
              std::cout << "bin1 " << bin1 << " bin2 " << bin2 << "\n";
              std::cout << "p1.position()[0] " << p1.position()[0] << " p2.position()[0]" << p2.position()[0] << "\n";
              std::cout << "FATAL ERROR: computeVirialX error" << "\n";
              exit(0);
          }
          
          p_xx_local.at(bin1) += vir_temp;
          p_xx_local.at(bin2) += vir_temp;
        }
      }
      
      // Calculate for atomistic particles.
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
         real w1, w2;
         // these are the two VP interacting
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
     
         Real3D dist = p1.position() - p2.position();
         w1 = p1.lambda();         
         w2 = p2.lambda();
         real w12 = w1*w2;  // AdResS

         if (p1.position()[0] > Li[0])
         {
             real p1_wrap = p1.position()[0] - Li[0];
             bin1 = floor (p1_wrap / Delta_x);    
         }         
         else if (p1.position()[0] < 0.0)
         {
             real p1_wrap = p1.position()[0] + Li[0];
             bin1 = floor (p1_wrap / Delta_x);    
         }
         else
         {
             bin1 = floor (p1.position()[0] / Delta_x);          
         }
 
         if (p2.position()[0] > Li[0])
         {
             real p2_wrap = p2.position()[0] - Li[0];
             bin2 = floor (p2_wrap / Delta_x);     
         }         
         else if (p2.position()[0] < 0.0)
         {
             real p2_wrap = p2.position()[0] + Li[0];
             bin2 = floor (p2_wrap / Delta_x);     
         }
         else
         {
             bin2 = floor (p2.position()[0] / Delta_x);          
         }
         
         // force between VP particles
       
         // force between AT particles
         if (w12 > ALMOST_ZERO) { // calculate AT force if both VP are outside CG region (HY-HY, HY-AT, AT-AT)
           int type1 = p1.type();
           int type2 = p2.type();
           Real3D force(0.0, 0.0, 0.0);
                 
           // AT forces
           const PotentialAT &potentialAT = getPotentialAT(p1.type(), p2.type());
           if (potentialAT._computeForce(force, p1, p2)) {
               force *= w12;
               real vir_temp = 0.5*dist[0]*force[0];
               p_xx_local.at(bin1) += vir_temp;
               p_xx_local.at(bin2) += vir_temp;
           }
         }
      }

      std::vector <real> p_xx_sum(size);
      for (i = 0; i < bins; ++i)
      {
          p_xx_sum.at(i) = 0.0;         
          boost::mpi::all_reduce(*mpiWorld, p_xx_local.at(i), p_xx_sum.at(i), std::plus<real>());
      }
      std::transform(p_xx_sum.begin(), p_xx_sum.end(), p_xx_sum.begin(),std::bind2nd(std::divides<real>(),Volume)); 
      for (i = 0; i < bins; ++i)
      {
          p_xx_total.at(i) += p_xx_sum.at(i);
      }
    
    }

    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the Verlet List");

      real w = 0.0;
      for (PairList::Iterator it(verletList->getPairs());                
           it.isValid(); ++it) {                                         
        Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.type();                                           
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          w = w + dist * force;
        }
      }

      for (PairList::Iterator it(verletList->getAdrPairs());
                 it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const PotentialCG &potential = getPotentialCG(type1, type2);

              Real3D force(0.0, 0.0, 0.0);
              if(potential._computeForce(force, p1, p2)) {
                Real3D dist = p1.position() - p2.position();
                w = w + dist * force;
              }
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      Tensor wlocal(0.0);
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it){
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D r21 = p1.position() - p2.position();
          wlocal += Tensor(r21, force);
        }
      }

      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it){
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const PotentialCG &potential = getPotentialCG(type1, type2);

              Real3D force(0.0, 0.0, 0.0);
              if(potential._computeForce(force, p1, p2)) {
                Real3D r21 = p1.position() - p2.position();
                wlocal += Tensor(r21, force);
              }
      }

      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }
    
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

    }
    
    
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListAdress does'n work"<<std::endl;
    }
    
    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListAdressInteractionTemplate< _PotentialAT, _PotentialCG >::
    getMaxCutoff() {
      real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotentialCG(i, j).getCutoff());
          cutoff = std::max(cutoff, getPotentialAT(i, j).getCutoff());
        }
      }
      return cutoff;
    }

  template < typename _PotentialAT, typename _PotentialCG >
  inline boost::python::list
  VerletListAdressInteractionTemplate< _PotentialAT, _PotentialCG >::getInteractionMatrix() {
    python::list params;

    for (int i = 0; i < ntypes; i++) {
      for (int j = 0; j < ntypes; j++) {
        params.append(python::make_tuple(i, j, getPotentialAT(i, j).getParams(), "AT"));
        params.append(python::make_tuple(i, j, getPotentialCG(i, j).getParams(), "CG"));
      }
    }
    return params;
  }

  }
}
#endif
