// ESPP_CLASS
#ifndef _INTERACTION_FIXEDQUADRUPLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDQUADRUPLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "FixedQuadrupleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace interaction {
    template < typename _DihedralPotential >
    class FixedQuadrupleListInteractionTemplate : public Interaction, SystemAccess {
    protected:
      typedef _DihedralPotential Potential;
    public:
      FixedQuadrupleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedQuadrupleList > _fixedquadrupleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), fixedquadrupleList(_fixedquadrupleList),
          potential(_potential)
      {
          if (! potential) {
              LOG4ESPP_ERROR(theLogger, "NULL potential");
          }

        //potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      void
      setFixedQuadrupleList(shared_ptr < FixedQuadrupleList > _fixedquadrupleList) {
        fixedquadrupleList = _fixedquadrupleList;
      }

      virtual ~FixedQuadrupleListInteractionTemplate() {};

      shared_ptr < FixedQuadrupleList > getFixedQuadrupleList() {
             return fixedquadrupleList;
      }

      void
      setPotential(shared_ptr < Potential> _potential) {
           if (_potential) {
              potential = _potential;
           } else {
              LOG4ESPP_ERROR(theLogger, "NULL potential");
           }
      }


      /*void
      setPotential(int type1, int type2, const Potential &potential) {
	potentialArray.at(type1, type2) = potential;
      }*/

      /*Potential &getPotential(int type1, int type2) {
        return potentialArray(0, 0);
      }*/

      shared_ptr < Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real xmin, real xmax,
          real ymin, real ymax, real zmin, real zmax);
      virtual real getMaxCutoff();
      virtual int bondType() { return Dihedral; }

    protected:
      int ntypes;
      shared_ptr < FixedQuadrupleList > fixedquadrupleList;
      //esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _DihedralPotential > inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    addForces() {

      LOG4ESPP_INFO(theLogger, "add forces computed by FixedQuadrupleList");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        Particle &p4 = *it->fourth;
        //const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D dist21, dist32, dist43; // 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

	    Real3D force1, force2, force3, force4;  // result forces

	    potential->_computeForce(force1, force2, force3, force4,
                                dist21, dist32, dist43);
        p1.force() += force1;
        p2.force() -= force2;
        p3.force() += force3;
        p4.force() += force4;
      }
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the quadruples");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      real e = 0.0;
      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;
        //const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D dist21, dist32, dist43; // 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

        e += potential->_computeEnergy(dist21, dist32, dist43);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the quadruples");

      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;
        //const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D dist21, dist32, dist43; 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        potential->_computeForce(force1, force2, force3, force4,
                                dist21, dist32, dist43);

        // TODO: formulas are not correct yet?

        w += dist21 * force1 + dist32 * force2;
      }
      
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return w;
    }

    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");
    
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;

      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;
        //const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D dist21, dist32, dist43; 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        potential->_computeForce(force1, force2, force3, force4,
                                dist21, dist32, dist43);

        // TODO: formulas are not correct yet

        wlocal += Tensor(dist21, force1) - Tensor(dist32, force2);
      }
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }


    // compute the pressure tensor localized between xmin, xmax, ymin, ymax, zmin, zmax
    // TODO !!!!! This doesn't work
    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor& w,
            real xmin, real xmax, real ymin, real ymax, real zmin, real zmax) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");
    
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      
      std::cout<<"Warning!!! computeVirialTensor in specified volume doesn't work for "
              "FixedQuadrupleListInteractionTemplate at the moment"<<std::endl;

      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        Real3D dist21, dist32, dist43; 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        potential->_computeForce(force1, force2, force3, force4,
                                dist21, dist32, dist43);

        // TODO: formulas are not correct yet

        wlocal += Tensor(dist21, force1) - Tensor(dist32, force2);
      }
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleListInteractionTemplate< _DihedralPotential >::
    getMaxCutoff() {
      /*real cutoff = 0.0;

      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }*/
      return potential->getCutoff();
    }
  }
}
#endif
