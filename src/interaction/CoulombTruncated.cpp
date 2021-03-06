/*
  Copyright (C) 2012,2013,2015
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
#include "CoulombTruncated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "FixedPairListAdressTypesInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHybridInteractionTemplate.hpp"
#include "Tabulated.hpp"

//For a Coulombic FixedPairList interaction, it's necessary to use FixedPairListTypesInteractionTemplate.hpp instead of FixedPairListInteractionTemplate.hpp
//so that we can use _computeForce and _computeEnergy which take both particles and distance vector as arguments
//because the Coulomb interaction needs access to both the charges (via the particles) and the minimum image distance (via the boundary conditions in the interaction template)

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< CoulombTruncated >
    VerletListCoulombTruncated;
    typedef class FixedPairListTypesInteractionTemplate< CoulombTruncated >
    FixedPairListTypesCoulombTruncated;

    typedef class FixedPairListAdressTypesInteractionTemplate<CoulombTruncated>
      FixedPairListAdressTypesCoulombTruncated;

    typedef class VerletListHybridInteractionTemplate<CoulombTruncated>
      VerletListHybridCoulombTruncated;

    typedef class VerletListAdressInteractionTemplate<CoulombTruncated, Tabulated>
        VerletListAdressCoulombTruncated;
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    CoulombTruncated::registerPython() {
      using namespace espressopp::python;

      class_< CoulombTruncated, bases< Potential > >
        ("interaction_CoulombTruncated", init< >())
        .def(init< real, real >())
        .add_property("prefactor", &CoulombTruncated::getPrefactor, &CoulombTruncated::setPrefactor)
      ;

      class_< VerletListCoulombTruncated, bases< Interaction > >
        ("interaction_VerletListCoulombTruncated", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListCoulombTruncated::setPotential)
          .def("getInteractionMatrix", &VerletListCoulombTruncated::getInteractionMatrix)
        .def("getPotential", &VerletListCoulombTruncated::getPotentialPtr)
        ;

      class_<VerletListHybridCoulombTruncated, bases<Interaction> >
          ("interaction_VerletListHybridCoulombTruncated", init<shared_ptr<VerletList>, bool>())
          .def("getVerletList", &VerletListHybridCoulombTruncated::getVerletList)
          .def("setPotential", &VerletListHybridCoulombTruncated::setPotential)
          .def("getPotential", &VerletListHybridCoulombTruncated::getPotentialPtr)
          .def("getInteractionMatrix", &VerletListHybridCoulombTruncated::getInteractionMatrix)
          .add_property("scale_factor", &VerletListHybridCoulombTruncated::scaleFactor,
                        &VerletListHybridCoulombTruncated::setScaleFactor)
          .add_property("max_force", &VerletListHybridCoulombTruncated::maxForce,
                        &VerletListHybridCoulombTruncated::setMaxForce)
          ;

      class_<VerletListAdressCoulombTruncated, bases<Interaction> >
          ("interaction_VerletListAdressCoulombTruncated",
                  init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
          .def("setPotentialAT", &VerletListAdressCoulombTruncated::setPotentialAT)
          .def("getInteractionMatrix", &VerletListAdressCoulombTruncated::getInteractionMatrix)
          .def("setPotentialCG", &VerletListAdressCoulombTruncated::setPotentialCG);
      ;

      class_< FixedPairListTypesCoulombTruncated, bases< Interaction > >
        ("interaction_FixedPairListTypesCoulombTruncated",
          init< shared_ptr<System>, shared_ptr<FixedPairList> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
        .def("setPotential", &FixedPairListTypesCoulombTruncated::setPotential)
          .def("getPotential", &FixedPairListTypesCoulombTruncated::getPotentialPtr)
        ;

      class_< FixedPairListAdressTypesCoulombTruncated, bases< Interaction > >
          ("interaction_FixedPairListAdressTypesCoulombTruncated",
           init< shared_ptr<System>, shared_ptr<FixedPairList>, bool >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, bool >())
          .def("setPotential", &FixedPairListAdressTypesCoulombTruncated::setPotential)
          .def("getPotential", &FixedPairListAdressTypesCoulombTruncated::getPotentialPtr)
          .def("setFixedPairList", &FixedPairListAdressTypesCoulombTruncated::setFixedPairList)
          .def("getFixedPairList", &FixedPairListAdressTypesCoulombTruncated::getFixedPairList)
          ;
    }

  }
}


