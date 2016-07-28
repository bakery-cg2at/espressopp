/*
  Copyright (C) 2015, 2016
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
#include "DihedralRB.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"
#include "FixedQuadrupleListTypesInteractionTemplate.hpp"
#include "FixedQuadrupleListAdressInteractionTemplate.hpp"

namespace espressopp {
namespace interaction {
void DihedralRB::registerPython() {
  using namespace espressopp::python;  //NOLINT

  class_<DihedralRB, bases<DihedralPotential> >(
      "interaction_DihedralRB", init< real, real, real, real, real, real >())
      .def("getParams", &DihedralRB::getParams);

  typedef class FixedQuadrupleListInteractionTemplate<DihedralRB> FixedQuadrupleListDihedralRB;
  class_ <FixedQuadrupleListDihedralRB, bases <Interaction> >
      ("interaction_FixedQuadrupleListDihedralRB",
          init<shared_ptr<System>, shared_ptr<FixedQuadrupleList>, shared_ptr<DihedralRB> >())
      .def(init<shared_ptr<System>, shared_ptr<FixedQuadrupleListAdress>,
                shared_ptr<DihedralRB> >())
      .def("setPotential", &FixedQuadrupleListDihedralRB::setPotential)
      .def("getPotential", &FixedQuadrupleListDihedralRB::getPotential)
      .def("getFixedQuadrupleList", &FixedQuadrupleListDihedralRB::getFixedQuadrupleList);

  typedef class FixedQuadrupleListTypesInteractionTemplate<DihedralRB>
    FixedQuadrupleListTypesDihedralRB;
  class_< FixedQuadrupleListTypesDihedralRB, bases< Interaction > >
    ("interaction_FixedQuadrupleListTypesDihedralRB",
       init< shared_ptr<System>, shared_ptr<FixedQuadrupleList> >())
      .def("setPotential", &FixedQuadrupleListTypesDihedralRB::setPotential)
      .def("getPotential", &FixedQuadrupleListTypesDihedralRB::getPotentialPtr)
      .def("setFixedQuadrupleList", &FixedQuadrupleListTypesDihedralRB::setFixedQuadrupleList)
      .def("getFixedQuadrupleList", &FixedQuadrupleListTypesDihedralRB::getFixedQuadrupleList);

  typedef class FixedQuadrupleListAdressInteractionTemplate <DihedralRB> FixedQuadrupleListAdressDihedralRB;
  class_ <FixedQuadrupleListAdressDihedralRB, bases <Interaction> >
  ("interaction_FixedQuadrupleListAdressDihedralRB",
	  init< shared_ptr<System>,
	        shared_ptr<FixedQuadrupleList>,
	        shared_ptr<DihedralRB>, bool >())
  .def(init< shared_ptr<System>, shared_ptr<FixedQuadrupleListAdress>,
      shared_ptr<DihedralRB>, bool>())
  .def("setPotential", &FixedQuadrupleListAdressDihedralRB::setPotential)
      .def("getPotential", &FixedQuadrupleListAdressDihedralRB::getPotential)
      .def("getFixedQuadrupleList", &FixedQuadrupleListAdressDihedralRB::getFixedQuadrupleList);


}

}  // end namespace interaction
}  // end namespace espressopp
