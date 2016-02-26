/*
  Copyright (C) 2015
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
#include "Resolution.hpp"
#include "storage/DomainDecomposition.hpp"
#include "integrator/DynamicResolution.hpp"

using namespace espressopp;  //NOLINT

namespace espressopp {
namespace analysis {

real Resolution::compute_real() const {
  return res_->resolution();
}

void Resolution::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<Resolution, bases<Observable> >
    ("analysis_Resolution", init< shared_ptr<System>,
                            shared_ptr<integrator::DynamicResolution> >())
    .add_property("value", &Resolution::compute_real);
}
}  // end namespace analysis
}  // end namespace espressopp
