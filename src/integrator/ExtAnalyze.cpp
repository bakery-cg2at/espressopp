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

#include "python.hpp"
#include "types.hpp"
#include "ExtAnalyze.hpp"
#include "SystemAccess.hpp"
//#include "analysis/AnalysisBase.hpp"
#include "ParticleAccess.hpp"

namespace espressopp {
  //using namespace analysis;
  namespace integrator {

    LOG4ESPP_LOGGER(ExtAnalyze::theLogger, "ExtAnalyze");

    ExtAnalyze::ExtAnalyze(shared_ptr<ParticleAccess > particle_access, int interval):
        Extension(particle_access->getSystem()), interval_(interval){
      LOG4ESPP_INFO(theLogger, "Analyze observable in integrator");
      particle_access_ = particle_access;
      type = Extension::ExtAnalysis;
    }

    void ExtAnalyze::disconnect(){
      _aftIntV.disconnect();
    }

    void ExtAnalyze::connect(){
      // connection to end of integrator
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&ExtAnalyze::perform_action, this), boost::signals2::at_front);
    }

    //void ExtAnalyze::performMeasurement() {
    void ExtAnalyze::perform_action() {
      if ((integrator->getStep() % interval_) == 0) {
          LOG4ESPP_INFO(theLogger, "performing measurement in integrator " << integrator->getStep() << " " << interval_ << " " << integrator->getStep() % interval_);
          particle_access_->perform_action();
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/
    void ExtAnalyze::registerPython() {
      using namespace espressopp::python;
      class_<ExtAnalyze, shared_ptr<ExtAnalyze>, bases<Extension> >
        ("integrator_ExtAnalyze", init< shared_ptr<ParticleAccess>, int >())
        .def("connect", &ExtAnalyze::connect)
        .def("disconnect", &ExtAnalyze::disconnect)
        .add_property("interval", &ExtAnalyze::interval, &ExtAnalyze::set_interval)
        ;
    }
  }
}
