/*
  Copyright (c) 2015-2016
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
#ifndef _DYNAMICRESOLUTION_HPP
#define _DYNAMICRESOLUTION_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "FixedVSList.hpp"
#include "Extension.hpp"

#include "boost/signals2.hpp"


namespace espressopp {
namespace integrator {

/*
 * This module implement dynamic resolution extension.
 */
class DynamicResolution : public Extension {
 public:
  DynamicResolution(shared_ptr<System> _system, shared_ptr<FixedVSList> _vslist, real _rate);

  ~DynamicResolution();
  real resolution() { return resolution_; }
  void set_resolution(real resolution) {
    resolution_ = resolution;
  }

  real rate() { return rate_; }
  void set_rate(real rate) { rate_ = rate; }

  bool active() { return active_; }
  void set_active(bool active);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  void updateWeights();

  void changeResolution();
  real rate_;
  bool active_;

  real resolution_;  /// Current value of resolution, between 0.0 - 1.0
  void connect();
  void disconnect();

  shared_ptr<FixedVSList> vs_list;

  // Signals
  boost::signals2::connection _aftIntV, _runInit;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};
}  // namespace integrator
}  // namespace espressopp
#endif
