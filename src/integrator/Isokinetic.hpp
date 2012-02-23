// ESPP_CLASS
#ifndef _INTEGRATOR_ISOKINETIC_HPP
#define _INTEGRATOR_ISOKINETIC_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace integrator {

    /** Langevin */

    class Isokinetic : public SystemAccess {

      public:

        Isokinetic(shared_ptr< System > system);

        void setTemperature(real temperature);

        real getTemperature();

        void setCoupling(int coupling);

        int getCoupling();

        ~Isokinetic();

        void rescaleVelocities();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        real temperature;  //!< desired user temperature
        int coupling; // how often to couple to the thermostat
        int couplecount;

        // not yet needed
        // shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
