#  Copyright (c) 2015-2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
*******************************************
**espressopp.integrator.DynamicResolution**
*******************************************

*DynamicResolution* extension allows changing the *lambda* parameter of particles, namely
the so called resolution. The module can be used to perform backmapping.

The resolution is changed by

.. math::

  \lambda(t) = \lambda_0 + \alpha t

where :math:`\lambda` is the current resolution of the particles, :math:`\alpha` is
the rate by which the resolution is changed during the simulation and :math:`\lambda_0`
is an initial resolution.

.. function:: espressopp.integrator.DynamicResolution(system, vs_list, rate)



"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_DynamicResolution


class DynamicResolutionLocal(ExtensionLocal, integrator_DynamicResolution):
    'The (local) AdResS'

    def __init__(self, _system, _vs_list, _rate):
        'Local construction of a verlet list for AdResS'
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_DynamicResolution, _system, _vs_list, _rate)

    def update_weights(self):
        if pmi.workerIsActive():
            self.cxxclass.update_weights(self)

if pmi.isController:
    class DynamicResolution(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.DynamicResolutionLocal',
            pmiproperty = [ 'resolution', 'rate', 'active' ],
            pmicall = ['update_weights']
            )
