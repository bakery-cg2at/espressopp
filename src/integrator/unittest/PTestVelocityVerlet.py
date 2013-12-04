#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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


import unittest
import espresso
import espresso.esutil
import espresso.unittest
import espresso.storage
import espresso.integrator
import espresso.interaction
import espresso.analysis
import espresso.bc
import MPI
import math
import logging

from espresso import Real3D

# Input values for system

N      = 10
cutoff = 2.5
skin   = 0.3

def calcNumberCells(size, nodes, cutoff):

    ncells = 1

    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1

    return ncells - 1

class TestVerletList(espresso.unittest.TestCase) :

    def test0Build(self) :

       system = espresso.System()

       rng  = espresso.esutil.RNG()

       SIZE = float(N)
       box  = Real3D(SIZE)
       bc   = espresso.bc.OrthorhombicBC(None, box)

       system.bc = bc
       system.rng = rng 
       system.skin = skin

       comm = espresso.MPI.COMM_WORLD

       nodeGrid = (1, 1, comm.size)
       cellGrid = [1, 1, 1]

       for i in range(3):
          cellGrid[i] = calcNumberCells(SIZE, nodeGrid[i], cutoff)

       print 'NodeGrid = %s'%(nodeGrid,)
       print 'CellGrid = %s'%cellGrid

       dd = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)

       system.storage = dd

       id = 0

       for i in range(N):
         for j in range(N):
           for k in range(N):
 
             m = (i + 2*j + 3*k) % 11
             r = 0.45 + m * 0.01
             x = (i + r) / N * SIZE
             y = (j + r) / N * SIZE
             z = (k + r) / N * SIZE
   
             dd.addParticle(id, Real3D(x, y, z))

             # not yet: dd.setVelocity(id, (1.0, 0.0, 0.0))

             id = id + 1

       dd.decompose()

       integrator = espresso.integrator.VelocityVerlet(system)

       print 'integrator.dt = %g, will be set to 0.005'%integrator.dt

       integrator.dt = 0.005
  
       print 'integrator.dt = %g, is now '%integrator.dt

       # now build Verlet List
       # ATTENTION: you have to add the skin explicitly here

       vl = espresso.VerletList(system, cutoff = cutoff + system.skin)

       potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = cutoff)

       # ATTENTION: auto shift was enabled

       print "potLJ, shift = %g"%potLJ.shift

       interLJ = espresso.interaction.VerletListLennardJones(vl)

       interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)

       # Todo

       system.addInteraction(interLJ)

       temp = espresso.analysis.Temperature(system)

       temperature = temp.compute()
       kineticEnergy = 0.5 * temperature * (3 * N * N * N)
       potentialEnergy = interLJ.computeEnergy()
       print 'Start: tot energy = %10.6f pot = %10.6f kin = %10.f temp = %10.6f'%(kineticEnergy + potentialEnergy,
                  potentialEnergy, kineticEnergy, temperature)

       nsteps = 10

       # logging.getLogger("MDIntegrator").setLevel(logging.DEBUG)

       for i in range(20):
          integrator.run(nsteps)
          temperature = temp.compute()
          kineticEnergy = 0.5 * temperature * (3 * N * N * N)
          potentialEnergy = interLJ.computeEnergy()
          print 'Step %6d: tot energy = %10.6f pot = %10.6f kin = %10.6f temp = %f'%(nsteps*(i+1), 
               kineticEnergy + potentialEnergy, potentialEnergy, kineticEnergy, temperature)

if __name__ == "__main__":

    unittest.main()
