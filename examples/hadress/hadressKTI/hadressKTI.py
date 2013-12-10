#!/usr/bin/env python
# -*- coding: utf-8 -*-

# relevant imports
import sys
import time
import espresso
import MPI

import Tetracryst # Preparation of tetrahedral crystal and constuctions of bonds in tetrahedral liquid

from espresso import Real3D, Int3D
from espresso.tools import decomp
from espresso.tools import timers

# timestep, cutoff, skin, AdResS specifications
timestep = 0.0005

rc = 4.5 # cutoff coarse-grained potential
rca = 1.122462048309373 # cutoff atomistic potential (cutoff (2^(1/6)), WCA)
skin = 0.4

# Parameters for the thermostat
gamma = 2.0
temp = 1.0

# Parameters for size of AdResS dimensions
ex_size = 500.0  # By choosing some random but large value here we make sure, that we have an "atomistic" region in the whole box.
                 # Although we do not perform an actual H-AdResS simulation we need the H-AdResS algorithms and the force calculation
                 # as it is performed in the atomistic region of an H-AdResS simulation.
hy_size = 10.0

# prepare tetrahedral liquid in crystal 
pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espresso.tools.readxyz("equilibrated_confKTI.xyz")

# Table for coarse-grained potential
tabCG = "table_potential.dat"

# number of CG particles
num_particlesCG = len(x)/4

# number of AT particles
num_particles = len(x)

# set up the system
sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)

# (H-)AdResS domain decomposition
system.storage = espresso.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

# prepare AT particles
allParticlesAT = []
allParticles = []
tuples = []
for pidAT in range(num_particles):
    allParticlesAT.append([pidAT, # add here these particles just temporarily! 
                         Real3D(x[pidAT], y[pidAT], z[pidAT]), # position
                         Real3D(vx[pidAT], vy[pidAT], vz[pidAT]), # velocity
                         Real3D(0, 0, 0),
                         1, 1.0, 1]) # type, mass, is AT particle

# create CG particles
for pidCG in range(num_particlesCG):
    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    cmp = espresso.tools.AdressSetCG(4, pidCG, allParticlesAT)
    # Preparation of tuples (tuples define, which atoms belong to which CG molecules)
    tmptuple = [pidCG+num_particles]
    for pidAT2 in range(4):
        pid = pidCG*4+pidAT2
        tmptuple.append(pid)
    
    # append CG particles    
    allParticles.append([pidCG+num_particles, # CG particle has to be added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         0, 4.0, 0]) # type, mass, is not AT particle
    # append AT particles
    for pidAT in range(4): 
        pid = pidCG*4+pidAT
        allParticles.append([pid, # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6]]) # is AT particle 
    # append tuple to tuplelist    
    tuples.append(tmptuple)
    

# add particles to system
system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

# add tuples to system
ftpl = espresso.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

# add bonds between AT particles
fpl = espresso.FixedPairListAdress(system.storage, ftpl)
bonds = Tetracryst.makebonds(len(x))
fpl.addBonds(bonds)

# decompose after adding tuples and bonds
print "Added tuples and bonds, decomposing now ..." 
system.storage.decompose()
print "done decomposing"

# AdResS Verlet list
vl = espresso.VerletListAdress(system, cutoff=rc, adrcut=rc,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2])

# non-bonded potentials
# LJ Capped WCA between AT and tabulated potential between CG particles
interNB = espresso.interaction.VerletListHadressLennardJones(vl, ftpl, KTI = True)  # Switch on KTI here!
potWCA  = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto', cutoff=rca)
potCG = espresso.interaction.Tabulated(itype=3, filename=tabCG, cutoff=rc) # CG
interNB.setPotentialAT(type1=1, type2=1, potential=potWCA) # AT
interNB.setPotentialCG(type1=0, type2=0, potential=potCG) # CG
system.addInteraction(interNB)

# bonded potentials
# Quartic potential between AT particles
potQuartic = espresso.interaction.Quartic(K=75.0, r0=1.0)
interQuartic = espresso.interaction.FixedPairListQuartic(system, fpl, potQuartic)
system.addInteraction(interQuartic)

# VelocityVerlet integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

# add AdResS extension
adress = espresso.integrator.Adress(system, vl, ftpl, True)
integrator.addExtension(adress)

# add Langevin thermostat extension
langevin = espresso.integrator.LangevinThermostat(system)
langevin.gamma = gamma
langevin.temperature = temp
langevin.adress = True # enable AdResS!
integrator.addExtension(langevin)

# distribute atoms and CGmolecules according to AdResS domain decomposition 
espresso.tools.AdressDecomp(system, integrator)

# system information
print ''
print 'number of AT particles =', num_particles
print 'number of CG particles =', num_particlesCG
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)

# Timer, Steps
start_time = time.clock()
  
# Set lambdas and derivates to zero   
for i in range(num_particles + num_particlesCG):
 system.storage.modifyParticle(i, 'lambda_adrd', 0.0)
 system.storage.modifyParticle(i, 'lambda_adr', 0.0)

### EQUILIBRATION ###
# Equilibration parameters
EQsteps = 1000
EQintervals = 100
EQnsteps = EQsteps/EQintervals

print ''
print 'Short equilibration'
print 'Equilibration steps =', EQsteps
print ''

# Print the data of the intial configuration
fmt = '%5d %8.4f %10.5f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n'
T = temperature.compute()
P = pressure.compute()
Ek = 0.5 * T * (3.0 * num_particles)
Ep = interNB.computeEnergy()
Eb = interQuartic.computeEnergy()
Eaa = interNB.computeEnergyAA()
Ecg = interNB.computeEnergyCG()
sys.stdout.write(' step      T         P         etotal     enonbonded    ebonded      ekinetic     eallatom        ecg \n')
sys.stdout.write(fmt % (0, T, P, Ek + Ep + Eb, Ep, Eb, Ek, Eaa, Ecg))

# Do equilibration
for s in range(1, EQintervals + 1):
  integrator.run(EQnsteps)
  EQstep = EQnsteps * s
  T = temperature.compute()
  P = pressure.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interNB.computeEnergy()
  Eb = interQuartic.computeEnergy()
  Eaa = interNB.computeEnergyAA()
  Ecg = interNB.computeEnergyCG()
  sys.stdout.write(fmt % (EQstep, T, P, Ek + Ep + Eb, Ep, Eb, Ek, Eaa, Ecg))

print ''
print 'Equilibration Done'
print ''

### KIRKWOOD TI ###
# TI parameters
bins = 1000
steps = 100
stepsequi = 50
intervals = 10

nstepsTI = steps/intervals  
lambdastep = 1.0/bins 

# Spefify output filename
namerawFile = 'KirkwoodTI_rawdata.dat'

print ''
print 'Starting Kirkwood TI'
print ''
print 'Kirkwood TI steps =', bins
print 'Kirkwood TI stepwidth =', lambdastep
print 'Integration steps for each lambda =', steps
print 'Equilibration steps after each lamda switch =', stepsequi
print 'Intervals for taking data and printing information to screen =', intervals
print ''

# Print the data of the starting configuration
fmt = '%5d %8.4f %10.5f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n'
T = temperature.compute()
P = pressure.compute()
Ek = 0.5 * T * (3.0 * num_particles)
Ep = interNB.computeEnergy()
Eb = interQuartic.computeEnergy()
Eaa = interNB.computeEnergyAA()
Ecg = interNB.computeEnergyCG()
sys.stdout.write(' step      T         P         etotal     enonbonded    ebonded      ekinetic     eallatom        ecg \n')
sys.stdout.write(fmt % (0, T, P, Ek + Ep + Eb, Ep, Eb, Ek, Eaa, Ecg))
print ''

#Output arrays
Energydiff = []
Pressurediff = []

# Kirkwood steps
for i in range(bins+1):
   
 # Changing Lambda
 print 'Kirkwood step: %d' %i
 print 'Lambda: %f' %(lambdastep*i)
 for p in range(num_particles + num_particlesCG):
   system.storage.modifyParticle(p, 'lambda_adr', lambdastep*i)

 # Equilibration
 integrator.run(stepsequi)
 step = i * (steps+stepsequi) + stepsequi
 T = temperature.compute()
 P = pressure.compute()
 Ek = 0.5 * T * (3.0 * num_particles)
 Ep = interNB.computeEnergy()
 Eb = interQuartic.computeEnergy()
 Eaa = interNB.computeEnergyAA()
 Ecg = interNB.computeEnergyCG()
 sys.stdout.write(fmt % (step, T, P, Ek + Ep + Eb, Ep, Eb, Ek, Eaa, Ecg))
   
 # Kirkwood Integration
 runningEdiff = 0.0
 runningP = 0.0   
 for s in range(1,intervals+1):
   integrator.run(nstepsTI)
   step = i * (steps+stepsequi) + s * nstepsTI + stepsequi
   T = temperature.compute()
   P = pressure.compute()
   Ek = 0.5 * T * (3.0 * num_particles)
   Ep = interNB.computeEnergy()
   Eb = interQuartic.computeEnergy()
   Eaa = interNB.computeEnergyAA()
   Ecg = interNB.computeEnergyCG()
   sys.stdout.write(fmt % (step, T, P, Ek + Ep + Eb, Ep, Eb, Ek, Eaa, Ecg))
     
   # Get the relevant energy and pressure differences
   runningEdiff += Ecg - Eaa
   runningP += P

 # Get the averages
 runningEdiff/=intervals
 runningP/=intervals
   
 # Append to output arrays
 Energydiff.append(runningEdiff) 
 Pressurediff.append(runningP) 

# Printing the raw output to file
print ''
print "Kirkwood TI done, printing raw data to %s\n" %namerawFile
form = '%12.8f %12.8f %12.8f\n'
rawFile = open (namerawFile, 'w')
rawFile.write('lambda  V_CG-V_AA  P(lambda)\n')
for i in range( bins+1 ):
 rawFile.write(form % ( lambdastep*i, Energydiff[i], Pressurediff[i] ))
rawFile.close()

# simulation information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

