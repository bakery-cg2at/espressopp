"""
************************************
**Storage** - Storage Object
************************************

This is the abstract base class for all storage objects.
It cannot be used directly. All derived classes implement at least
the following methods:

* `decompose()`
* `addParticle(pid, pos)`
* `getParticle(pid)`
* `addAdrParticle(pid, pos, last_pos)`
* `setFixedTuples(fixed_tuple_list)`
* `addParticles(particle_list, *properties)`
* `modifyParticle(pid, property, value, decompose='yes')`

The property 'system' returns the System object of the storage.

Examples:
>>> s.storage.addParticles([[1, espresso.Real3D(3,3,3)], [2, espresso.Real3D(4,4,4)]],'id','pos')
>>> s.storage.decompose()

"""


from espresso import pmi
import MPI
import logging
from espresso import toReal3DFromVector, ParticleLocal, Particle
from espresso.Exceptions import ParticleDoesNotExistHere

class StorageLocal(object):
    """Abstract local base class for storing particles"""

    logger = logging.getLogger("Storage")

    def addParticle(self, pid, *args):
        """Add a particle locally if it belongs to the local domain."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.addParticle(
                self, pid, toReal3DFromVector(*args)
                )
            
    def addAdrATParticle(self, pid, *args):
        """Add a particle locally if it belongs to the local domain."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.addAdrATParticle(
                self, pid, toReal3DFromVector(*args)
                )
    
    def setFixedTuples(self, fixedtuples):
        if pmi.workerIsActive():
            self.cxxclass.setFixedTuples(self, fixedtuples)
    
    def getParticle(self, pid):
        """Get the local particle. If it is not on this node, any
        attempt to access the particle will raise an exception."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return ParticleLocal(pid, self)

    def addParticles(self, particleList, *properties):
        """
        This routine adds particles with certain properties to the storage
        where only one processor adds the particle in its local storage.

        :param particleList: list of list particles (and properties) to be added
        :param properties: property strings

        Each particle in the list must be itself a list where each entry corresponds
        to the property specified in properties.
        
        Example: addParticles([[id, pos, type, ... ], ...], 'id', 'pos', 'type', ...)

        Improvement: make list of supported properties more flexible in use
        """
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():

            index_id   = -1
            index_pos  = -1
            index_v    = -1
            index_f    = -1
            index_q    = -1
            index_type = -1
            index_mass = -1
            index_adrAT= -1 # adress AT particle if 1
            
            last_pos = toReal3DFromVector([-99,-99,-99])

            if properties == None:
                # default properities = (id, pos)
                index_id = 0
                index_pos = 1
                nindex = 2
            else:
                nindex = 0
                for val in properties:
                    if val.lower() == "id": index_id = nindex
                    elif val.lower() == "pos": index_pos = nindex
                    elif val.lower() == "type": index_type = nindex
                    elif val.lower() == "mass": index_mass = nindex
                    elif val.lower() == "v": index_v = nindex
                    elif val.lower() == "f": index_f = nindex
                    elif val.lower() == "q": index_q = nindex
                    elif val.lower() == "adrat": index_adrAT = nindex
                    else: raise "unknown particle property: %s"%val
                    nindex += 1

            if index_id < 0  : raise "particle property id is mandatory"
            if index_pos < 0 : raise "particle property pos is mandatory"

            for particle in particleList:

                # verify that each particle has enough entries, avoids index errors
                if len(particle) != nindex:
                    raise "particle has %d entries, but %d expected"%(len(particle), nindex)

                id = particle[index_id]
                pos = particle[index_pos]

                if index_adrAT >= 0:
                    if particle[index_adrAT] == 0:
                        #print "%d:  addParticle %d, last_pos=pos %d, %d, %d"%(pmi._MPIcomm.rank,id,pos[0], pos[1], pos[2])
                        storedParticle = self.cxxclass.addParticle(self, id, pos)
                        last_pos = pos
                    else:
                        #print "%d:  addAdrATparticle %d, last_pos %d, %d, %d"%(pmi._MPIcomm.rank,id,last_pos[0], last_pos[1], last_pos[2])
                        storedParticle = self.cxxclass.addAdrATParticle(self, id, pos, last_pos)
                else:
                    #print "%d:  addParticle %d, last_pos=pos %d, %d, %d"%(pmi._MPIcomm.rank,id,pos[0], pos[1], pos[2])
                    storedParticle = self.cxxclass.addParticle(self, id, pos)
                    last_pos = pos
                    
                if storedParticle != None:
                    self.logger.debug("Processor %d stores particle id = %d"%(pmi.rank, id))
                    self.logger.debug("particle property indexes: id=%i pos=%i type=%i mass=%i v=%i f=%i q=%i"%(index_id,index_pos,index_type,index_mass,index_v,index_f,index_q))

                    # only the owner processor writes other properties

                    if index_v >= 0:
                        storedParticle.v = particle[index_v]

                    if index_f >= 0:
                        storedParticle.f = particle[index_f] 

                    if index_q >= 0:
                        storedParticle.q = particle[index_q]

                    if index_type >= 0:
                        storedParticle.type = particle[index_type]

                    if index_mass >= 0:
                        storedParticle.mass = particle[index_mass]
 
    def modifyParticle(self, pid, property, value):
        """
        This routine allows to modify any properties of an already existing particle
        where only one processor modifies the particle in its local storage.
        
        Example: modifyParticle(pid, 'pos', Real3D(new_x, new_y, new_z))

        """
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            particle = self.getParticle(pid)
            if particle:
              try:
                if not particle.isGhost:
                  print "particle pid=%i rank=%i" % (pid, pmi.rank)
                  if   property.lower() == "id"   : raise "particles pid cannot be modified !"
                  elif property.lower() == "pos"  : particle.pos  = value
                  elif property.lower() == "type" : particle.type = value
                  elif property.lower() == "mass" : particle.mass = value
                  elif property.lower() == "v"    : particle.v    = value
                  elif property.lower() == "f"    : particle.f    = value
                  elif property.lower() == "q"    : particle.q    = value
                  else: raise UnknownParticleProperty( 'unknown particle property: %s' % property)
              except ParticleDoesNotExistHere:
                self.logger.debug("ParticleDoesNotExistHere pid=% rank=%i" % (pid, pmi.rank))
                pass
            
if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "decompose", "addParticles", "setFixedTuples", "modifyParticle" ],
            pmiproperty = [ "system" ]
            )

        def addParticle(self, pid, *args):
            pmi.call(self.pmiobject, 'addParticle', pid, *args)
            return Particle(pid, self)
        
        def addAdrATParticle(self, pid, *args):
            pmi.call(self.pmiobject, 'addAdrATParticle', pid, *args)
            return Particle(pid, self)
        
        #def setFixedTuples(self, tuples):
        #    pmi.call(self.pmiobject, 'setFixedTuples', tuples)

        def getParticle(self, pid):
            return Particle(pid, self)

