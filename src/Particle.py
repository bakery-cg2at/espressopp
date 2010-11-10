import _espresso
import esutil
import pmi
from espresso import toReal3DFromVector, MPI

# Controller Particle:
# * requests are directly forwarded

# Parallel Particle:
# * throw exception if particle does not exist locally
# * otherwise do it

# Parallel Ghost Particle:
# * throw exception if particle does not exist locally
# * can not be written
# * should throw exception if data is not available

# The class _TmpParticle wraps the C++ internal pointer to a particle
# _TmpParticle should not be used, as it might die easily and will
# cause a SegFault when used after it has died.

class ParticleLocal(object):
    """The local particle.

    Throws an exception:
    * when the particle does not exists locally

    TODO: Should throw an exception:
    * when a ghost particle is to be written 
    * when data is to be read from a ghost that is not available
    """
    def __init__(self, pid, storage):
        self.pid = pid
        self.storage = storage

    def __getTmp(self):
        tmp = self.storage.lookupLocalParticle(self.pid)
        if tmp is None:
            # TODO: Exception
            raise 'Particle not here!'
        else:
            return tmp

    # Defining __getattr__ will make sure that you can use any
    # property defined in _TmpParticle
    def __getattr__(self, key):
        return getattr(self.__getTmp(), key)

#     def __setattr__(self, key, value):
#         return setattr(self.__getTmp(), key, value)

    # The following properties are modified between Python and C++
    @property
    def f(self): return self.__getTmp().f
    @f.setter
    def f(self, val): self.__getTmp().f = toReal3DFromVector(val)

    @property
    def v(self): return self.__getTmp().v
    @v.setter
    def v(self, val): self.__getTmp().v = toReal3DFromVector(val)

    @property
    def pos(self): return self.__getTmp().pos
    @pos.setter
    def pos(self, val): self.__getTmp().pos = toReal3DFromVector(val)

    def getLocalData(self, key):
        tmp = self.storage.lookupRealParticle(self.pid)
        if tmp is not None:
            return getattr(tmp, key)
        else:
            return None

    def locateParticle(self):
        tmp = self.storage.lookupRealParticle(self.pid)
        return (tmp is not None)

if pmi.isController:
    class Particle(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.ParticleLocal',
            pmiproperty = [ "id", "storage" ]
            )

        @property
        def node(self):
            value, node = pmi.reduce(pmi.MAXLOC, self, 'locateParticle')
            return node

        def __getattr__(self, key):
            value, node = pmi.reduce(pmi.MAXLOC, self, 'getLocalData', key)
            return value
