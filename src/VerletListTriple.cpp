#include "VerletListTriple.hpp"
#include "Real3D.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllTriplesIterator.hpp"

namespace espresso {

  using namespace espresso::iterator;

  LOG4ESPP_LOGGER(VerletListTriple::theLogger, "VerletListTriple");

/*-------------------------------------------------------------*/

  // cut is a cutoff (without skin)
  VerletListTriple::VerletListTriple(shared_ptr<System> system, real _cut, bool rebuildVL):SystemAccess(system){
    LOG4ESPP_INFO(theLogger, "construct VerletListTriple, cut = " << _cut);
  
    if (!system->storage) {
      throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletListTriple::rebuild, this));
  }
  
  real VerletListTriple::getVerletCutoff(){
    return cutVerlet;
  }
  
  void VerletListTriple::connect(){
    // make a connection to System to invoke rebuild on resort
    connectionResort = getSystem()->storage->onParticlesChanged.connect( 
            boost::bind(&VerletListTriple::rebuild, this));
  }

  void VerletListTriple::disconnect(){
    // disconnect from System to avoid rebuild on resort
    connectionResort.disconnect();
  }

  /*-------------------------------------------------------------*/
  
  void VerletListTriple::rebuild(){
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    
    vlTriples.clear();

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
    for (CellListAllTriplesIterator it(cl); it.isValid(); ++it) {
      checkTriple(*it->first, *it->second, *it->third);
      LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() <<
              ", " << it->second->id()<< " and "<<it->third->id());
    }

    builds++;
    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlTriples.size());
  }
  

  /*-------------------------------------------------------------*/
  
  void VerletListTriple::checkTriple(Particle& pt1, Particle& pt2, Particle& pt3){
    Real3D d1 = pt1.position() - pt2.position();
    Real3D d2 = pt2.position() - pt3.position();
    
    real distsq1 = d1.sqr();
    real distsq2 = d2.sqr();

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id()
                   << " @ " << pt1.position() 
		   << " - p2: " << pt2.id() << " @ " << pt2.position()
		   << " - p3: " << pt3.id() << " @ " << pt3.position()
		   << " -> distsq1 = " << distsq1
		   << " -> distsq2 = " << distsq2);

    if (distsq1>cutsq || distsq2>cutsq) return;

    // see if it's in the exclusion list (both directions)
    /*
    if (exList.count(boost::make_tuple(pt1.id(), pt2.id(), pt3.id())) == 1) return;
    if (exList.count(boost::make_tuple(pt3.id(), pt1.id(), pt2.id())) == 1) return;
    if (exList.count(boost::make_tuple(pt2.id(), pt3.id(), pt1.id())) == 1) return;
    
    if (exList.count(boost::make_tuple(pt2.id(), pt1.id(), pt3.id())) == 1) return;
    if (exList.count(boost::make_tuple(pt1.id(), pt3.id(), pt2.id())) == 1) return;
    if (exList.count(boost::make_tuple(pt3.id(), pt2.id(), pt1.id())) == 1) return;
    */
    
    vlTriples.add(pt1, pt2, pt3); // add triple to Verlet List
  }
  
  /*-------------------------------------------------------------*/
  
  int VerletListTriple::totalSize() const{
    System& system = getSystemRef();
    int size = localSize();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
  }

  int VerletListTriple::localSize() const{
    System& system = getSystemRef();
    return vlTriples.size();
  }

  python::tuple VerletListTriple::getTriple(int i) {
    if (i <= 0 || i > vlTriples.size()) {
      std::cout << "Warning! VerletList pair " << i << " does not exists" << std::endl;
      return python::make_tuple();
    }
    else{
      return python::make_tuple(vlTriples[i-1].first->id(),
              vlTriples[i-1].second->id(), vlTriples[i-1].third->id());
    }
  }


  bool VerletListTriple::exclude(longint pid1, longint pid2, longint pid3) {
    //exList.insert(boost::make_tuple(pid1, pid2, pid3));
    return true;
  }
  

  /*-------------------------------------------------------------*/
  
  VerletListTriple::~VerletListTriple()
  {
    LOG4ESPP_INFO(theLogger, "~VerletListTriple");
  
    if (!connectionResort.connected()){
      connectionResort.disconnect();
    }
  }
  
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VerletListTriple::registerPython() {
    using namespace espresso::python;

    bool (VerletListTriple::*pyExclude)(longint pid1, longint pid2, longint pid3)
          = &VerletListTriple::exclude;

    class_<VerletListTriple, shared_ptr<VerletListTriple> >
      ("VerletListTriple", init< shared_ptr<System>, real, bool >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletListTriple::getBuilds, &VerletListTriple::setBuilds)
      .def("totalSize", &VerletListTriple::totalSize)
      .def("localSize", &VerletListTriple::localSize)
      .def("getTriple", &VerletListTriple::getTriple)
      .def("exclude", pyExclude)
      .def("rebuild", &VerletListTriple::rebuild)
      .def("connect", &VerletListTriple::connect)
      .def("disconnect", &VerletListTriple::disconnect)
    
      .def("getVerletCutoff", &VerletListTriple::getVerletCutoff)
      ;
  }

}
