#ifndef _DOMAINDECOMPOSITION_HPP
#define _DOMAINDECOMPOSITION_HPP
#include "Storage.hpp"
#include "CellGrid.hpp"
#include "NodeGrid.hpp"

namespace espresso {

  class NodeGridMismatch: public std::runtime_error
  {
  public:
    NodeGridMismatch();
  };

  class DomainDecomposition: public Storage {
  public:
    DomainDecomposition(System *,
			const boost::mpi::communicator &,
			int nodeGrid[3],
			int cellGrid[3],
			bool useVList);

    virtual ~DomainDecomposition() {}

    void exchangeAndSortParticles();
    void sendGhostData();
    void collectGhostForces();

    virtual Cell *mapPositionToCellClipping(const real pos[3]);
    virtual Cell *mapPositionToCellChecked(const real pos[3]);

    const NodeGrid &getNodeGrid() const { return nodeGrid; }
    const CellGrid &getCellGrid() const { return cellGrid; }

  protected:
    /// init global Verlet list
    void initCellInteractions();
    /// set the grids and allocate space accordingly
    void createCellGrid(const int _nodeGrid[3], const int _cellGrid[3]);
    /// sort cells into local/ghost cell arrays
    void markCells();

    /// spatial domain decomposition of nodes
    NodeGrid nodeGrid;

    /// spatial domain decomposition on node in cells
    CellGrid cellGrid;

    /// expected capacity of send/recv buffers for neighbor communication
    size_t exchangeBufferSize;

    static LOG4ESPP_DECL_LOGGER(logger);
  };
}
#endif
