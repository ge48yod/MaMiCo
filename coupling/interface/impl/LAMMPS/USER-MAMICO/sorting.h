// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef LMP_MAMICO_SORTING_H_
#define LMP_MAMICO_SORTING_H_

#include "ghost_atoms.h"
#include "mamico_cell.h"

#include "coupling/CouplingMDDefinitions.h"
// MY NEW CHANGES , Maybe NOt needed
#include "coupling/indexing/IndexingService.h"
#include "coupling/indexing/IndexTypes.h"
#include "coupling/indexing/Operations.h"
#include "tarch/la/Vector.h"

#include "lammps/atom.h"
#include "lammps/lammps.h"
#include "lammps/memory.h"

namespace LAMMPS_NS {

/** this class is used for sorting particles into the MamicoCell structures and
 * to provide access to those.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class Sorting {
public:
  enum PrintType { PRINT_ALL_CELLS = 0, PRINT_INNER_CELLS = 1, PRINT_GHOST_CELLS = 2 };

  /** initialise "numberCells" mamico cells. This number must be big enough to
   * represent all local Mamico-coupling cells (incl.
   *  ghost cells) on every process.
   */
  Sorting(int numberCells, LAMMPS_NS::LAMMPS* lmp)
      : _lmp(lmp), _numberCells((unsigned int)numberCells), _mamicoCells(new MamicoCell[numberCells]), _ghostAtoms(lmp) {
    if (_mamicoCells == NULL) {
      std::cout << "ERROR Sorting: _mamicoCells==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  ~Sorting() {
    if (_mamicoCells != NULL)
      delete[] _mamicoCells;
    _numberCells = 0;
    _lmp = NULL;
  }

  /** prints the molecules in the mamico cells for debugging purposes */
  void printMolecules(LAMMPS_NS::Sorting<dim>::PrintType printType) {
    // MY BAD CHANGE
    // auto& idxs = IDXS;
    // tarch::la::Vector<3, unsigned int> loop(0);
    // const tarch::la::Vector<3, unsigned int> end =
    //     coupling::initRange<dim>(I11::numberCellsInDomain + tarch::la::Vector<dim, unsigned int>(2));

    // for (loop[2] = 0; loop[2] < end[2]; loop[2]++) {
    //   for (loop[1] = 0; loop[1] < end[1]; loop[1]++) {
    //     for (loop[0] = 0; loop[0] < end[0]; loop[0]++) {
    //       // decied whether to output this line or not
    //       bool decide = true;
    //       switch (printType) {
    //       case PRINT_ALL_CELLS:
    //         decide = true;
    //         break;
    //       case PRINT_INNER_CELLS:
    //         decide = true;
    //         for (unsigned int d = 0; d < dim; d++) {
    //           decide = decide && (loop[d] > 0) && (loop[d] < end[d] - 1);
    //         }
    //         break;
    //       case PRINT_GHOST_CELLS:
    //         decide = true;
    //         for (unsigned int d = 0; d < dim; d++) {
    //           decide = decide && (loop[d] > 0) && (loop[d] < end[d] - 1);
    //         }
    //         decide = !decide;
    //         break;
    //       default:
    //         std::cout << "ERROR printMolecules(): This case should never be reached!" << std::endl;
    //         exit(EXIT_FAILURE);
    //         break;
    //       }

    //       if (decide) {
    //         const unsigned int cellIndex = coupling::indexing::convertToScalar(I03{tarch::la::Vector<3, int>{coupling::initDimVector<dim>(loop)}});
    //         coupling::interface::MoleculeIterator<LAMMPS_NS::MamicoCell, dim>* it =
    //             coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getMoleculeIterator(
    //                 _mamicoCells[cellIndex]);

    //         for (it->begin(); it->continueIteration(); it->next()) {
    //           const coupling::interface::Molecule<dim>& molecule = it->getConst();
    //           std::cout << "Rank " << idxs.getRank() << ", cell " << coupling::initDimVector<dim>(loop) << ", molecule "
    //                     << molecule.getPosition() << std::endl;
    //         }
    //         delete it;
    //       }
    //     }


    for (auto idx : I02()) {
      // decied whether to output this line or not
      bool decide;
      switch (printType) {
      case PRINT_ALL_CELLS:
        decide = true;
        break;
      case PRINT_INNER_CELLS:
        decide = I10::contains(idx);
        break;
      case PRINT_GHOST_CELLS:
        decide = !I10::contains(idx);
        break;
      default:
        std::cout << "ERROR printMolecules(): This case should never be reached!" << std::endl;
        exit(EXIT_FAILURE);
        break;
      }

      if (decide) {
        coupling::interface::MoleculeIterator<LAMMPS_NS::MamicoCell, dim>* it =
            coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getMoleculeIterator(
                _mamicoCells[idx.get()]);

        for (it->begin(); it->continueIteration(); it->next()) {
          const coupling::interface::Molecule<dim>& molecule = it->getConst();
          std::cout << "Rank " << IDXS.getRank() << ", cell " << I03{idx} << ", molecule " << molecule.getPosition() << std::endl;
        }
        delete it;
      }
    }
  }

  /** returns the cell at index "index" */
  MamicoCell& getMamicoCell(unsigned int index) { return _mamicoCells[index]; }

  void updateAllCells() {
// remove all atoms from the cell lists
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Flag and reset cells..." << std::endl;
#endif
    flagAndResetCells();
// extract ghost atoms from lammps
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Extract ghost atoms..." << std::endl;
#endif
    _ghostAtoms.extractGhostAtoms();
// update all non-ghost cells
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Update non-ghost cells..." << std::endl;
#endif
    updateNonGhostCells(false);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    printMolecules(PRINT_INNER_CELLS);
#endif
// update all ghost cells
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Update ghost cells..." << std::endl;
#endif
    updateGhostCells();
  }

  /** removes all atoms from the non-ghost cells and sorts the "nlocal" atoms
   * from this process into the local cells.
   *  If "clearCellLists is set "false", the cell lists are not cleared before
   * the update is carried out; by default, the
   *  cells are emptied.
   */
  void updateNonGhostCells(bool clearCellLists = true) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    const tarch::la::Vector<dim, double> localOffset = getLocalOffset();
    const tarch::la::Vector<dim, double> localSize = getLocalSize();
    // reset all non-ghost cells

    if (clearCellLists) {
      //MY BAD CHANGE
      // const tarch::la::Vector<dim, unsigned int> localCells = I11::numberCellsInDomain + tarch::la::Vector<dim, unsigned int>(1);
      // tarch::la::Vector<3, unsigned int> end = coupling::initRange<dim>(localCells);
      // tarch::la::Vector<3, unsigned int> start(0);
      // for (unsigned int d = 0; d < dim; d++)
      //   start[d] = 1;
      // tarch::la::Vector<3, unsigned int> loop(0);
      // for (loop[2] = start[2]; loop[2] < end[2]; loop[2]++) {
      //   for (loop[1] = start[1]; loop[1] < end[1]; loop[1]++) {
      //     for (loop[0] = start[0]; loop[0] < end[0]; loop[0]++) {
      //       tarch::la::Vector<dim, unsigned int> index = coupling::initDimVector<dim>(loop);
      //       _mamicoCells[coupling::indexing::convertToScalar(I03{tarch::la::Vector<dim, int>(index)})].clear();
      //     }
      //   }


      for (auto idx : I10()) {
        _mamicoCells[idx.get()].clear();
      }
    }

    // loop over all local atoms and sort them into the cells
    const int nLocalAtoms = _lmp->atom->nlocal;
    for (int n = 0; n < nLocalAtoms; n++) {
      // extract atom position
      tarch::la::Vector<dim, double> position(0.0);
      for (unsigned int d = 0; d < dim; d++) {
        position[d] = _lmp->atom->x[n][d];
      }

      // determine global cell index for this atom; convert global to local
      // vector index
      //MY BAD CHANGE
//       if (isInLocalDomain(position, localOffset, localSize)) {
//         auto baseIdx = idxs.getCellIndex(position);
//         I01 globalVectorIndex{coupling::indexing::convertToVector<dim>(baseIdx)};
// #if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
//         std::cout << "Rank " << rank << ": Sort molecule at position " << position << " into global cell index " << globalVectorIndex.get() << std::endl;
//         // check if this is a global non-ghost cell and throw error otherwise
//         bool isInnerCell = true;
//         for (unsigned int d = 0; d < dim; d++) {
//           isInnerCell = isInnerCell && (globalVectorIndex.get()[d] > 0) && (globalVectorIndex.get()[d] < (int)I09::numberCellsInDomain[d] + 1);
//         }
//         if (!isInnerCell) {
//           std::cout << "ERROR Sorting::updateNonGhostCells: Molecule is not "
//                        "sorted into an inner cell!"
//                     << std::endl;
//           exit(EXIT_FAILURE);
//         }
// #endif
//         // convert global to local vector index
//         I03 localVectorIndex{globalVectorIndex.get()};
// // add atom to cell
// #if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
//         std::cout << "Sort molecule into local cell " << localVectorIndex.get() << ", corresponding to " << coupling::indexing::convertToScalar(localVectorIndex) << std::endl;
// #endif
//         _mamicoCells[coupling::indexing::convertToScalar(localVectorIndex)].addAtom(n);


      if (isInLocalDomain(position, localOffset, localSize)) {
        I00 idx = IDXS.getCellIndex(position);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
        std::cout << "Rank " << rank << ": Sort molecule at position " << position << " into global cell index " << I01{idx} << std::endl;
        // check if this is a global non-ghost cell and throw error otherwise
        if (!I04::contains(idx)) {
          std::cout << "ERROR Sorting::updateNonGhostCells: Molecule is not "
                       "sorted into an inner cell!"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
#endif
// add atom to cell
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
        std::cout << "Sort molecule into local cell " << I03{idx} << ", corresponding to " << I02{idx} << std::endl;
#endif
        _mamicoCells[I02{idx}.get()].addAtom(n);
      }
    }
  }

  /** forward call to _ghostAtoms */
  double** const getGhostAtomPositions() const { return _ghostAtoms.getGhostAtomPositions(); }

private:
  /** sets the ghost flag in all local mamico cells and removes all particle ids
   * from the cells */
  void flagAndResetCells() {
    //MY BAD CHANGE
//     // auto& idxs = IDXS;
//     const tarch::la::Vector<dim, unsigned int> localNumberCells = tarch::la::Vector<dim, unsigned int>(2) + I11::numberCellsInDomain;
//     tarch::la::Vector<3, unsigned int> loop(0);
//     tarch::la::Vector<3, unsigned int> end = coupling::initRange<dim>(localNumberCells);

//     // loop over all local mamico cells
//     unsigned int counter = 0;
//     for (loop[2] = 0; loop[2] < end[2]; loop[2]++) {
//       for (loop[1] = 0; loop[1] < end[1]; loop[1]++) {
//         for (loop[0] = 0; loop[0] < end[0]; loop[0]++) {
//           // determine ghost flag of cell
//           bool isInnerCell = true;
//           for (unsigned int d = 0; d < dim; d++) {
//             isInnerCell = isInnerCell && (loop[d] > 0) && (loop[d] < end[d] - 1);
//           }
// #if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
//           std::cout << "Flag cell " << loop << " to be ghost cell: " << isInnerCell << std::endl;
// #endif
//           _mamicoCells[counter].setGhostCell(!isInnerCell);
//           // reset vector with ids in cell
//           _mamicoCells[counter].clear();
//           counter++;
//         }
//       }

    for (auto idx : I02()) {
      bool isGhostCell = !I10::contains(idx);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "Flag cell " << I03{idx} << " to be ghost cell: " << isGhostCell << std::endl;
#endif
      _mamicoCells[idx.get()].setGhostCell(isGhostCell);
      // reset vector with ids in cell
      _mamicoCells[idx.get()].clear();
    }
  }

  /** sorts the ghost atoms from ghostX into the ghost cells */
  void updateGhostCells() {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    // compute bounding box for ghost particles -> if atoms are outside this
    // bounding box, we do not want to consider them anymore
    const tarch::la::Vector<dim, double> localOffset = getLocalOffset();
    const tarch::la::Vector<dim, double> localSize = getLocalSize();

    double** ghostX = _ghostAtoms.getGhostAtomPositions();
    const int nghost = _ghostAtoms.getNGhost();

    for (int n = 0; n < nghost; n++) {
      // extract atom position
      tarch::la::Vector<dim, double> position(0.0);
      for (unsigned int d = 0; d < dim; d++) {
        position[d] = ghostX[n][d];
      }

      // if the current position is in the region of interest, compute cell
      // index and add it to mamico cell
      if (isInLocalDomain(position, localOffset, localSize)) {
        //MY BAD CHANGE
//         // determine global cell index for this atom; convert global to local
//         // vector index
//         auto baseIdx = idxs.getCellIndex(position);
//         I01 globalVectorIndex{coupling::indexing::convertToVector<dim>(baseIdx)};
//         I03 localVectorIndex{globalVectorIndex.get()};
// #if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
//         std::cout << "Rank " << rank << ": Sort molecule at position " << position << " into global cell index "
//                   << globalVectorIndex.get() << std::endl;
//         // further check if this is a ghost cell and throw and error otherwise
//         bool isGhostCell = false;
//         for (unsigned int d = 0; d < dim; d++) {
//           isGhostCell = isGhostCell || (localVectorIndex.get()[d] == 0) || (localVectorIndex.get()[d] == (int)I11::numberCellsInDomain[d] + 1);
//         }
//         if (!isGhostCell) {
//           std::cout << "ERROR Sorting::updateGhostCells: Molecule is not "
//                        "sorted into a ghost cell!"
//                     << std::endl;
//           exit(EXIT_FAILURE);
//         }
// #endif
//         // add atom to cell
//         _mamicoCells[coupling::indexing::convertToScalar(localVectorIndex)].addAtom(n);


        // determine global cell index for this atom
        I02 idx = IDXS.getCellIndex(position);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
        std::cout << "Rank " << rank << ": Sort molecule at position " << position << " into global cell index " << I01{idx} << std::endl;
        // further check if this is a ghost cell and throw and error otherwise
        if (I11::contains(idx)) {
          std::cout << "ERROR Sorting::updateGhostCells: Molecule is not "
                       "sorted into a ghost cell!"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
#endif
        // add atom to cell
        _mamicoCells[idx.get()].addAtom(n);
      }
    }
  }

  tarch::la::Vector<dim, unsigned int> getThisProcess() const {
    // return {0}; // FIXME: Not yet implemented
    //MY NEW CHANGE
    return coupling::indexing::IndexingService<dim>::getInstance().getThisProcess();
  }
  tarch::la::Vector<dim, unsigned int> getAverageLocalNumberCouplingCells() const {
    // return {0}; // FIXME: Not yet implemented
    //MY NEW CHANGE
    return coupling::indexing::IndexingService<dim>::getInstance()
    .getAverageLocalNumberCouplingCells();
  }

  /** returns the offset of the local MD domain, incl. a ghost layer of mamico
   * cells */
  tarch::la::Vector<dim, double> getLocalOffset() const {
    // init local offset to global MD offset (very lower left of ghost layer)
    const tarch::la::Vector<dim, double> meshsize = IDXS.getCouplingCellSize();
    tarch::la::Vector<dim, double> localOffset = IDXS.getGlobalMDDomainOffset() - meshsize;

    // shift lower offset to the correct process
    const tarch::la::Vector<dim, unsigned int> thisProcess = getThisProcess();
    const tarch::la::Vector<dim, unsigned int> avgNumberCells = getAverageLocalNumberCouplingCells();
    for (unsigned int d = 0; d < dim; d++) {
      localOffset[d] += meshsize[d] * thisProcess[d] * avgNumberCells[d];
    }

    return localOffset;
  }

  /** returns the local domain size of this process, including a ghost layer of
   * mamico cells */
  tarch::la::Vector<dim, double> getLocalSize() const {
    const tarch::la::Vector<dim, double> meshsize = IDXS.getCouplingCellSize();
    const tarch::la::Vector<dim, unsigned int> localNumberCells = I11::numberCellsInDomain;
    tarch::la::Vector<dim, double> localSize(0.0);
    for (unsigned int d = 0; d < dim; d++) {
      localSize[d] = meshsize[d] * (localNumberCells[d] + 2);
    }
    return localSize;
  }

  /** returns true if the position vector lies inside the local domain,
   * specified by localOffset and localSize */
  bool isInLocalDomain(const tarch::la::Vector<dim, double>& position, const tarch::la::Vector<dim, double>& localOffset,
                       const tarch::la::Vector<dim, double>& localSize) const {
    bool isInside = true;
    const tarch::la::Vector<dim, double> upperOffset = localSize + localOffset;
    for (unsigned int d = 0; d < dim; d++) {
      isInside = isInside && (position[d] >= localOffset[d]) && (position[d] < upperOffset[d]);
    }
    return isInside;
  }

  LAMMPS_NS::LAMMPS* _lmp;     // points to our lammps instance
  unsigned int _numberCells;   // stores the number of allocated mamico cells
  MamicoCell* _mamicoCells;    // stores all mamico cells
  GhostAtoms<dim> _ghostAtoms; // stores and manages ghost atoms
};

} // namespace LAMMPS_NS

#endif // LMP_MAMICO_SORTING_H_
