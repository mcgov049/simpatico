#ifndef DDMD_CELL_LIST_H
#define DDMD_CELL_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Cell.h"
#include <ddMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/space/Grid.h>
#include <util/containers/DArray.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /**
   * A cell list used only to identify nearby atom pairs.
   *
   * An CellList divides the domain owned by a processor, plus a frame
   * containing ghost particles, into a grid of cells, such that the length
   * of each cell in each of Cartesian direction is greater than a specified 
   * cutoff distance. 
   *
   * All operations of this class are local (no MPI).
   *
   * Building a CellList (usage):
   * \code
   *
   *    AtomStorage storage;
   *    CellList cellList;
   *    Vector   lower;        // Vector of lower bounds (local atoms)
   *    Vector   upper;        // Vector of upper bounds (local atoms)
   *    double   cutoff;       // minimum cell dimension
   *    int      atomCapacity  // max number of atoms on this processor
   *
   *    // Bounds on lower and upper used here to allocate memory.
   *    cellList.allocate(atomCapacity, lower, upper, cutoff);
   *  
   *    // Make the actual grid.
   *    cellList.makeGrid(lower, upper, cutoff);
   *
   *    // Place all atoms and ghosts
   *    cellList.clear();
   *    AtomStorage::AtomIterator  atomIter;
   *    for (storage.begin(atomIter); atomIter.notEnd(); ++atomIter) {
   *       cellList.placeAtom(*atomIter);
   *    }
   *    AtomStorage::GhostIterator ghostIter;
   *    for (storage.begin(ghostIter); ghostIter.notEnd(); ++ghostIter){
   *       cellList.placeAtom(*ghostIter);
   *    }
   *    
   *    // Build cell list
   *    cellList.build();
   *
   *
   * \endcode
   *
   * The atomCapacity parameter should be set equal to the sum of atomCapacity 
   * and the ghostCapacity of the associated atomStorage, which is the maximum 
   * total number of atoms that can exist on this processor.
   *
   * See Cell documentation for an example of how to iterate over local cells 
   * and neighboring atom pairs. 
   * 
   * \ingroup DdMd_Neighbor_Module
   */
   class CellList
   {

   public:

      // Public methods

      /**
      * Constructor.
      */
      CellList();

      /**
      * Destructor.
      *
      * Deallocates array of Cell objects, if necessary.
      */
      virtual ~CellList();

      /**
      * Allocate memory for this CellList.
      *
      * This function:
      *
      *   - Allocates an array of atomCapacity Tag objects.
      *   - Allocates an array of atomCapacity Atom* objects.
      *   - Allocates an array of Cell objects sized for this boundary.
      *
      * \param atomCapacity dimension of global array of atoms
      * \param lower        lower bound used to allocate array of cells.
      * \param upper        upper bound used to allocate array of cells.
      * \param cutoff       minimum dimension of a cell in any direction
      */
      void allocate(int atomCapacity, const Vector& lower, const Vector& upper, 
                    double cutoff);

      /**
      * Make the cell grid.
      *
      * The number of cells in each direction is chosen such that the dimension
      * of each cell in each direction is greater than or equal to the cutoff
      * parameter. To calculate nonbonded pair interaction energies, the cutoff
      * parameter should thus be equal to or greater than the maximum range of
      * nonbonded interactions.
      *
      * \param lower    lower bound of local atom coordinates.
      * \param upper    upper bound of local atom coordinates.
      * \param cutoff   minimum dimension of a cell in any direction
      */
      void makeGrid(const Vector& lower, const Vector& upper, double cutoff);

      /**
      * Determine the appropriate cell for an Atom, based on its position.
      *
      * This method does not place the atom in a cell, but retains a record
      * of the cell index that is used to place atoms in the build() method.
      *
      * This method quietly does nothing if the atom is outside the expanded 
      * domain for nonbonded ghosts, which extends one cutoff length beyond
      * the domain boundaries the domain boundaries in each direction.
      *
      * \param atom  Atom object to be added.
      */
      void placeAtom(Atom &atom);

      /**
      * Build the cell list.
      *
      * This method must be called after completing a loop over atoms,
      * in which placeAtom() is called for each atom. The build() method
      * uses information gathered in this loop to build and fill all of
      * the cells. 
      */
      void build();

      /**
      * Reset the cell list to its empty state (no Atoms).
      */
      void clear();

      /**
      * Return Grid object by const reference.
      */
      const Grid& grid() const;

      /**
      * Dimension of each cell in direction i.
      *
      * \param i Cartesian index i = 0, 1, 2
      */
      double cellLength(int i) const;

      /**
      * Return the index of the cell that contains a position Vector. 
      *
      * Returns a null value of -1 if the position is outside the 
      * expanded domain for nonbonded ghost atoms. 
      *
      * \param position position Vector, inside the boundary
      */
      int cellIndexFromPosition(const Vector& position) const;

      /**
      * Return pointer to first local cell in linked list.
      */
      const Cell* begin() const;

      /**
      * Return a specified cell by const reference.
      * 
      * \param i cell index
      */
      const Cell& cell(int i) const;

      /**
      * Get total number of atoms (local and ghost) in this CellList.
      */
      int nAtom() const;

      /**
      * Get number of atoms that were rejected (not placed in cells)
      */
      int nReject() const;

      /**
      * Maximum number of atoms for which space is allocated.
      */
      int atomCapacity() const;

      /**
      * Number of cells for which space has been allocated.
      */
      int cellCapacity() const;

      /**
      * Has memory been allocated for this CellList?
      */
      bool isAllocated() const;

      /**
      * Has this CellList been built?
      *
      * IsBuilt is set true by the build() method and false by clear().
      */
      bool isBuilt() const;

      /**
      * Return true if valid, or throw Exception.
      *
      * \return true if valid, otherwise throw an Exception.
      */
      bool isValid() const;

   private:

      /*
      * Struct containing an Atom* pointer and cell id for that atom.
      *
      * The CellList::placeAtom(&Atom) method calculates the cell rank
      * for an atom, and stores the rank and a pointer in an AtomTag.
      * It also increments the number of atoms in the relevant cell.
      * The number of atoms in each cell is known only after a loop
      * over all atoms. The CellList::build() method then sorts the
      * Atom* pointers into contiguous blocks of the handles_ array.
      */
      struct AtomTag 
      {
         Atom* handle;
         int cellRank;
      };

      /// Array of offsets to neighbors.
      Cell::OffsetArray offsets_;

      /// Grid for cells.
      Grid   grid_;

      /// Array of atom tags (dimension atomCapacity_)
      DArray<AtomTag> tags_;

      /// Array of Atom handles, sorted by cell (dimension atomCapacity_).
      DArray<Atom*>   handles_;

      /// Array of Cell objects.
      DArray<Cell>  cells_;

      /// Lower coordinate bounds (local atoms).
      Vector  lower_; 

      /// Upper coordinate bounds (local atoms).
      Vector  upper_; 

      /// Length of each cell in grid
      Vector  cellLengths_; 

      /// Lower bound for nonbonded ghosts.
      Vector  lowerOuter_; 

      /// Upper coordinate bound for nonbonded ghosts.
      Vector  upperOuter_; 

      /// Pointer to first local cell (to initialize iterator).
      Cell* begin_;

      /// Total number of atoms in cell list.
      int nAtom_;

      /// Number of atoms that were not placed in cells.
      int nReject_;

      /// Has this CellList been built?
      bool isBuilt_;

      /*
      * Calculate required dimensions for cell grid.
      *
      * Called internally by allocate and makeGrid. Does not link cells or
      * calculate offsets to neighbors.
      *
      * \param lower    lower bound used to allocate array of cells.
      * \param uppper   upper bound used to allocate array of cells.
      * \param cutoff   minimum dimension of a cell in any direction
      */
      void setGridDimensions(const Vector& lower, const Vector& upper, 
                             double cutoff);

      /**
      * Return true if atomId is valid, i.e., if 0 <= 0 < atomCapacity.
      */
      bool isValidAtomId(int atomId)
      { return ( (0 <= atomId) && (atomId < tags_.capacity()) ); }

   }; 

   // Public inline method definitions:

   /*
   * Identify the cell for an Atom, based on its position.
   */
   inline int CellList::cellIndexFromPosition(const Vector& position) const
   {
      IntVector r;
      for (int i = 0; i < Dimension; ++i) {
         if (position[i] < lowerOuter_[i]) {
            return -1;
         }
         if (position[i] > upperOuter_[i]) {
            return -1;
         }
         if (position[i] < lower_[i]) {
            r[i] = 0;
         } else {
            r[i] = int( (position[i] - lower_[i])/ cellLengths_[i] ) + 1;
         }
         assert(r[i] < grid_.dimension(i));
      }
      return grid_.rank(r);
   }

   /*
   * Add an Atom to the appropriate cell, based on its position.
   */
   inline void CellList::placeAtom(Atom &atom)
   {
      // Preconditon
      assert(nAtom_ < tags_.capacity());

      int rank = cellIndexFromPosition(atom.position());
      if (rank >= 0) {
         tags_[nAtom_].cellRank = rank;
         tags_[nAtom_].handle = &atom;
         cells_[rank].incrementCapacity();
         ++nAtom_;
      } else {
         ++nReject_;
      }
   }

   inline const Grid& CellList::grid() const
   {  return grid_; }

   inline double CellList::cellLength(int i) const
   {  return cellLengths_[i]; }

   inline const Cell& CellList::cell(int i) const
   {  return cells_[i]; }

   inline const Cell* CellList::begin() const
   {  return begin_; }

   inline bool CellList::isAllocated() const
   {  return (cells_.capacity() > 0); }

}
#endif