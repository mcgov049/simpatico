#ifndef MCMD_DIBLOCK_GRAND_MOVE_H
#define MCMD_DIBLOCK_GRAND_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>   // base class
#include <mcMd/chemistry/Molecule.h>

namespace McMd
{

   using namespace Util;

   class Diblock;
   class McSystem;

   /**
   * A move that inserts of deltes a Diblock molecule.
   *
   * \ingroup McMd_McMove_Module
   */
   class DiblockGrandMove : public SystemMove
   {
   
   public:
   
      /**
      * Constructor. 
      */
      DiblockGrandMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();

      /**
      * Generate position of molecule to insert
      */
      Molecule& Insert();
   
   protected:
   
      /// Pointer to instance of Diblock.
      Diblock* speciesPtr_;

   private:

      /// Number of atoms per molecule of this species
      int nAtom_;

      /// Array of positions for new molecule
      DArray<Vector> newPositions_;
   
      /// Float for chemical potential.
      double mu_;

      /// Integer index for molecular species.
      int speciesId_;

   };

}      
#endif
