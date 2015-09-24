/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DiblockGrandMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/species/Diblock.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Bond.h>
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   DiblockGrandMove::DiblockGrandMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1),
      mu_(0.0)
   {  setClassName("DiblockGrandMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void DiblockGrandMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "mu", mu_);

      nAtom_ = simulation().species(speciesId_).nAtom();
      newPositions_.allocate(nAtom_);
      // Cast the Species Diblock
      speciesPtr_ = dynamic_cast<Diblock*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be Diblock");
      }
  
   }
 
   /*
   * Load state from an archive.
   */
   void DiblockGrandMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<double>(ar, "mu", mu_);

      nAtom_ = simulation().species(speciesId_).nAtom();
      newPositions_.allocate(nAtom_);
      // Cast the Species to Diblock
      speciesPtr_ = dynamic_cast<Diblock*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Species is not a Diblock");
      }
   }

   /*
   * Save state to an archive.
   */
   void DiblockGrandMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;  
      ar & mu_;
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool DiblockGrandMove::move() 
   {
      double     MolEnergy;
      int        iAtom;
      Molecule*  molPtr;
      Atom*      atomPtr;
      double     volume_;
      int        nMol;
      bool accept = true;

      incrementNAttempt();
      //Decide whether to remove or add
      int AddRm = random().uniformInt(0,2); //1=Add, 0=Remove
      nMol = system().nMolecule(speciesId_);
      volume_ = system().boundary().volume();

      MolEnergy = 0;
      if(AddRm==0) 
      {
        if (nMol == 0) 
        {
           bool accept = false;
           return accept;
        }
        Molecule& molecule = system().randomMolecule(speciesId_);
        molPtr = &molecule;
        for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
          atomPtr = &molPtr->atom(iAtom);
          #ifndef INTER_NOPAIR
          MolEnergy += system().pairPotential().atomEnergy(*atomPtr);
          #endif
//          MolEnergy += system().bondPotential().atomEnergy(*atomPtr);
          #ifdef INTER_EXTERNAL
          MolEnergy += system().externalPotential().atomEnergy(*atomPtr);
          #endif
          #ifdef INTER_TETHER
          MolEnergy += system().atomTetherEnergy(*atomPtr);
          #endif
        }
        accept = random().metropolis((nMol/volume_)*boltzmann(mu_ - MolEnergy));
        if (accept) {
          system().removeMolecule(molecule);
          simulation().returnMolecule(molecule);
          for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
               atomPtr = &molPtr->atom(iAtom);
               system().pairPotential().deleteAtom(*atomPtr);
           }
          incrementNAccept();
        }
      }
      else
      {
        Molecule& molecule = Insert();
        molPtr = &molecule;
        for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
          atomPtr = &molPtr->atom(iAtom);
          #ifndef INTER_NOPAIR
          MolEnergy += system().pairPotential().atomEnergy(*atomPtr);
          #endif
//          MolEnergy += system().bondPotential().atomEnergy(*atomPtr);
          #ifdef INTER_EXTERNAL
          MolEnergy += system().externalPotential().atomEnergy(*atomPtr);
          #endif
          #ifdef INTER_TETHER
          MolEnergy += system().atomTetherEnergy(*atomPtr);
          #endif
        }
        accept = random().metropolis((volume_/(nMol+1))*boltzmann(MolEnergy - mu_));
        if (accept) {
          incrementNAccept();
        }
        else {
          system().removeMolecule(molecule);
          simulation().returnMolecule(molecule);
          for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
               atomPtr = &molPtr->atom(iAtom);
               system().pairPotential().deleteAtom(*atomPtr);
           }
        }

      }


      return accept;
   }
   
   Molecule& DiblockGrandMove::Insert()
   {
      Vector     StartPos, bondVector;
      Molecule*  molPtr;
      Atom*      atomPtr;
      Atom*      oldAtomPtr;
      double     beta;
      double     bondL;

      molPtr = &(simulation().getMolecule(speciesId_));
      Molecule& molecule = *molPtr;
      system().addMolecule(molecule);
      beta = energyEnsemble().beta();
      int bondTypeId = molPtr->bond(0).typeId();
      bondL = system().bondPotential().randomBondLength(&random(),beta,bondTypeId);
      Vector Bounds = system().boundary().lengths();
      for (int i=0; i<3; i++)
      {
         StartPos[i] = random().uniform(0,Bounds[i]);
      }

      oldAtomPtr = &molPtr->atom(0);
      oldAtomPtr->position() = StartPos;
      system().pairPotential().addAtom(*oldAtomPtr);
      
      for (int iAtom=1; iAtom<nAtom_; iAtom++)
      {
         atomPtr = &molPtr->atom(iAtom);
         random().unitVector(bondVector);
         bondVector *= bondL;
         atomPtr->position() = oldAtomPtr->position();
         atomPtr->position() += bondVector;
         boundary().shift(atomPtr->position());
         system().pairPotential().addAtom(*atomPtr);
         oldAtomPtr = atomPtr;
      }
      
      return molecule;
   }

}
