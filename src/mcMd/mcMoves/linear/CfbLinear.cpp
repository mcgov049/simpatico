/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include "CfbLinear.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/species/Linear.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   CfbLinear::CfbLinear(McSystem& system) :
      SystemMove(system),
      speciesId_(-1),
      nTrial_(-1)
   {
      #ifndef INTER_BOND
      UTIL_THROW("CfbLinear requires that bonds be enabled");
      #endif
   }

   /*
   * Destructor
   */
   CfbLinear::~CfbLinear()
   {}

   /*
   * Read and validate parameters speciesId and nTrial.
   */
   void CfbLinear::processParameters()
   {
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      // Check that Species is a subclass of Linear
      Species& species = simulation().species(speciesId_);
      Linear* linearPtr = dynamic_cast<Linear*>(&species);
      if (linearPtr == 0) {
         UTIL_THROW("Species is not Linear in CfbLinear");
      }

      #ifdef INTER_ANGLE
      hasAngles_ = system().hasAnglePotential() && species.nAngle() > 0;
      #endif
      #ifdef INTER_DIHEDRAL
      hasDihedrals_ = system().hasAnglePotential() && species.nDihedral() > 0;
      #endif
      #ifdef INTER_EXTERNAL
      hasExternal_ = system().hasExternalPotential();
      #endif

      // Read and validate parameter nTrial.
      if (nTrial_ <=0) {
         UTIL_THROW("Invalid parameter value, nTrial <= 0");
      }
      if (nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid parameter value, nTrial > MaxTrial_");
      }
   }

   /*
   * Read and validate parameters speciesId and nTrial.
   */
   void CfbLinear::readParameters(std::istream& in)
   {
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nTrial", nTrial_);
      processParameters();
   }

   /*
   * Load and validate parameters speciesId and nTrial.
   */
   void CfbLinear::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "nTrial", nTrial_);
      processParameters();
   }

   /*
   * Save parameters speciesId and nTrial.
   */
   void CfbLinear::save(Serializable::OArchive& ar)
   {
      ar & speciesId_;
      ar & nTrial_;
   }

   /*
   * Configuration bias algorithm for deleting one atom from chain end.
   */
   void
   CfbLinear::deleteAtom(Molecule& molecule, int atomId, 
                         int sign, double &rosenbluth, double &energy)
   {
      // sign == -1, shift == 0 -> atomId = 0, 1, 2, ... 
      // sign == +1, shift == 1 -> atomId = nAtom-1, nAtom-2, ...
      assert(sign == 1 || sign == -1);
      int shift = (sign == 1) ? 1 : 0;
      Atom& atom0 = molecule.atom(atomId);
      Atom& atom1 = molecule.atom(atomId - sign);
      Vector& pos0 = atom0.position();
      Vector& pos1 = atom1.position();

      // Calculate bond length and energy
      Vector v1; // Vector between atoms 0 and 1
      Vector u1; // Unit vector parallel to v1
      double r1, bondEnergy;
      int bondTypeId = molecule.bond(atomId - shift).typeId();
      r1 = boundary().distanceSq(pos1, pos0, v1); 
      // Here r1 = bond length squared, v1 = pos1 - pos0
      bondEnergy = system().bondPotential().energy(r1, bondTypeId);
      r1 = sqrt(r1); // bond length
      u1 = v1;
      u1 /= r1;      // unit vector

      #ifndef INTER_NOPAIR
      McPairPotential& pairPotential = system().pairPotential();
      energy = pairPotential.atomEnergy(atom0);
      #else
      energy = 0.0;
      #endif

      #ifdef INTER_ANGLE
      Vector v2, u2;
      double r2, cosTheta;
      Vector* pos2Ptr;
      AnglePotential& anglePotential = system().anglePotential();
      int angleTypeId;
      bool hasAngle = false;
      if (hasAngles_) {
         int atom2Id = atomId - 2*sign;
         if (shift) {
            assert(sign == 1);
            hasAngle = (atom2Id >= 0);
         } else {
            assert(sign == -1);
            hasAngle = (atom2Id < molecule.nAtom());
         }
         if (hasAngle) {
            Atom& atom2 = molecule.atom(atom2Id);
            pos2Ptr = &(atom2.position());
            r2 = boundary().distanceSq(*pos2Ptr, pos1, v2); // v2=p2-p1
            r2 = sqrt(r2); // r2 = |v2| = length of bond 2
            u2 = v2;
            u2 /= r2; // unit vector
            cosTheta = u1.dot(u2);
            angleTypeId = molecule.bond(atomId - 2*shift).typeId();
            energy += anglePotential.energy(cosTheta, angleTypeId);
         }
      }
      #endif

      #ifdef INTER_EXTERNAL
      ExternalPotential& externalPotential = system().externalPotential();
      if (hasExternal_) {
         energy += externalPotential.atomEnergy(atom0);
      }
      #endif

      // Here: energy = total energy - bond energy
      // Rosenbluth factor = exp(-beta*(pair + angle + external))
      rosenbluth = boltzmann(energy);

      // Add bond energy to total energy in current position
      energy += bondEnergy;

      // Loop over nTrial - 1 additional trial positions:
      double trialEnergy;
      for (int iTrial = 0; iTrial < nTrial_ - 1; ++iTrial) {

         // Generate trial position
         random().unitVector(u1);
         v1 = u1;
         v1 *= r1;
         pos0.subtract(pos1, v1);
         boundary().shift(pos0);

         // Compute trial energy (excluding bond energy)
         #ifndef INTER_NOPAIR
         trialEnergy = system().pairPotential().atomEnergy(atom0);
         #else
         trialEnergy = 0.0;
         #endif
         #ifdef INTER_ANGLE
         if (hasAngle) {
            cosTheta = u1.dot(u2);
            trialEnergy += anglePotential.energy(cosTheta, angleTypeId);
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (hasExternal_) {
            trialEnergy += externalPotential.atomEnergy(atom0);
         }
         #endif

         rosenbluth += boltzmann(trialEnergy);
      }

   }

   /*
   * Configuration bias algorithm for adding one atom to a chain end.
   */
   void
   CfbLinear::addAtom(Molecule& molecule, Atom& atom0, Atom& atom1, int atomId, 
                      int sign, double &rosenbluth, double &energy)
   {
      // sign == -1, shift == 0 -> atomId = 0, 1, 2, ... 
      // sign == +1, shift == 1 -> atomId = nAtom-1, nAtom-2, ...
      assert(sign == 1 || sign == -1);
      int shift = (sign == 1) ? 1 : 0;

      assert(atom0.typeId() == molecule.species().atomTypeId(atomId));
      Vector& pos0 = atom0.position();
      Vector& pos1 = atom1.position();

      // Generate a random bond length r1
      double beta, r1;
      int bondTypeId, iTrial;
      bondTypeId = molecule.bond(atomId - shift).typeId();
      beta = energyEnsemble().beta();
      r1 = system().bondPotential().randomBondLength(&random(), beta, bondTypeId);

      #ifndef INTER_NOPAIR
      McPairPotential& pairPotential = system().pairPotential();
      #endif

      #ifdef INTER_ANGLE
      // Calculate vector v2 = pos2 - pos1, r2 = |v2|
      Vector v2, u2;
      double r2, cosTheta;
      Vector* pos2Ptr;
      AnglePotential& anglePotential = system().anglePotential();
      int angleTypeId;
      bool hasAngle = false;
      if (hasAngles_) {
         int atom2Id = atomId - 2*sign;
         if (shift) {
            assert(sign == 1);
            hasAngle = (atom2Id >= 0);
         } else {
            assert(sign == -1);
            hasAngle = (atom2Id < molecule.nAtom());
         }
         if (hasAngle) {
            angleTypeId = molecule.bond(atomId - 2*shift).typeId();
            Atom* atom2Ptr = &atom1 - sign;
            pos2Ptr = &(atom2Ptr->position());
            r2 = boundary().distanceSq(*pos2Ptr, pos1, v2);
            r2 = sqrt(r2);
            u2 = v2;
            u2 /= r2;  // unit vector
         }
      }
      #endif
      #ifdef INTER_EXTERNAL
      ExternalPotential& externalPotential 
                                      = system().externalPotential();
      #endif

      // Loop over nTrial trial positions:
      Vector v1, u1;
      Vector trialPos[MaxTrial_];
      double trialProb[MaxTrial_], trialEnergy[MaxTrial_];
      rosenbluth = 0.0;
      for (iTrial = 0; iTrial < nTrial_; ++iTrial) {

         // Generate trial bond vector v1 and position pos0
         random().unitVector(u1);
         v1 = u1;
         v1 *= r1;
         trialPos[iTrial].subtract(pos1, v1);
         boundary().shift(trialPos[iTrial]);
         pos0 = trialPos[iTrial];

         // Compute trial energy (excluding bond energy)
         #ifndef INTER_NOPAIR
         trialEnergy[iTrial] = pairPotential.atomEnergy(atom0);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif
         #ifdef INTER_ANGLE
         if (hasAngle) {
            cosTheta = u1.dot(u2);
            trialEnergy[iTrial] += 
                        anglePotential.energy(cosTheta, angleTypeId);
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (hasExternal_) {
            trialEnergy[iTrial] +=
                        externalPotential.atomEnergy(atom0);
         }
         #endif

         // Compute unnormalized probability, increment rosenbluth 
         trialProb[iTrial] = boltzmann(trialEnergy[iTrial]);
         rosenbluth += trialProb[iTrial];
      }

      // Normalize trial probabilities
      for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
         trialProb[iTrial] = trialProb[iTrial]/rosenbluth;
      }

      // Choose a trial 
      iTrial = random().drawFrom(trialProb, nTrial_);

      // Calculate total energy for chosen trial, including bond energy 
      energy = trialEnergy[iTrial];
      energy += system().bondPotential().energy(r1*r1, bondTypeId);

      // Set position to chosen value
      pos0 = trialPos[iTrial];
   }

}
