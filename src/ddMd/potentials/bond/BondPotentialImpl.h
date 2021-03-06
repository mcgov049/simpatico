#ifndef DDMD_BOND_POTENTIAL_IMPL_H
#define DDMD_BOND_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "BondPotential.h" // base class

namespace DdMd
{

   using namespace Util;

   class Simulation;
   template <int N> class GroupStorage;

   /**
   * Implementation template for a BondPotential.
   *
   * \ingroup DdMd_Bond_Module
   */
   template <class Interaction>
   class BondPotentialImpl : public BondPotential
   {

   public:

      /** 
      * Constructor.
      *
      * This is the primary constructor, which is used in a simulation.
      *
      * \param simulation  parent Simulation
      */
      BondPotentialImpl(Simulation& simulation);

      /** 
      * Default constructor.
      *
      * This constructor is provided to simplify unit testing.
      */
      BondPotentialImpl();

      /** 
      * Destructor.
      */
      virtual ~BondPotentialImpl();

      /// \name Initialization and Serialization
      //@{
      
      /**
      * Set the maximum number of atom types.
      *
      * This should be called by the user iff the object is instantiated
      * with the default constructor. It is called internally by the main 
      * constructor BondPotentialImpl(Simulation& ).
      *
      * \param nAtomType  maximum number of atom types (max index + 1)
      */
      virtual void setNBondType(int nAtomType);
  
      /**
      * Read potential energy parameters.
      * 
      * This method reads the bond potential Interaction parameter block. 
      *
      * \pre nBondType must be set
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \pre nBondType must be set
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      //@}
      /// \name Interaction interface
      //@{

      /**
      * Return pair energy for a single pair.
      */
      virtual double bondEnergy(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual double bondForceOverR(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual 
      double randomBondLength(Random* random, double beta, int bondTypeId) 
             const;

      /**
      * Modify a bond parameter, identified by a string.
      *
      * \param name       parameter variable name
      * \param bondTypeId type index of first atom
      * \param value      new value of parameter
      */
      void set(std::string name, int bondTypeId, double value);

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter variable name
      * \param bondTypeId  type index of first atom
      */
      double get(std::string name, int bondTypeId) const;

      /**
      * Return pair interaction class name (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return bond interaction by const reference.
      */
      const Interaction& interaction() const;

      /**
      * Return bond interaction by reference.
      */
      Interaction& interaction();

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add the bond forces for all atoms.
      */
      virtual void computeForces();

      /**
      * Compute the total bond energy for all processors
      * 
      * Call on all processors (MPI reduce operation).
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator);
      #else
      virtual void computeEnergy();
      #endif

      /**
      * Compute the covalent bond stress.
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeStress(MPI::Intracomm& communicator);
      #else
      virtual void computeStress();
      #endif

      /**
      * Compute bond forces and stress.
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeForcesAndStress(MPI::Intracomm& communicator);
      #else
      virtual void computeForcesAndStress();
      #endif

      //@}

   private:

      /**
      * Total bond energy on all processors.
      */
      double energy_;
 
      /**
      * Pointer to Interaction (evaluates the bond potential function).
      */ 
      Interaction* interactionPtr_;

      /*
      * Initialized false, reset true in (read|load)Parameters.  
      */
      bool isInitialized_;

   };

}

#include "BondPotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
//#include <util/accumulators/setToZero.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl(Simulation& simulation)
    : BondPotential(simulation),
      interactionPtr_(0),
      isInitialized_(false)
   {  
      interactionPtr_ = new Interaction(); 
      setNBondType(simulation.nBondType());
   }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl()
    : BondPotential(),
      interactionPtr_(0),
      isInitialized_(false)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::~BondPotentialImpl() 
   {
      if (interactionPtr_) {
         delete interactionPtr_;
      }
   }

   /*
   * Set the maximum number of atom types.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::setNBondType(int nBondType)
   {  interaction().setNBondType(nBondType); }

   /*
   * Read bond interaction parameters.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::readParameters(std::istream& in)
   {
      UTIL_CHECK(!isInitialized_);
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   template <class Interaction>
   void 
   BondPotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      UTIL_CHECK(!isInitialized_);
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().loadParameters(ar);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {  interaction().save(ar); }

   /*
   * Return bond energy for a single pair.
   */
   template <class Interaction>
   double 
   BondPotentialImpl<Interaction>::bondEnergy(double rsq, int bondTypeId) const
   {  return interaction().energy(rsq, bondTypeId); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction>
   double 
   BondPotentialImpl<Interaction>::bondForceOverR(double rsq, int bondTypeId)
   const
   {  return interaction().forceOverR(rsq, bondTypeId); }

   /*
   * Return random bond length.
   */
   template <class Interaction> double 
   BondPotentialImpl<Interaction>::
      randomBondLength(Random* random, double beta, int bondTypeId) const
   {  return interaction().randomBondLength(random, beta, bondTypeId); }

   /*
   * Return bond potential interaction class name.
   */
   template <class Interaction>
   std::string BondPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

   /*
   * Modify a bond parameter, identified by a string.
   */
   template <class Interaction>
   inline
   void BondPotentialImpl<Interaction>::set(std::string name, int bondTypeId, 
                                            double value)
   {  interactionPtr_->set(name, bondTypeId, value); }

   /*
   * Get a parameter value, identified by a string.
   */
   template <class Interaction>
   inline
   double BondPotentialImpl<Interaction>::get(std::string name, int bondTypeId) const
   {  return interactionPtr_->get(name, bondTypeId); }

   /*
   * Get Interaction by reference.
   */
   template <class Interaction>
   inline Interaction& BondPotentialImpl<Interaction>::interaction()
   {  return *interactionPtr_; }

   /*
   * Get Interaction by const reference.
   */
   template <class Interaction>
   inline const Interaction& BondPotentialImpl<Interaction>::interaction() const
   {  return *interactionPtr_; }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::computeForces()
   {  
      Vector f;
      double rsq;
      GroupIterator<2> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      int type, isLocal0, isLocal1;

      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         // Set f = r0 - r1, minimum image separation between atoms
         rsq = boundary().distanceSq(atom0Ptr->position(), 
                                     atom1Ptr->position(), f);
         // Set force = (r0-r1)*(forceOverR)
         f *= interactionPtr_->forceOverR(rsq, type);
         if (isLocal0) {
            atom0Ptr->force() += f;
         }
         if (isLocal1) {
            atom1Ptr->force() -= f;
         }
      }
   }

   /*
   * Compute total bond energy on all processors, store result on master.
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void 
   BondPotentialImpl<Interaction>::computeEnergy(MPI::Intracomm& communicator)
   #else
   void BondPotentialImpl<Interaction>::computeEnergy()
   #endif
   { 

      // If energy is already set, do nothing and return.
      if (isEnergySet()) return;

      Vector f;
      double rsq;
      double bondEnergy;
      double localEnergy = 0.0;
      GroupIterator<2> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      int type, isLocal0, isLocal1;

      // Loop over bonds to compute localEnergy on this processor.
      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);

         // Calculate minimum image square distance between atoms
         rsq = boundary().distanceSq(atom0Ptr->position(), 
                                     atom1Ptr->position());
         bondEnergy = interactionPtr_->energy(rsq, type);

         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         if (isLocal0 && isLocal1) {
            localEnergy += bondEnergy;
         } else {
            assert(isLocal0 || isLocal1);
            localEnergy += 0.5*bondEnergy;
         }
      }

      // Add localEnergy from all processors, set energy to sum on master.
      reduceEnergy(localEnergy, communicator);
   }

   /*
   * Compute total pair stress (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void BondPotentialImpl<Interaction>::computeStress(MPI::Intracomm& communicator)
   #else
   void BondPotentialImpl<Interaction>::computeStress()
   #endif
   {
      // Do nothing and return if stress is already set.
      if (isStressSet()) return;
 
      Tensor localStress;
      Vector dr;
      Vector f;
      double rsq, forceOverR;
      GroupIterator<2> iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type;
      int    isLocal0, isLocal1;

      localStress.zero();

      // Iterate over bonds
      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         rsq = boundary().distanceSq(atom0Ptr->position(), 
                                        atom1Ptr->position(), dr);
         f = dr;
         assert(isLocal0 || isLocal1);
         if (isLocal0 && isLocal1) {
            forceOverR = interactionPtr_->forceOverR(rsq, type);
         } else {
            forceOverR = 0.5*interactionPtr_->forceOverR(rsq, type);
         }
         f *= forceOverR;
         incrementPairStress(f, dr, localStress);
      }

      // if (reverseUpdateFlag()) { } else { } 

      // Normalize by volume 
      localStress /= boundary().volume();

      // Add localStress from all processors, set stress to sum on master
      reduceStress(localStress, communicator);
   }

   /*
   * Compute total pair stress (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void BondPotentialImpl<Interaction>::computeForcesAndStress(MPI::Intracomm& communicator)
   #else
   void BondPotentialImpl<Interaction>::computeForcesAndStress()
   #endif
   {
      // Do nothing and return if stress is already set.
      if (isStressSet()) {
         computeForces();
         return;
      }
 
      Tensor localStress;
      Vector dr;
      Vector f;
      double rsq;
      GroupIterator<2> iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type;
      int    isLocal0, isLocal1;

      localStress.zero();

      // Iterate over bonds
      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         rsq = boundary().distanceSq(atom0Ptr->position(), 
                                        atom1Ptr->position(), dr);
         f  = dr;
         f *= interactionPtr_->forceOverR(rsq, type);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         assert(isLocal0 || isLocal1);
         if (isLocal0) {
            atom0Ptr->force() += f;
         }
         if (isLocal1) {
            atom1Ptr->force() -= f;
         }
         if (!(isLocal0 && isLocal1)) {
            f *= 0.5;
         }
         incrementPairStress(f, dr, localStress);
      }

      // if (reverseUpdateFlag()) { } else { } 

      // Normalize by volume 
      localStress /= boundary().volume();

      // Add localStress from all processors, set stress to sum on master
      reduceStress(localStress, communicator);
   }

}
#endif
