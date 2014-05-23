#ifndef MDPP_PROCESSOR_H
#define MDPP_PROCESSOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class
#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function

#include <mdPp/chemistry/Atom.h>              // member (template argument)
#include <mdPp/chemistry/Group.h>             // member (template argument)
#include <mdPp/processor/GroupStorage.h>      // member 
#include <mdPp/configIos/ConfigIoFactory.h>   // member 
#include <mdPp/analyzers/AnalyzerManager.h>   // member 

namespace MdPp 
{

   class ConfigIo;

   using namespace Util;

   /**
   * A post-processor for analyzing outputs of MD simulations.
   */
   class Processor : public ParamComposite 
   {

   public:

      typedef ArrayIterator<Atom> AtomIterator;
      typedef ArrayIterator<Group <2> > BondIterator;

      /**
      * Constructor
      */
      Processor();

      /**
      * Destructor
      */
      ~Processor();

      using ParamComposite::readParam;

      /**
      * Open, read, and close parameter file.
      */
      void readParam(const char* filename);

      /**
      * Read parameters.
      */
      void readParameters(std::istream& in);

      /**
      * Set ConfigIo style.
      */
      void setConfigIo(std::string configIoName);

      /**
      * Return the ConfigIo (create default if necessary).
      */
      ConfigIo& configIo();
   
      /**
      * Read a single configuration file.
      */
      void readConfig(std::ifstream& in);

      /**
      * Open, read and close a configuration file.
      */
      void readConfig(const char* filename);
   
      /**
      * Open, read and close a configuration file.
      */
      void readConfig(const std::string& filename);
   
      /**
      * Write a single configuration file.
      */
      void writeConfig(std::ofstream& out);

      /**
      * Open, write and close a configuration file.
      */
      void writeConfig(const std::string& filename);

      /**
      * Read and analyze a sequence of numbered configuration files.
      *
      * This function reads and analyzes a sequence of configuration files 
      * that were generated by running a previous simulation. The function 
      * reads files with names of the form inputPrefix() + n for integer 
      * suffixes min <= n <= max. 
      *
      * \param min  integer suffix of first configuration file name
      * \param max  integer suffix of last configuration file name
      * \param fileBaseName root name for dump files (without integer suffix)
      */  
      void analyzeDumps(int min, int max, std::string fileBaseName);

      /**
      * Analyze a trajectory file.
      */
      void analyzeTrajectory(std::string& filename);

      /**
      * Clear all atoms and groups.
      */
      void clear();
  
      /**
      * Get the Boundary by non-const reference
      */
      Boundary& boundary();

      // Atom container interface

      /**
      * Return pointer to location for new atom.
      *
      * \param  global id for new atom
      * \return pointer to location of new atom
      */
      Atom* newAtomPtr();

      /**
      * Finalize addition of atom (allows lookup by id).
      */
      void addAtom();

      /**
      * Get a pointer to an atom by global id.
      */
      Atom* atomPtr(int id);

      /**
      * Initialize an iterator for atoms.
      */
      void initAtomIterator(AtomIterator& iter);

      /**
      * Get atom capacity (maximum id + 1).
      */ 
      int atomCapacity() const;

      /**
      * Get number of atoms.
      */ 
      int nAtom() const;

      // Group storage interface

      #ifdef INTER_BOND
      GroupStorage<2>& bonds();
      #endif

      #ifdef INTER_ANGLE
      GroupStorage<3>& angles();
      #endif

      #ifdef INTER_DIHEDRAL
      GroupStorage<4>& dihedrals();
      #endif

      // etc. for angles and dihedrals

      // Accessors, for use in Analyzer and ConfigIo classes.

   private:
     
      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// Array of atom objects, added in order read from file.
      DSArray<Atom> atoms_;

      /// Pointers to atoms indexed by ids. Missing atoms are null pointers.
      DArray<Atom*> atomPtrs_;

      #ifdef INTER_BOND
      /// Array of bond objects, added in order read from file.
      GroupStorage<2> bonds_;
      #endif

      #ifdef INTER_ANGLE
      /// Array of angle objects, added in order read from file.
      GroupStorage<3> angles_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Array of dihedral objects, added in order read from file.
      GroupStorage<4> dihedrals_;
      #endif

      /// Pointer to current ConfigIo object.
      ConfigIo* configIoPtr_;

      /// Factory for generating ConfigIo at run time.
      ConfigIoFactory configIoFactory_;

      /// Manager for analyzers
      AnalyzerManager analyzerManager_;

      /// Pointer to new atom.
      Atom* newAtomPtr_;

      /// Maximum allowed atom id + 1 (used to allocate arrays).
      int atomCapacity_;

      #ifdef INTER_BOND
      /// Maximum number of bonds (used to allocate array).
      int bondCapacity_;
      #endif

      #ifdef INTER_ANGLE
      /// Maximum number of angles (used to allocate array).
      int angleCapacity_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Maximum number of dihedrals (used to allocate array).
      int dihedralCapacity_;
      #endif

      /// String identifier for ConfigIo class name
      std::string configIoName_;

      /// Name of configuration or trajectory input file
      std::string configFileName_;

   };

   // inline functions

   /*
   * Return the Boundary by reference.
   */
   inline Boundary& Processor::boundary() 
   {  return boundary_; }

   /*
   * Return number of atoms.
   */
   inline int Processor::nAtom() const
   {  return atoms_.size(); }

   /*
   * Get atom capacity (maximum id + 1).
   */ 
   inline
   int Processor::atomCapacity() const
   { return atoms_.capacity(); }

   /*
   * Return a pointer to an atom with a specific id.
   */
   inline Atom* Processor::atomPtr(int id)
   {  return atomPtrs_[id]; }

   /*
   * Initialize an iterator for atoms.
   */
   inline 
   void Processor::initAtomIterator(Processor::AtomIterator& iter)
   {  atoms_.begin(iter); }

   #ifdef INTER_BOND
   inline GroupStorage<2>& Processor::bonds()
   {  return bonds_; }
   #endif

   #ifdef INTER_ANGLE
   inline GroupStorage<3>& Processor::angles()
   {  return angles_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline GroupStorage<4>& Processor::dihedrals()
   {  return dihedrals_; }
   #endif

}
#endif