#ifndef MCMD_EE_SIMULATION_H
#define MCMD_EE_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/Simulation.h>   // base class
#include <mcMd/mcSimulation/McSystem.h>   // class member
#include <util/global.h>

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   using namespace Util;

   class McMove;
   class McMoveManager;
   class EeAnalyzerManager;

   /**
   * An Extended Ensemble Monte-Carlo simulation of one McSystem.
   *
   * \ingroup McMd_Simulation_Module
   */
   class EeSimulation : public Simulation
   {

   public:

      #ifdef UTIL_MPI
      /**
      * Constructor.
      */
      EeSimulation(MPI::Intracomm& communicator);
      #endif

      /**
      * Constructor.
      */
      EeSimulation();

      /**
      * Destructor.
      */
      virtual ~EeSimulation();

      /**
      * Read parameters from a specific stream.
      *
      * \param in parameter file input stream.
      */
      virtual void readParam(std::istream &in);

      /**
      * Read parameters from the default parameter stream.
      *
      * Default parameter istream is std::cin in serial mode 
      * (ifndef UTIL_MPI) and the file "n/param" for 
      * processor n in parallel mode (ifdef UTIL_MPI).
      */
      void readParam();

      /**
      * Read and execute commands from a specific input stream.
      * 
      * \param in command file input stream. 
      */
      void readCommands(std::istream& in);

      /**
      * Read and execute commands from the default parameter file.
      * 
      * This method opens the file with the name commandFile read
      * by the FileMaster.
      */
      void readCommands();

      /**
      * Run an MC simulation of specified length.
      * 
      * This method implements the main MC loop. The step counter iStep_
      * is incremented until it reaches endStep. Each step involves a
      * random selection and attempt of one Markov MC move. Upon exit,
      * iStep_ = endStep.
      *
      * If isContinuation is false, the step counter iStep_ is initialized 
      * to zero, and analyzers and mcmoves are set to default initial 
      * states before entering the main loop. If isContinuation is true, 
      * no such initialization is done for iStep_, analyzers, or the
      * MC moves.
      *  
      * \param endStep        Final value of MC step counter iStep_.
      * \param isContinuation Is this a continuation of a previous run?
      */
      void simulate(int endStep, bool isContinuation = false);

      /**
      * Read and analyze a sequence of configuration files.
      *
      * This method reads and analyzes a sequence of configuration files,
      * which were normally generated by running a previous simulation using 
      * DumpConfig, and applies the sample() method of every Analyzer to
      * each such configuration. 
      *
      * The method reads files with names of the form inputPrefix() + n for 
      * integer suffixes min <= n <= max. This is consistent with the output
      * format used by the DumpConfig class.
      *
      * In serial mode, the inputPrefix should be given as a path relative
      * to the directory in which the program is executed. The inputPrefix 
      * of the simulation FileMaster is not prepended to the dump Prefix.  
      *
      * In parallel mode, for processor with MPI rank m, the path "m/" is 
      * prepended to the dumpPrefix, so that all files associated with this
      * processor are in this directory, but no inputPrefix is added after
      * the string "m/".
      *
      * \param min        integer suffix of first configuration file name
      * \param max        integer suffix of last configuration file name
      * \param dumpPrefix root name for dump files (without integer suffix)
      */  
      void analyze(int min, int max, std::string dumpPrefix);

      /**
      * Write restart files.
      */
      void writeRestart(const std::string& filename);

      /**
      * Read restart files.
      */
      void readRestart(const std::string& filename);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, unsigned int version);

      /**
      * Get the McSystem by reference.
      */
      McSystem& system();

      /**
      * Get the McSystem by const refererence.
      */
      const McSystem& system() const;

      /**
      * Get the McMove factory by reference.
      */
      Factory<McMove>& mcMoveFactory();

      /**
      * Return true if valid, or throw an Exception. 
      */
      virtual bool isValid() const;

   protected:

      /**
      * Get the McMoveManager by reference.
      */
      McMoveManager& mcMoveManager();

   private:
   
      /// System.
      McSystem       system_;
   
      /// Pointer to Manager for Monte Carlo moves.
      McMoveManager* mcMoveManagerPtr_;

      /// Pointer to Manager for analyzers.
      EeAnalyzerManager* mcAnalyzerManagerPtr_;

      /// Pointer to parameter file passed to readParam(istream&)
      std::istream*   paramFilePtr_;

      /// Has readParam been called?
      bool isInitialized_;

      /// Is this EeSimulation in the process of restarting?
      bool isRestarting_;

   }; 

   // Inline Methods

   /* 
   * Get the McSystem.
   */
   inline McSystem& EeSimulation::system()
   { return system_; }

   /* 
   * Get a const ref to the McSystem.
   */
   inline const McSystem& EeSimulation::system() const
   { return system_; }

   /* 
   * Get the McMoveManager (protected).
   */
   inline McMoveManager& EeSimulation::mcMoveManager()
   {
      assert(mcMoveManagerPtr_);  
      return *mcMoveManagerPtr_; 
   }

}    
#endif
