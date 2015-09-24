/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McNMolOutput.h"                  
#include <mcMd/mcSimulation/mc_potentials.h> // include all MC potentials
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McNMolOutput::McNMolOutput(McSystem& system) :
      SystemAnalyzer<McSystem>(system)
   {  setClassName("McNMolOutput"); }

   /*
   * Read file name and open output file.
   */
   void McNMolOutput::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }
 
   /*
   * Load state from an archive, and open output file.
   */
   void McNMolOutput::loadParameters(Serializable::IArchive& ar)
   {  
      Analyzer::loadParameters(ar);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }

   /*
   * Save state to an archive.
   */
   void McNMolOutput::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Evaluate energy and output to outputFile_.
   */
   void McNMolOutput::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         int nMol = system().nMolecule(speciesId_);
         outputFile_ << nMol << std::endl;
      }
   }
 
   /* 
   * Summary
   */
   void McNMolOutput::output() 
   {
      // Close *.dat file
      outputFile_.close();

      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_ << std::endl;
      outputFile_ << std::endl;

      outputFile_ << "File format:" << std::endl;
      outputFile_    << std::endl;

      outputFile_.close();
   }
   
}
