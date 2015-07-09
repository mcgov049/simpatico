/*MMOD
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaxWallTime.h"
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>
time_t t_start = time(0);

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MaxWallTime::MaxWallTime(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("MaxWallTime"); } 

   /*
   * Read interval and outputFileName. 
   */
   void MaxWallTime::readParameters(std::istream& in) 
   {
      readInterval(in);
      read<double>(in, "max_time", MAXT);
      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void MaxWallTime::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      #if 0
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
      isInitialized_ = true;
   }

   /*
   * Read interval and outputFileName. 
   */
   void MaxWallTime::clear() 
   {  nSample_ = 0; }

   /*
   * Dump configuration to file
   */
   void MaxWallTime::sample(long iStep) 
   {
      if (isAtInterval(iStep))  
      {
        time_t Tnow = time(0);
        double T  = Tnow - t_start;
        if (T>MAXT) {UTIL_THROW("Out of time");}
       }

         ++nSample_;
   }

}
