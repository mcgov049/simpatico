#ifndef DDMD_MAX_WALL_TIME_H
#define DDMD_MAX_WALL_TIME_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>
#include <util/mpi/MpiLoader.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write (scalar) pressure to file.
   *
   * \sa \ref ddMd_analyzer_OutputRegionPressure_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class MaxWallTime : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      MaxWallTime(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~MaxWallTime()
      {} 
   
      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Dump configuration to file
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);


   private:
 
      /// Has readParam been called?
      long    isInitialized_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;

      double MAXT;
   
   };

}
#endif 
