#ifndef GROUP_DISTRIBUTOR_H
#define GROUP_DISTRIBUTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/chemistry/Group.h>            // member
#include <ddMd/communicate/Buffer.h>         // member data type
#include <util/containers/DArray.h>          // member

namespace DdMd
{

   class AtomStorage;
   template <int N> class GroupStorage;
   class Domain;
   class Buffer;

   using namespace Util;

   /**
   * Class template for distributing Group<N> objects among processors.
   *
   * A GroupDistributor is used to distribute items among processors during
   * startup, when the master process must read a configuration file.  The
   * loops to distribute Groups (bonds, angles, or dihedrals) requires access
   * to an DArray<int> atomOwners[i] in which atomOwners[i] is the rank of 
   * the processor that owns atom i.
   *
   * Usage:
   * \code
   * 
   *    GroupDistributor  distributor;
   *    GroupStorage<N>   groupStorage;
   *    DArray<int>       atomOwners;
   *    Group<N>*         ptr;
   *    std::ifstream file
   *
   *    if (rank = 0) {  // If master processor
   *
   *       distributor.initSendBuffer();
   *
   *       // Read from file
   *       for (i = 0; i < nGroup; ++i) {
   *
   *           ptr = distributor.newPtr();
   *  
   *           // Read properties of Group<N> *ptr from file
   *           file >> ptr->position();
   *           // ...
   *
   *           // Cache active Group<N> for sending.
   *           distributor.add(groupStorage, atomOwners);
   *       }
   *
   *       // Send all remaining.
   *       distributor.send();
   *
   *    } else { // If not master processor
   *
   *       distributor.receive(groupStorage);
   *
   *    }
   * \endcode
   * The add(groupStorage) method can send items if required by
   * memory limits, and the send() method then sends all remaining
   * groups.
   *
   */
   template <int N>
   class GroupDistributor : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      GroupDistributor();

      /**
      * Destructor.
      */
      ~GroupDistributor();

      /**
      * Set pointers to AtomStorage.
      */
      void associate(Domain& domain, 
                     AtomStorage& atomStorage, 
                     GroupStorage<N>& groupStorage, 
                     Buffer& buffer);

      /**
      * Initialize Buffer for sending.
      */
      void initSendBuffer();

      /**
      * Set cacheCapacity, allocate memory and initialize object.
      *
      * This method sets cacheCapacity from file, and then goes through
      * the same initialization steps as readParam.
      *
      * \param cacheCapacity max number of groups cached for sending
      */
      void setParam(int cacheCapacity = -1);

      /**
      * Read cacheCapacity, allocate memory and initialize object.
      *
      * This method and setParam both allocate all required memory and 
      * initialize the GroupDistributor. Memory for a temporary read 
      * cache is alocated only on the master.
      * 
      * Preconditions: The initialize method must have been called, the 
      * Buffer must be allocated, and the Domain must be initialized.
      *
      * \param in input stream from which parameter is read.
      */
      virtual void readParam(std::istream& in);

      /**
      * Returns pointer an address available for a new Group<N>.
      *
      * This method should be called only by the master processor. It
      * returns the address within the internal cache for a new Group<N>
      * object.  Each call to newPtr() must be followed by a matching
      * call to add(). 
      *
      * \return address for a new Group<N> object.
      */
      Group<N>* newPtr(); 
     
      /**
      * Add a group to the cache for sending, send if necessary.
      *
      * If the addition of this atom to the cache would make the cache
      * full, this method broadcasts the cache to all processors.
      *
      * \param  storage GroupStorage<N> object on master processor.
      * \return rank of processor that owns the active atom.
      */
      void add();

      /**
      * Send all atoms that have not be sent previously.
      *
      * This method should be called only by the master processor, after
      * completion of a loop over all atoms to be sent.
      */
      void send();

      /**
      * Receive all atoms sent by master processor.
      *
      * This should be called by all processes except the master.
      */ 
      void receive();

   protected:

      /// Maximum number of items that can be cached on master.
      int  cacheCapacity_;

      /**
      * Allocate memory and initialize object.
      *
      * This method allocates all required memory and leaves the object ready
      * for use. It is called by setParam and readParam.
      */
      void allocate();

   private:

      /// Array that holds cached Group objects to be sent to other processors.
      /// Allocated only on the master processor.
      DArray< Group<N> > cache_;

      /// Pointer to space for a new local Group<N>. Null when inactive.
      Group<N>*   newPtr_;
      
      /// Pointer to associated Buffer object.
      Domain* domainPtr_;

      /// Pointer to associated Buffer object.
      AtomStorage* atomStoragePtr_;

      /// Pointer to associated GroupStorage<N> object.
      GroupStorage<N>* groupStoragePtr_;

      /// Pointer to associated Buffer object.
      Buffer* bufferPtr_;

      /// Type of object to send.
      Buffer::BlockDataType sendType_;

      /// Current size of cache_ (defined only on master)
      int cacheSize_;

      /// Total number of groups sent(defined only on master)
      int nSentTotal_;

      /// Total number of atoms in groups recieved (defined on all)
      int nAtomRecv_;

   };

}
#endif
