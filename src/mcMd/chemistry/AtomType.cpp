#ifndef MCMD_ATOM_TYPE_CPP
#define MCMD_ATOM_TYPE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomType.h"
#include <util/param/Parameter.h>

namespace McMd
{

   using namespace Util;

   // Constructor
   AtomType::AtomType()
    : mass_(1.0),
      #ifdef INTER_COULOMB
      charge_(0.0),
      #endif
      name_(),
      id_(-1) 
   {}

   // Destructor
   AtomType::~AtomType()
   {}

   // Set the index.
   void AtomType::setId(int id)
   {  id_ = id; }

   /* 
   * Input a AtomType from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, AtomType &atomType)
   {
      in >> atomType.name_;
      in >> atomType.mass_;
      #ifdef INTER_COULOMB
      in >> atomType.charge_;
      #endif
      return in;
   }
   
   /* 
   * Output a AtomType to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const AtomType &atomType) 
   {
      out.width(Parameter::Width);
      out << atomType.name_;
      out.setf(std::ios::scientific);
      out.width(Parameter::Width);
      out.precision(Parameter::Precision);
      out << atomType.mass_;
      #ifdef INTER_COULOMB
      out.width(Parameter::Width);
      out << atomType.charge_;
      #endif
      return out;
   }

} 

#ifdef UTIL_MPI
namespace Util
{

   template <>
   void send<McMd::AtomType>(MPI::Comm& comm, McMd::AtomType& data, int dest, int tag)
   {
      std::string  name = data.name();
      send<std::string>(comm, name, dest, tag);

      double  mass = data.mass();
      send<double>(comm, mass, dest, tag);

      #ifdef INTER_COULOMB
      double charge = data.charge();
      send<double>(comm, charge, dest, tag);
      #endif
   }

   template <>
   void recv<McMd::AtomType>(MPI::Comm& comm, McMd::AtomType& data, int source, int tag)
   {
      std::string name;
      recv<std::string>(comm, name, source, tag);
      data.setName(name);

      double  mass; 
      recv<double>(comm, mass, source, tag);
      data.setMass(mass);

      #ifdef INTER_COULOMB
      double  charge; 
      recv<double>(comm, charge, source, tag);
      data.setCharge(charge);
      #endif
   }

   template <>
   void bcast<McMd::AtomType>(MPI::Intracomm& comm, McMd::AtomType& data, int root)
   {
      std::string  name;
      double  mass; 
      #ifdef INTER_COULOMB
      double  charge; 
      #endif
      int  rank = comm.Get_rank();
      if (rank == root) {
         name = data.name();
         mass = data.mass();
         #ifdef INTER_COULOMB
         charge = data.charge();
         #endif
      }
      bcast<std::string>(comm, name, root);
      bcast<double>(comm, mass, root);
      #ifdef INTER_COULOMB
      bcast<double>(comm, charge, root);
      #endif
      if (rank != root) {
         data.setName(name);
         data.setMass(mass);
         #ifdef INTER_COULOMB
         data.setCharge(charge);
         #endif
      }
   }

   /**
   * Initialize AtomType MPI Datatype.
   */
   MPI::Datatype MpiTraits<McMd::AtomType>::type = MPI::BYTE;
   bool MpiTraits<McMd::AtomType>::hasType = false;

}
#endif

#endif
