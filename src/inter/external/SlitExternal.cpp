#ifndef INTER_SLIT_EXTERNAL_CPP
#define INTER_SLIT_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SlitExternal.h"

#include <iostream>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   SlitExternal::SlitExternal() 
    : boundaryPtr_(0),
      nAtomType_(0),
      isInitialized_(false)
   { setClassName("SlitExternal"); }
   
   /* 
   * Copy constructor.
   */
   SlitExternal::SlitExternal(const SlitExternal& other)
    : epsilon_(other.epsilon_),
      sigma_(other.sigma_),
      sigmaCb_(other.sigmaCb_),
      cutoff_(other.cutoff_),
      coeff_(other.coeff_),
      boundaryPtr_(other.boundaryPtr_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {}
   
   /* 
   * Assignment operator.
   */
   SlitExternal& SlitExternal::operator = (const SlitExternal& other)
   {
      boundaryPtr_   = other.boundaryPtr_;
      nAtomType_     = other.nAtomType_;
      isInitialized_ = other.isInitialized_;
      epsilon_       = other.epsilon_;
      sigma_         = other.sigma_;
      sigmaCb_       = other.sigmaCb_;
      cutoff_        = other.cutoff_;
      coeff_         = other.coeff_;
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void SlitExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > SlitExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   /* 
   * Set pointer to the Boundary.
   */
   void SlitExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /**
   * Set epsilon to the externalParameter.
   */
   void SlitExternal::setExternalParameter(double externalParameter)
   {
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("SlitExternal potential is not initialized");
      }

      epsilon_ = externalParameter;
   }

   /* 
   * Read potential parameters from file.
   */
   void SlitExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
   
      // Read parameters
      read<double>(in, "epsilon", epsilon_);
      read<double>(in, "sigma", sigma_);
      read<double>(in, "cutoff", cutoff_);

      // Initialize dependent variables (particle number density = 0.7)
      sigmaCb_ = sigma_ * sigma_ * sigma_;
      coeff_   = 4.0 * epsilon_ * acos(-1.0) / 45.0 * 0.7 * sigmaCb_;

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void SlitExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_; 
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before loadParameters");
      }
      loadParameter<double>(ar, "epsilon", epsilon_);
      loadParameter<double>(ar, "sigma", sigma_);
      loadParameter<double>(ar, "cutoff", cutoff_);
      ar >> sigmaCb_;
      ar >> coeff_;
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void SlitExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << epsilon_;
      ar << sigma_;
      ar << cutoff_;
      ar << sigmaCb_;
      ar << coeff_;
   }

   /**
   * Return epsilon for externalParameter.
   */
   double SlitExternal::externalParameter() const
   {
     return epsilon_;
   }

   /**
   * Return name string "SlitExternal" for this evaluator class.
   */
   std::string SlitExternal::className() const
   {  return std::string("SlitExternal"); }
 
} 
#endif
