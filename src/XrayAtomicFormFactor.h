/*
 * AtomicFormFactorIntegrand.h
 *
 *  Created on: Oct 5, 2012
 *      Author: bholland
 *
 * Function for determining the X-ray atomic form factor for a particular
 * AtomType
 *
 * The members 'a', 'b' and 'c' are constants used to calculate the form factor from 
 * a (presumably) approximated analytical expression, see: Cromer and Mann (1968) Acta
 * Crystallogr A, 24: 321-324
 * 
 * Arrays have been replaced by 'vector' objects, more versatile and safer, and not much
 * slower (although more memory, but not a big deal).
 *
 */

#ifndef XRAYATOMICFORMFACTOR_H_
#define XRAYATOMICFORMFACTOR_H_

//C++ includes
#include <vector>

//system includes
#include "AtomType.h"
#include "AtomicInfo.h"

class XrayAtomicFormFactor {
   
   std::vector<double> constants; //get this from AtomicInfo
   double ua_constant;
   AtomType type;
   bool isUnitedAtom;
   
   bool atomicFF_calculated;
   double atomicFF;

   public:
      
      XrayAtomicFormFactor();
      XrayAtomicFormFactor(AtomType,const AtomicInfo&);
      virtual ~XrayAtomicFormFactor();
      
      AtomType getAtomType() const;
      
      double calcAtomicFF(double);
      void setAtomicFF_calculated(bool);
};



#endif /* ATOMICFORMFACTORINTEGRAND_H_ */ 
