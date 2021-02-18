/*

 For details see XrayAtomicFormFactor.h

 */

#include "XrayAtomicFormFactor.h"

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

const double FOUR_PI_SQUARED = 157.91367041742973790135185599802;
const unsigned int ITERATIONS = 4; //Summation iterations for the analytic form factor expression

//Default constructor - not intended for use at the moment
XrayAtomicFormFactor::XrayAtomicFormFactor() {

	isUnitedAtom = false;
	atomicFF_calculated = false;
	atomicFF = 0;
	ua_constant = 0;
	type = UNKNOWN;
}

//Constructor that sets the constants using the AtomType and the atomic information
XrayAtomicFormFactor::XrayAtomicFormFactor(AtomType type, const AtomicInfo& info) {

	atomicFF_calculated = false;
	atomicFF = 0;

   this->type = type;
   isUnitedAtom = info.isUnitedAtom(type);
   
   //set the constant(s) based on whether a united or regular atom
   if(isUnitedAtom) {

      ua_constant = info.UAFormFactorConstant(type);
   }
   else {

      constants = info.atomicFormFactorConstants(type);
      ua_constant = 0;
   }
}

//Destructor
XrayAtomicFormFactor::~XrayAtomicFormFactor() {}

//returns the AtomType of the form factor
AtomType XrayAtomicFormFactor::getAtomType() const {return type;}

//function that calculates the atomic form factor at 'q'
//There are 4 'a', 4 'b' and 1 'c' constant, thus the array indexes below (i.e. a: 0-3, b: 4-7, c: 8)
double XrayAtomicFormFactor::calcAtomicFF(double q) {

	//if already calculated just return the value
	if(atomicFF_calculated)
		return atomicFF;

	//otherwise need to calculate
	//keep track of total
	double fi = 0.0;

	if(isUnitedAtom) {

		fi += ua_constant;

	} else {

		for (uint i = 0; i < ITERATIONS; ++i)
			fi += constants.at(i)*exp(-constants.at(i+4)*q*q/FOUR_PI_SQUARED);

		fi += constants.at(8);
	}

    //return
	atomicFF = fi;
	return atomicFF;
}

void XrayAtomicFormFactor::setAtomicFF_calculated(bool val) {
    atomicFF_calculated = val;
}
