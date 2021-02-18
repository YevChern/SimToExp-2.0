/*
 * HSVColour.cpp
 *
 *  Created on: May 30, 2013
 *      Author: bholland
 */

#include "HSVColour.h"
#include "MathUtils.h"

const uint MAX_HSV = 255;

//Default constructor, sets all to zero
HSVColour::HSVColour() {

	hue_ = 0;
	saturation_ = 0;
	value_ = 0;
}

//Constructor that sets all values - if any value outside of range, just set to zero
HSVColour::HSVColour(double hue, double saturation, double value) {

	if(hue >= 0 and hue <= 1)
		hue_ = hue;
	else
		hue_ = hue;

	if(saturation >= 0 and saturation <= 1)
		saturation_ = saturation;
	else
		saturation_ = 0;

	if(value >= 0 and value <= 1)
		value_ = value;
	else
		value_ = 0;
}

//Destructor
HSVColour::~HSVColour() {}

//getters - return signed double for compatability with other functions; still guaranteed to be positive
double HSVColour::hue() const {return hue_;}
double HSVColour::saturation() const {return saturation_;}
double HSVColour::value() const {return value_;}

//getters that return the colour values in the standard 256 integer fashion
uint HSVColour::intHue() const {return math::round(hue_ * MAX_HSV);}
uint HSVColour::intSaturation() const {return math::round(saturation_ * MAX_HSV);}
uint HSVColour::intValue() const {return math::round(value_ * MAX_HSV);}

//setters
void HSVColour::setHue(double hue) {hue_ = hue;}
void HSVColour::setSaturation(double saturation) {saturation_ = saturation;}
void HSVColour::setValue(double value) {value_ = value;}

//less than operator for sorting - in order of importance (to me anyway): hue, value, saturation
bool HSVColour::operator <(const HSVColour& other) const {

	if(hue_ < other.hue())
		return true;

	//if they happen to be equal, next compare values
	else if(hue_ == other.hue()) {

		if(value_ < other.value())
			return true;

		//if still equal, compare saturations
		else if(value_ == other.value()) {

			if(saturation_ < other.saturation())
				return true;
		}
	}

	return false;
}
