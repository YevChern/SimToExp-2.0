/*
 * HSVColour.h
 *
 *  Created on: May 30, 2013
 *      Author: bholland
 */

#ifndef HSVCOLOUR_H_
#define HSVCOLOUR_H_

typedef unsigned int uint;

class HSVColour {

	double hue_; //colour on a colour wheel
	double saturation_; //intensity of the colour, greys have value of zero
	double value_; //brightness of the colour

public:

	HSVColour();
	HSVColour(double,double,double);
	virtual ~HSVColour();

    double hue() const;
    double saturation() const;
    double value() const;

    uint intHue() const;
    uint intSaturation() const;
    uint intValue() const;

    void setHue(double);
    void setSaturation(double);
    void setValue(double);

    bool operator< (const HSVColour&) const;
};

#endif /* HSVCOLOUR_H_ */
