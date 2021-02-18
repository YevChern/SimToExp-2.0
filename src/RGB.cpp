
#include "RGB.h"
#include "MathUtils.h"

using namespace math;

const unsigned int MAX_RGB = 255;

//Default constructor, sets to black;
RGB::RGB() {

	red = 0;
	blue = 0;
	green = 0;
}

//constructor that sets the colours as provided.  Any values NOT between 0 and 1 are set to zero.
RGB::RGB(double red, double green, double blue) {
   
   if(red >= 0 and red <= 1)
      this->red = red;
   else
      this->red = 0;
   
   if(green >= 0 and green <= 1)
      this->green = green;
   else
      this->green = 0;
   
   if(blue >= 0 and blue <= 1)
      this->blue = blue;
   else
      this->blue = 0;
}

//Destructor
RGB::~RGB() {}

//getters - return signed double for compatability with other functions; still guaranteed to be positive
double RGB::getRed() const {return red;}
double RGB::getGreen() const {return green;}
double RGB::getBlue() const {return blue;}

//getters that return the colour values in the standard 256 integer fashion
uint RGB::getIntRed() const {return round(red * MAX_RGB);}
uint RGB::getIntGreen() const {return round(green * MAX_RGB);}
uint RGB::getIntBlue() const {return round(blue * MAX_RGB);}

//setters
void RGB::setRed(double red) {this->red = red;}
void RGB::setGreen(double green) {this->green = green;}
void RGB::setBlue(double blue) {this->blue = blue;}
