/*
 *  RGB.h
 *
 *  Created on: Oct 12, 2012
 *      Author: bholland
 * 
 *  Simple class that defines an RGB colour using a range between 0 and 1
 *  for red, green and blue. 
 *
 */ 

#ifndef RGB_H_
#define RGB_H_ 

typedef unsigned int uint;

class RGB {
   
   double red;
   double green;
   double blue;
   
   public:
      
      RGB();
      RGB(double,double,double);
      virtual ~RGB();
   
      double getRed() const;
      double getGreen() const;
      double getBlue() const;
      
      uint getIntRed() const;
      uint getIntGreen() const;
      uint getIntBlue() const;
      
      void setRed(double);
      void setGreen(double);
      void setBlue(double);
};

#endif /* RGB_H_ */
