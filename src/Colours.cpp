/*
 * Colours.cpp
 *
 *  Created on: Dec 28, 2012
 *      Author: bryan
 */

#include "Colours.h"

namespace colours {

RGB red() {return RGB (1, 0, 0);}
RGB pink() {return RGB (1, 0.5, 0.5);}
RGB maroon() {return RGB (0.5, 0, 0);}

RGB green() {return RGB (0, 1, 0);}
RGB olive() {return RGB (0.5, 0.5, 0);}
RGB darkKhaki() {return RGB (0.7412, 0.7176, 0.4196);}

RGB blue() {return RGB (0, 0, 1);}
RGB lightBlue() {return RGB (0.6784, 0.8471, 0.9020);}
RGB cobaltBlue() {return RGB (0, 0.2784, 0.6706);}
RGB navyBlue() {return RGB (0, 0, 0.5);}

RGB indigo() {return RGB (0.2941, 0, 0.5098);}

RGB yellow() {return RGB (1, 1, 0);}

RGB black() {return RGB (0, 0, 0);}
RGB white() {return RGB (1, 1, 1);}
RGB grey() {return RGB (0.5, 0.5, 0.5);}

QColor qRed() {return rgbToQt(red());}
QColor qPink() {return rgbToQt(pink());}
QColor qMaroon() {return rgbToQt(maroon());}
QColor qGreen() {return rgbToQt(green());}
QColor qOlive() {return rgbToQt(olive());}
QColor qDarkKhaki() {return rgbToQt(darkKhaki());}
QColor qBlue() {return rgbToQt(blue());}
QColor qLightBlue() {return rgbToQt(lightBlue());}
QColor qCobaltBlue() {return rgbToQt(cobaltBlue());}
QColor qNavyBlue() {return rgbToQt(navyBlue());}
QColor qIndigo() {return rgbToQt(indigo());}
QColor qYellow() {return rgbToQt(yellow());}
QColor qBlack() {return rgbToQt(black());}
QColor qWhite() {return rgbToQt(white());}
QColor qGrey() {return rgbToQt(grey());}

//Convert the floating point RGB to floating version of QColor
QColor rgbToQt(const RGB& rgb) {

	QColor colour;
	colour.setRgbF(rgb.getRed(), rgb.getGreen(), rgb.getBlue());
	return colour;
}

} /* namespace colours */
