/*
 * Colours.h
 *
 *  Created on: Dec 28, 2012
 *      Author: Bryan Holland
 *
 * Namespace for preset colours; outputs both QColor and RGB objects, the former
 * being very useful when using Qt - some of these colours are available in the Qt presets,
 * but it is a very limited set of colours.
 */

#ifndef COLOURS_H_
#define COLOURS_H_

#include <QtGui/QColor>

#include "RGB.h"

namespace colours {

	//reds
	RGB red();
	RGB pink();
	RGB maroon();

	//greens
	RGB green();
	RGB olive();
	RGB darkKhaki();

	//blues
	RGB blue();
	RGB lightBlue();
	RGB cobaltBlue();
	RGB navyBlue();

	//purples
	RGB indigo();

	//yellows
	RGB yellow();
	RGB black();
	RGB white();
	RGB grey();

	//Qt versions
	QColor qRed();
	QColor qPink();
	QColor qMaroon();
	QColor qGreen();
	QColor qOlive();
	QColor qDarkKhaki();
	QColor qBlue();
	QColor qLightBlue();
	QColor qCobaltBlue();
	QColor qNavyBlue();
	QColor qIndigo();
	QColor qYellow();
	QColor qBlack();
	QColor qWhite();
	QColor qGrey();

	QColor rgbToQt(const RGB&);

} /* namespace colours */
#endif /* COLOURS_H_ */
