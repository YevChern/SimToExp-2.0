/*
 * QUtil.cpp
 *
 *  Created on: Nov 26, 2012
 *      Author: bholland
 */

#include <cfloat>

#include <QtGui/QPalette>

#include "Util.h"

namespace qutil {

//Private utility function for calculating a bounding box for a set of points - NOTE: might want to consider moving this to a utility class
QRectF getBoundingRect(const QVector<QPointF>& points) {

	//if nothing, simply return a rect with zero size
	if(points.empty())
		return QRectF(QPointF(0,0), QPointF(0,0));

	double x_min = DBL_MAX, x_max = DBL_MIN;
	double y_min = DBL_MAX, y_max = DBL_MIN;

	int points_size = points.size();
	for(int i = 0; i < points_size; ++i) {

		QPointF cur_point = points[i];

		//do a check for all four
		if(x_min > cur_point.x()) {
			x_min = cur_point.x();
		}

		if(x_max < cur_point.x()) {
			x_max = cur_point.x();
		}

		if(y_min > cur_point.y()) {
			y_min = cur_point.y();
		}

		if(y_max < cur_point.y()) {
			y_max = cur_point.y();
		}
	}

	//return the rectangle bordering all the points
	return QRectF (QPointF(x_min, y_max), QPointF(x_max, y_min));
}

//function for changing the background colour of a widget - specifically it changes the
//colour of the QPallete::Base
void setBackgroundColour(QWidget* widget, QColor colour) {

	QPalette palette = widget->palette();
	palette.setColor(QPalette::Base, colour);
	widget->setPalette(palette);
}

} /* namespace qutil */
