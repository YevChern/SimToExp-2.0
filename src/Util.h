/*
 * QUtil.h
 *
 *  Created on: Nov 26, 2012
 *      Author: bholland
 *
 *  Class that provides utility functions for Qt based code
 */

#ifndef QUTIL_H_
#define QUTIL_H_

#include <QtCore/QRect>
#include <QtCore/QVector>
#include <QtCore/QPointF>
#include <QtGui/QColor>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QWidget>
#else
    #include <QtGui/QWidget>
#endif

namespace qutil {

	QRectF getBoundingRect(const QVector<QPointF>&);
	void setBackgroundColour(QWidget*, QColor);

} /* namespace qutil */
#endif /* QUTIL_H_ */
