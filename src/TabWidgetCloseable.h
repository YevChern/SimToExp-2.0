/*
 * QTabWidgetCloseable.h
 *
 *  Created on: Dec 5, 2012
 *      Author: bholland
 */

#ifndef QTABWIDGETCLOSEABLE_H_
#define QTABWIDGETCLOSEABLE_H_

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QTabWidget>
    #include <QtWidgets/QPushButton>
#else
    #include <QtGui/QTabWidget>
    #include <QtGui/QPushButton>
#endif


class TabWidgetCloseable : public QTabWidget {

public:

	TabWidgetCloseable(QWidget* parent = 0);
	virtual ~TabWidgetCloseable();

	void hideCloseButton(int);

	QPushButton* closeButton(int);
};

#endif /* QTABWIDGETCLOSEABLE_H_ */
