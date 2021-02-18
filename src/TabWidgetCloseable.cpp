/*
 * QTabWidgetCloseable.cpp
 *
 *  Created on: Dec 5, 2012
 *      Author: bholland
 */
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QTabBar>
#else
    #include <QtGui/QTabBar>
#endif

#include "TabWidgetCloseable.h"

//constructor that passes parent to super constructor
TabWidgetCloseable::TabWidgetCloseable(QWidget* parent) : QTabWidget(parent) {}

//destructor
TabWidgetCloseable::~TabWidgetCloseable() {}

//wrapper that hides the close button for tabs that shoulnd't be closed
void TabWidgetCloseable::hideCloseButton(int index) {

   tabBar()->tabButton(index, QTabBar::RightSide)->resize(0,0);
}

//returns the close button that is in the tab bar
QPushButton* TabWidgetCloseable::closeButton(int index) {

	return (QPushButton*) tabBar()->tabButton(index, QTabBar::RightSide);
}
