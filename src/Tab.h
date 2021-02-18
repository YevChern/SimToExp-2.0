/*
 * QTab.h
 *
 *  Created on: Dec 7, 2012
 *      Author: bholland
 *
 * Simple class that inherits from QWidget (which is normally used for tabs), but I want
 * to include a few extra slots/functions.
 *
 * Specifically, the tab provides the ability to send a signal to the QTabWidget that
 * it wants to close, and includes the index.
 *
 * The enum types here as specific to SIMtoEXP but could be changed for any situation
 */

#ifndef QTAB_H_
#define QTAB_H_

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QWidget>
#else
    #include <QtGui/QWidget>
#endif

//describes what the tab is used for
enum TabName {Scattering, AtomicFF, SimSpecifics, Unknown};

class Tab : public QWidget {

	Q_OBJECT;

	int tab_index;
	TabName name;

public:

	Tab(QWidget* parent = 0);
	Tab(TabName, QWidget* parent = 0);
	virtual ~Tab();

	void setTabIndex(int);
	int getTabIndex() const;
	TabName getTabName() const;

public slots:

	void closeTab();

signals:

	void needToClose(int);
};

#endif /* QTAB_H_ */
