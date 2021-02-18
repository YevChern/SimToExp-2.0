/*
 * QTab.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: bholland
 */

#include "Tab.h"

//Default constructor; set index to -1 meaning no index has been set
Tab::Tab(QWidget* parent) : QWidget(parent) {

	tab_index = -1;
	name = Unknown;
}

//Constructor that also sets the tab name/type
Tab::Tab(TabName name, QWidget* parent) : QWidget(parent) {

	tab_index = -1;
	this->name = name;
}

//Destructor
Tab::~Tab() {}

//getter
int Tab::getTabIndex() const {return tab_index;}
TabName Tab::getTabName() const {return name;}

//setter
void Tab::setTabIndex(int index) {tab_index = index;}

//SLOT - for receiving a closing signal, just emits its own signal with index information
void Tab::closeTab() {

	emit needToClose(tab_index);
}
