/*
 * ScatteringTabFrame.h
 *
 *  Created on: Dec 7, 2012
 *      Author: bholland
 *
 * Contains all of the layout and data/info for the Scattering Density Tab that opens from
 * Tools -> Scattering Lengths
 */

#ifndef SCATTERINGTABFRAME_H_
#define SCATTERINGTABFRAME_H_
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QWidget>
    #include <QtWidgets/QFrame>
    #include <QtWidgets/QHBoxLayout>
    #include <QtWidgets/QGridLayout>
    #include <QtWidgets/QPushButton>
    #include <QtWidgets/QLabel>
    #include <QtWidgets/QLineEdit>
#else
    #include <QtGui/QWidget>
    #include <QtGui/QFrame>
    #include <QtGui/QHBoxLayout>
    #include <QtGui/QGridLayout>
    #include <QtGui/QPushButton>
    #include <QtGui/QLabel>
    #include <QtGui/QLineEdit>
#endif

#include <QtCore/QList>
#include <QtCore/QMap>
#include <QtGui/QColor>
#include <QtCore/QSet>

#include "AtomicInfo.h"

typedef QList<QLabel*> Labels;
typedef QMap<AtomType,QLineEdit*> LineEditMap;
typedef QSet<QLineEdit*> LineEditSet;

class ScatteringTabFrame : public QFrame {

	Q_OBJECT

	AtomicInfo* info; //object that stores all the particle specifics

	//the layouts used
	QGridLayout main_layout;
	QVBoxLayout button_layout;
	QGridLayout info_layout;
	QVBoxLayout text_layout;

	//buttons
	QPushButton apply_changes_button;
	QPushButton reset_defaults_button;

	//containers for data output / editing
	Labels labels;
	LineEditMap name_edits;
	LineEditMap id_edits;
	LineEditMap nelec_edits;
	LineEditMap nsl_edits;
	Labels messages;

	//For tracking edited containers
	LineEditSet edited_names;
	LineEditSet edited_ids;
	LineEditSet edited_nelecs;
	LineEditSet edited_nsls;


public:

	ScatteringTabFrame(AtomicInfo*, QWidget* parent = 0);
	virtual ~ScatteringTabFrame();

	QGridLayout* mainLayout();

//slots here are private since no the functions are only used internally
private slots:

	void applyChanges();
	void resetDefaults();
	void changeWidgetColour();
	void nameEdited();
	void idEdited();
	void nelecEdited();
	void nslEdited();

private:

	void createContainters();
	void createLayout();
	void createConnections();
};

#endif /* SCATTERINGTABFRAME_H_ */
