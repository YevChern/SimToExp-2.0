/*
 * ScatteringTabFrame.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: bholland
 */

#include <vector>

#include "ScatteringTabFrame.h"
#include "AtomType.h"
#include "Util.h"
#include "Colours.h"

#include <QtCore/QtDebug>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QMessageBox>
#endif

using namespace std;

typedef QList<QLineEdit*> LineEditList;

const QString NAME_STRING = "Particle name";
const QString ID_STRING = "Particle ID";
const QString NELEC_STRING = "Electron number";
const QString NSL_STRING = "Neutron scattering length";

const int LABELS_ROW = 1;
const int NAMES_COLUMN = 0;
const int ID_COLUMN = 1;
const int NELEC_COLUMN = 2;
const int NSL_COLUMN = 3;

const int BUTTONS_COLUMN = 0;

const QString CHANGES_MESSAGE = "Once changes have been made, press the 'Apply changes' button to finalize; accepted values will have a green background. Changes applied will only alter data loaded henceforth.  To cancel changes without applying them, simply close the tab.";
const QString DEFAULT_MESSAGE = "To return to the default scattering values, simply press the 'Reset defaults' button. The default values will NOT be applied to the currently loaded data.  Also, this will NOT affect user defined particles.";
const QString VALUES_REFERENCE = "The default values are taken from: ";

const int NEUTRON_SL_PRECISION = 3;
const int MESSAGE_MARGIN = 15; //width margin around the text in pixels
const int MESSAGE_WIDTH = 230; //since can't seem to get it to fill properly, set a max width
const int BUTTON_HEIGHT = 50;
const int BUTTON_TOP_SPACE = 30;
const int BUTTON_SPACE_BETWEEN = 25;
const int BUTTON_INFO_SPACING = 15;
const int BLANK_ROW_HEIGHT = 20;

//Only constructor as the frame requires the AtomicInfo pointer to function properly
ScatteringTabFrame::ScatteringTabFrame(AtomicInfo* info, QWidget* parent) : QFrame(parent) {

	this->info = info;
	createContainters();
	createLayout();
	createConnections();
}

//Destructor
ScatteringTabFrame::~ScatteringTabFrame() {

	//need to delete the container contents
	qDeleteAll(labels);
	qDeleteAll(name_edits);
	qDeleteAll(id_edits);
	qDeleteAll(nelec_edits);
	qDeleteAll(nsl_edits);
	qDeleteAll(messages);
}

//getter for layout
QGridLayout* ScatteringTabFrame::mainLayout() {return &main_layout;}

//private function that sets up the layout to be presented in the tab
void ScatteringTabFrame::createLayout() {

	//put a spacer for the first row
	main_layout.setRowMinimumHeight(0, BLANK_ROW_HEIGHT);

	//start with buttons, go in their own vertical layout
	apply_changes_button.setText("Apply changes");
	apply_changes_button.setFixedHeight(BUTTON_HEIGHT);

	reset_defaults_button.setText("Reset to default values");
	reset_defaults_button.setFixedHeight(BUTTON_HEIGHT);

	button_layout.addSpacing(BUTTON_TOP_SPACE);
	button_layout.addWidget(&apply_changes_button);
	button_layout.addSpacing(BUTTON_SPACE_BETWEEN);
	button_layout.addWidget(&reset_defaults_button);
	button_layout.addStretch(1);

	main_layout.addLayout(&button_layout, 1, BUTTONS_COLUMN);
	main_layout.setColumnMinimumWidth(BUTTONS_COLUMN+1, BUTTON_INFO_SPACING);

	//Now the info go in a grid layout
	//The labels
	for(int i = 0; i < labels.size(); i++)
		info_layout.addWidget(labels[i], LABELS_ROW, i);

	//The line edits
	int row = LABELS_ROW + 1;
	LineEditMap::iterator itr;
	for(itr = name_edits.begin(); itr != name_edits.end(); itr++) {
		info_layout.addWidget(itr.value(), row, NAMES_COLUMN);
		row++;
	}

	//reset row and do other line edits
	row = LABELS_ROW + 1;
	for(itr = id_edits.begin(); itr != id_edits.end(); itr++) {
		info_layout.addWidget(itr.value(), row, ID_COLUMN);
		row++;
	}

	row = LABELS_ROW + 1;
	for(itr = nelec_edits.begin(); itr != nelec_edits.end(); itr++) {
		info_layout.addWidget(itr.value(), row, NELEC_COLUMN);
		row++;
	}

	row = LABELS_ROW + 1;
	for(itr = nsl_edits.begin(); itr != nsl_edits.end(); itr++) {
		info_layout.addWidget(itr.value(), row, NSL_COLUMN);
		row++;
	}

	main_layout.addLayout(&info_layout, 1, BUTTONS_COLUMN+2);

	//Add the text on the far right
	for(int i = 0; i < messages.size(); i++) {

		messages[i]->setWordWrap(true);
		messages[i]->setMargin(MESSAGE_MARGIN);
		messages[i]->setMaximumWidth(MESSAGE_WIDTH);
		text_layout.addWidget(messages[i]);
	}
	text_layout.addStretch(1);

	main_layout.addLayout(&text_layout, 1, BUTTONS_COLUMN+3);
	this->setLayout(&main_layout);
}

//SLOT - changes the background colour of the widget sending the signal
void ScatteringTabFrame::changeWidgetColour() {

	//assume the object is a  - this will break if it isn't but the object should
	//be a widget to send a signal, so should be okay
	QWidget* sender = (QWidget*) QObject::sender();
	qutil::setBackgroundColour(sender, colours::qLightBlue());
}

//SLOT - simply adds the calling widget to a list of edited QLineEdits
void ScatteringTabFrame::nameEdited() {

	QLineEdit* sender = (QLineEdit*) QObject::sender();
	edited_names.insert(sender);
}

//SLOT - simply adds the calling widget to a list of edited QLineEdits
void ScatteringTabFrame::idEdited() {

	QLineEdit* sender = (QLineEdit*) QObject::sender();
	edited_ids.insert(sender);
}

//SLOT - simply adds the calling widget to a list of edited QLineEdits
void ScatteringTabFrame::nelecEdited() {

	QLineEdit* sender = (QLineEdit*) QObject::sender();
	edited_nelecs.insert(sender);
}

//SLOT - simply adds the calling widget to a list of edited QLineEdits
void ScatteringTabFrame::nslEdited() {

	QLineEdit* sender = (QLineEdit*) QObject::sender();
	edited_nsls.insert(sender);
}

//SLOT - for connecting to button for applying changes. Actually applies the changes
//to the AtomicInfo object
void ScatteringTabFrame::applyChanges() {

	bool elec_err_msg_displayed = false;
	bool nsl_err_msg_displayed = false;

	//create a temporary container for each type of the edited widgets, and clear
	//the class members.  This is done so that at any point if a change fails, the pointer
	//can be added back to the appropriate container for next time
	LineEditList names_tmp = edited_names.values();
	edited_names.clear(); //hopefully doesn't 'delete' the objects

	LineEditList ids_tmp = edited_ids.values();
	edited_ids.clear(); //hopefully doesn't 'delete' the objects

	LineEditList nelecs_tmp = edited_nelecs.values();
	edited_nelecs.clear(); //hopefully doesn't 'delete' the objects

	LineEditList nsls_tmp = edited_nsls.values();
	edited_nsls.clear(); //hopefully doesn't 'delete' the objects

	//just go through lists of edited pointers, if empty do nothing
	for(int i = 0; i < names_tmp.size(); i++) {

		//get the AtomType
		AtomType type = name_edits.key(names_tmp[i]);

		//Use this to change the value in AtomicInfo
		string name = names_tmp[i]->text().toStdString();
		info->setUserDefinedName(type, name);
		qutil::setBackgroundColour(names_tmp[i], colours::qDarkKhaki());
	}

	//Need to do some checks for ids, since they must be unique
	for(int i = 0; i < ids_tmp.size(); i++) {

		//get the AtomType
		AtomType type = id_edits.key(ids_tmp[i]);

		//Use this to change the value in AtomicInfo
		QString id_string = ids_tmp[i]->text();

		//check if already defined and return an error message if so
		if(not info->setUserDefinedID(type, id_string.toStdString())) {

			QString msg = "Particle ID '";
			msg += id_string;
			msg += "' is already in use.\nNote: only one of upper/lower letters case can be used.";
			QMessageBox::warning(this, "Data error:", msg);
			edited_ids.insert(ids_tmp[i]); //add back

		} else {

			//if changed, change background
			qutil::setBackgroundColour(ids_tmp[i], colours::qDarkKhaki());
		}

	}

	//Number of electrons must be > 0
	for(int i = 0; i < nelecs_tmp.size(); i++) {

		//first check the validity of the input
		bool ok;
		double nelec = nelecs_tmp[i]->text().toDouble(&ok);
		if(not ok or nelec < 0) {

			//only display the message once for a set of errors
			if(not elec_err_msg_displayed) {

				QString msg = "The number of electrons must be a real number >= 0.";
				QMessageBox::warning(this, "Input error:", msg);
				elec_err_msg_displayed = true;
			}
			edited_nelecs.insert(nelecs_tmp[i]); //add back

		} else {

			//if ok, change the value
			AtomType type = nelec_edits.key(nelecs_tmp[i]);
			info->setNumberElectrons(type, nelec);
			qutil::setBackgroundColour(nelecs_tmp[i], colours::qDarkKhaki());
		}
	}

	//The neutron scattering length can be any value
	for(int i = 0; i < nsls_tmp.size(); i++) {

		bool ok;
		double nsl = nsls_tmp[i]->text().toDouble(&ok);
		if(not ok) {

			//only display error message once for any number of errors
			if(not nsl_err_msg_displayed) {

				QString msg = "The neutron scattering length must be a real number.";
				QMessageBox::warning(this, "Input error:", msg);
				nsl_err_msg_displayed = true;
			}
			edited_nsls.insert(nsls_tmp[i]); //add back

		} else {

			//if ok, change the value
			AtomType type = nsl_edits.key(nsls_tmp[i]);
			info->setNeutronSL(type, nsl);
			qutil::setBackgroundColour(nsls_tmp[i], colours::qDarkKhaki());
		}
	}
}

//SLOT - for connecting to button for resetting the values to the defaults
void ScatteringTabFrame::resetDefaults() {

	//First need to reset the values in the AtomicInfo pointer
	info->resetScatteringDefaults();

	//Then need to re-populate the QLineEdits
	vector<AtomType> types = info->allAtomTypes();
	for(uint i = 0; i < types.size(); i++) {

		//nelecs
		QLineEdit* cur_edit = nelec_edits[types[i]];
		if(types[i] < USER_1)
			qutil::setBackgroundColour(cur_edit, colours::qWhite());

		cur_edit->setText(QString::number(info->numberElectrons(types[i])));

		//nsls
		cur_edit = nsl_edits[types[i]];
		if(types[i] < USER_1)
			qutil::setBackgroundColour(cur_edit, colours::qWhite());

		cur_edit->setText(QString::number(info->neutronSL(types[i]), 'e', NEUTRON_SL_PRECISION));
	}
}

//private function that uses the AtomicInfo to create the line edits and set their text
void ScatteringTabFrame::createContainters() {

	//first make the labels
	labels.append(new QLabel(NAME_STRING));
	labels.append(new QLabel(ID_STRING));
	labels.append(new QLabel(NELEC_STRING));
	labels.append(new QLabel(NSL_STRING));

	//create the line edits
	vector<AtomType> types = info->allAtomTypes();
	for(uint i = 0; i < types.size(); i++) {

		//create name with info, and make it
		QLineEdit* cur_name = new QLineEdit;
		QString name_string(info->nameAndSymbol(types[i]).first.c_str());
		cur_name->setText(name_string);
		name_edits.insert(types[i], cur_name);
		//connect all line edits to a colour changing slot if they are edited at all
		connect(cur_name, SIGNAL(textEdited(QString)), this, SLOT(changeWidgetColour()));
		connect(cur_name, SIGNAL(textEdited(QString)), this, SLOT(nameEdited()));

		//get the id
		QLineEdit* cur_id = new QLineEdit;
		cur_id->setText(QString(info->nameAndSymbol(types[i]).second.c_str()));
		id_edits.insert(types[i], cur_id);
		connect(cur_id, SIGNAL(textEdited(QString)), this, SLOT(changeWidgetColour()));
		connect(cur_id, SIGNAL(textEdited(QString)), this, SLOT(idEdited()));

		QLineEdit* cur_nelec = new QLineEdit;
		cur_nelec->setText(QString::number(info->numberElectrons(types[i])));
		nelec_edits.insert(types[i], cur_nelec);
		connect(cur_nelec, SIGNAL(textEdited(QString)), this, SLOT(changeWidgetColour()));
		connect(cur_nelec, SIGNAL(textEdited(QString)), this, SLOT(nelecEdited()));

		QLineEdit* cur_nsl = new QLineEdit;
		cur_nsl->setText(QString::number(info->neutronSL(types[i]), 'e', NEUTRON_SL_PRECISION));
		nsl_edits.insert(types[i], cur_nsl);
		connect(cur_nsl, SIGNAL(textEdited(QString)), this, SLOT(changeWidgetColour()));
		connect(cur_nsl, SIGNAL(textEdited(QString)), this, SLOT(nslEdited()));

		//make read only if not user defined (i.e. can change the name for 'user defined' particles
		if(types[i] < USER_1) {

			cur_name->setReadOnly(true);
			qutil::setBackgroundColour(cur_name, QColor(Qt::lightGray));
			cur_id->setReadOnly(true);
			qutil::setBackgroundColour(cur_id, QColor(Qt::lightGray));
		}
	}

	//add the message labels
	messages.append(new QLabel(CHANGES_MESSAGE));
	messages.append(new QLabel(DEFAULT_MESSAGE));
	messages.append(new QLabel(VALUES_REFERENCE));
}

//private functions for dealing with any connections that need to be made upon construction
void ScatteringTabFrame::createConnections() {

	connect(&apply_changes_button, SIGNAL(clicked(bool)), this, SLOT(applyChanges()));
	connect(&reset_defaults_button, SIGNAL(clicked(bool)), this, SLOT(resetDefaults()));
}
