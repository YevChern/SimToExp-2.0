/*
 * TrajectoryInfoWindow.cpp
 *
 *  Created on: Apr 12, 2013
 *      Author: bholland
 */

//for debugging
#include <iostream>

//Qt includes
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QGridLayout>
    #include <QtWidgets/QHBoxLayout>
    #include <QtWidgets/QFileDialog>
    #include <QtWidgets/QMessageBox>
    #include <QtWidgets/QDesktopWidget>
#else
    #include <QtGui/QGridLayout>
    #include <QtGui/QHBoxLayout>
    #include <QtGui/QFileDialog>
    #include <QtGui/QMessageBox>
    #include <QtGui/QDesktopWidget>
#endif
#include <QtCore/QList>

//project includes
#include "TrajectoryInfoWindow.h"

//make a list of trajectory types to go in drop down menu
void createTrajectoryLists();

QList<TrajectoryType> traj_types;
QList<QString> traj_titles;
QList<QString> traj_opt_paths;

const QString gromacs_title = QString("Gromacs XTC");
const QString gromacs_opt_path = QString("'GRO' path:");

const QString namd_title = QString("CHARMM/NAMD DCD");
const QString namd_opt_path = QString("'PSF' path:");

//Default constructor
TrajectoryInfoWindow::TrajectoryInfoWindow(QWidget* parent) : QWidget(parent) {

	info.setDefaultValues();
	setupWindow();
}

//Destructor
TrajectoryInfoWindow::~TrajectoryInfoWindow() {}

//Getters
TrajectoryInfo TrajectoryInfoWindow::trajectoryInfo() {return info;}
const TrajectoryInfo& TrajectoryInfoWindow::trajectoryInfo() const {return info;}

//private function that does the work of creating the window
void TrajectoryInfoWindow::setupWindow() {

	//first need to set the widget as a window since it has a parent
	this->setWindowFlags(Qt::Window);

	//delete the window upon closing

	createTrajectoryLists();

	//create all of the boxes
	createPathBox();
	createGranularityBox();
	createAxisBox();
	createHeavyFractionBox();
	createNormalBox();
	createPlanarBox();

	//use a grid layout and add all of the boxes to it
	QGridLayout* grid = new QGridLayout(this); //should get deleted when 'this' is deleted (?)

	grid->addWidget(&path_box, 0, 0, 1, -1);
	grid->addWidget(&granularity_box, 1, 0, 1, -1);
	grid->addWidget(&axis_box, 2, 0, 1, -1);
	grid->addWidget(&heavy_fraction_box, 3, 0, 1, -1);
	grid->addWidget(&normal_box, 4, 0, 1, -1);
	grid->addWidget(&planar_box, 5, 0, 1, -1);

	//also add the 'command' buttons
	done_button.setText("Analyze trajectory");
	cancel_button.setText("Cancel");
	connect(&done_button, SIGNAL(clicked()), this, SLOT(dataEntryComplete()));
	connect(&cancel_button, SIGNAL(clicked()), this, SLOT(cancel()));

	QHBoxLayout* hbox = new QHBoxLayout;
	hbox->addWidget(&done_button);
	hbox->addWidget(&cancel_button);
	hbox->addStretch(1);
	grid->addLayout(hbox, 10, 0, 1, -1);

	setWindowTitle(tr("Open and analyze a trajectory file"));
	setAttribute(Qt::WA_DeleteOnClose, true); //don't want it around in memory after
	resize(500, 400);

	//center the window in the screen
	QRect frect = frameGeometry();
	frect.moveCenter(QDesktopWidget().availableGeometry().center());
	move(frect.topLeft());

	//want the focus only on this window, so disable the parent
	setEnabled(true);
}

//private function that sets up the path group box
void TrajectoryInfoWindow::createPathBox() {

	//populate the combo box with the trajectory types available
	int titles_size = traj_titles.size();
	for(int i = 0; i < titles_size; ++i) {
		traj_chooser.insertItem(i, traj_titles[i]);
	}

	//Create a vertical layout for drop down / paths
	QVBoxLayout* vbox = new QVBoxLayout(&path_box);
	QLabel* traj_label = new QLabel(tr("Trajectory type:"), &path_box);
	QHBoxLayout* traj_hbox = new QHBoxLayout;
	traj_hbox->addWidget(traj_label);
	traj_hbox->addWidget(&traj_chooser);
	vbox->addLayout(traj_hbox);

	//these will get deleted when 'path_box' does
	QLabel* path_label = new QLabel(tr("File path:"), &path_box);
	browse_button.setText("Browse");
	browse_button.setObjectName(QString::fromUtf8("browse_button"));

	//Optionally available file path - label will change with trajectory type, but
	//start with first option
	opt_path_label.setText(traj_opt_paths[0]);
	opt_browse_button.setText("Browse");
    opt_browse_button.setObjectName(QString::fromUtf8("opt_browse_button"));

	QHBoxLayout* path_hbox = new QHBoxLayout;
	path_hbox->addWidget(path_label);
	path_hbox->addWidget(&path_line_edit);
	path_hbox->addWidget(&browse_button);
	vbox->addLayout(path_hbox);

	QHBoxLayout* opt_path_hbox = new QHBoxLayout;
	opt_path_hbox->addWidget(&opt_path_label);
	opt_path_hbox->addWidget(&opt_path_line_edit);
	opt_path_hbox->addWidget(&opt_browse_button);
	vbox->addLayout(opt_path_hbox);

	//make necessary connections
	connect(&path_line_edit, SIGNAL(textChanged(QString)), this, SLOT(updateFilePath(QString)));
	connect(&opt_path_line_edit, SIGNAL(textChanged(QString)), this, SLOT(updateOptFilePath(QString)));
	connect(&browse_button, SIGNAL(clicked()), this, SLOT(openFileBrowser()));
	connect(&opt_browse_button, SIGNAL(clicked()), this, SLOT(openFileBrowser()));
	connect(&traj_chooser, SIGNAL(currentIndexChanged(int)), this, SLOT(updateTrajInfo(int)));
}

//private function that sets up the granularity group box
void TrajectoryInfoWindow::createGranularityBox() {

	granularity_box.setTitle(tr("Simulation granularity:"));

	atomistic_button.setText(tr("&Atomistic"));
	cg_button.setText(tr("&MARTINI"));

	atomistic_button.setChecked(true); //assume atomistic

	QHBoxLayout* hbox = new QHBoxLayout(&granularity_box);
	hbox->addWidget(&atomistic_button);
	hbox->addWidget(&cg_button);

	granularity_box.setLayout(hbox);
}

//private function that sets up the axis group box
void TrajectoryInfoWindow::createAxisBox() {

	//all dynamically created objects are children of the window, should be deleted automatically when window
	//goes out of scope or deleted
	axis_box.setTitle(tr("Axis normal to surface:"));

	x_button.setText(tr("&x"));
	y_button.setText(tr("&y"));
	z_button.setText(tr("&z"));

	z_button.setChecked(true); //assume z-axis to start

	QHBoxLayout* hbox = new QHBoxLayout(&axis_box);
	hbox->addWidget(&x_button);
	hbox->addWidget(&y_button);
	hbox->addWidget(&z_button);

	axis_box.setLayout(hbox);
}

//private function that sets up the heavy fraction box
void TrajectoryInfoWindow::createHeavyFractionBox() {

	//create a QString to put in the subscript unicode char for '2'
	QString label = "Fraction of heavy water (D";
	label.append(QChar(0x2082));
	label.append("O):");

	heavy_fraction_label.setText(label);
	heavy_fraction_lineEdit.setText(QString::number(info.heavyWaterFraction()));

	QHBoxLayout* hbox = new QHBoxLayout(&heavy_fraction_box);
	hbox->addWidget(&heavy_fraction_label);
	hbox->addWidget(&heavy_fraction_lineEdit);

	heavy_fraction_box.setLayout(hbox);
}

//private function that sets up the normal group box
void TrajectoryInfoWindow::createNormalBox() {

	normal_box.setTitle(tr("Calculate probabilities normal to surface:"));
	normal_box.setCheckable(true);

	//set text and starting values
	QString grid_label = "Grid size (";
	grid_label.append(QChar(0x00C5)); //unicode symbol for Angstrom
	grid_label.append(")");
	normal_grid_label.setText(grid_label);
	normal_grid_lineEdit.setText(QString::number(info.normalGridSize()));
	normal_avg_button.setText(tr("Ensemble distribution"));
	normal_avg_button.setChecked(true);
	normal_time_button.setText(tr("As function of time"));
	normal_num_label.setText(tr("Number of samples per point"));
//	normal_num_lineEdit.setReadOnly(true); //only make writable if the time_button is checked
	//disable both 'num' variable to begin - only available for 'time_button'
	normal_num_label.setEnabled(false);
	normal_num_lineEdit.setEnabled(false);

	//layout for the box
	QVBoxLayout* vbox = new QVBoxLayout(&normal_box);

	//for grid label and data entry
	QHBoxLayout* grid_box = new QHBoxLayout;
	grid_box->addWidget(&normal_grid_label);
	grid_box->addWidget(&normal_grid_lineEdit);
	grid_box->addStretch(1);

	//add vertically
	vbox->addLayout(grid_box);
	vbox->addWidget(&normal_avg_button);
	vbox->addWidget(&normal_time_button);

	//make a horizontal box for the line edit
	QHBoxLayout* hbox = new QHBoxLayout;
	hbox->addStretch(1);
	hbox->addWidget(&normal_num_label);
	hbox->addWidget(&normal_num_lineEdit);
	vbox->addLayout(hbox);

	//connections for this box
	connect(&normal_time_button, SIGNAL(toggled(bool)), this, SLOT(normalLineEditWritable(bool)));
	connect(&normal_box, SIGNAL(toggled(bool)), this, SLOT(calcNormal(bool)));
}

//private function that sets up the planar group box
void TrajectoryInfoWindow::createPlanarBox() {

	planar_box.setTitle(tr("Calculate probabilities parallel to surface:"));
	planar_box.setCheckable(true);
	planar_box.setChecked(false);

	//set text and starting values
	QString grid_label = "Radial grid size (";
	grid_label.append(QChar(0x00C5)); //unicode symbol for Angstrom
	grid_label.append(")");
	radial_grid_label.setText(grid_label);
	radial_grid_lineEdit.setText(QString::number(info.radialGridSize()));
	angular_grid_label.setText(tr("Angular grid size (rads)"));
	angular_grid_lineEdit.setText(QString::number(info.angularGridSize()));
	planar_avg_button.setText(tr("Ensemble distribution"));
	planar_avg_button.setChecked(true);
	planar_time_button.setText(tr("As function of time"));
	planar_num_label.setText(tr("Number of samples per point"));
	planar_num_label.setEnabled(false);
	planar_num_lineEdit.setEnabled(false); //only make writable if the time_button is checked

	QVBoxLayout* vbox = new QVBoxLayout(&planar_box);

	//for grid label and data entry
	QHBoxLayout* radial_grid_box = new QHBoxLayout;
	radial_grid_box->addWidget(&radial_grid_label);
	radial_grid_box->addWidget(&radial_grid_lineEdit);
	radial_grid_box->addStretch(1);

	//for grid label and data entry
	QHBoxLayout* angular_grid_box = new QHBoxLayout;
	angular_grid_box->addWidget(&angular_grid_label);
	angular_grid_box->addWidget(&angular_grid_lineEdit);
	angular_grid_box->addStretch(1);

	vbox->addLayout(radial_grid_box);
	vbox->addLayout(angular_grid_box);
	vbox->addWidget(&planar_avg_button);
	vbox->addWidget(&planar_time_button);

	//make a horizontal box for the line edit
	QHBoxLayout* hbox = new QHBoxLayout;
	hbox->addStretch(1);
	hbox->addWidget(&planar_num_label);
	hbox->addWidget(&planar_num_lineEdit);
	vbox->addLayout(hbox);

	//connections for this box
	connect(&planar_time_button, SIGNAL(toggled(bool)), this, SLOT(planarLineEditWritable(bool)));
	connect(&planar_box, SIGNAL(toggled(bool)), this, SLOT(calcPlanar(bool)));
}

//private slot that updates the member 'file_path' whenever the line edit is edited
void TrajectoryInfoWindow::updateFilePath(const QString& text) {

	info.setFilePath(text);
}

//private slot that updates the member 'opt_file_path' whenever the line edit changes
void TrajectoryInfoWindow::updateOptFilePath(const QString& text) {

	info.setOptFilePath(text);
}

//private slot that updates the optional path label when the combo box is changed
void TrajectoryInfoWindow::updateTrajInfo(int index) {

	opt_path_label.setText(traj_opt_paths[index]);
	info.setTrajectoryType(traj_types[index]);
}

//private slot that opens the standard file browser dialog
void TrajectoryInfoWindow::openFileBrowser() {

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open trajectory file"), "", "Trajectory files (*.xtc *.gro *.dcd *.psf *.pdb) ;; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//update the text in the line edit with the chosen filename
	//need to figure out which browse button was used
	QObject* sender = QObject::sender();
	if(sender->objectName() == "browse_button")
		path_line_edit.setText(fileName);
	else if(sender->objectName() == "opt_browse_button")
		opt_path_line_edit.setText(fileName);
}

//private slot that closes the window after user information has been entered. It leaves everything else,
//so all of the information stored in the window object is still accessible.
//Emits: dataEntered() signal for main window to know what has happened
void TrajectoryInfoWindow::dataEntryComplete() {

	//set all of the class members so they can be read
	//file paths are updated automatically, so don't need to worry about them

	//first check to make sure at least one box is checked, otherwise just go back
	if(not (normal_box.isChecked() or planar_box.isChecked())) {

		QMessageBox::warning(this, "Input error:", "At least one of the 'normal/parallel' boxes must be set.");
		return;
	}

	//granularity buttons
	if(atomistic_button.isChecked()) {
		info.setSimGranularity(ATOMISTIC);
	}
	else {
		info.setSimGranularity(MARTINI);
	}

	//normal axis buttons
	if(x_button.isChecked()) {
		info.setNormalAxis(X_AXIS);
	}
	else if(y_button.isChecked()) {
		info.setNormalAxis(Y_AXIS);
	}
	else {
		info.setNormalAxis(Z_AXIS);
	}

	//heavy water fraction
	bool fraction_ok;
	double heavy_fraction = heavy_fraction_lineEdit.text().toDouble(&fraction_ok);
	if(not fraction_ok or heavy_fraction < 0 or heavy_fraction > 1) {

		QMessageBox::warning(this, "Input error:", "Fraction of heavy water must be a real value, 0 <= f <= 1.");
		return;
	}

	info.setHeavyWaterFraction(heavy_fraction);

	//then check all grid values, if any of them are invalid AND the box is checked, just return (with msg of course)
	//need to do this here to avoid having to reset a bunch of values below if the planar values or invalid
	bool normal_ok, radial_ok, angular_ok;

	double normal_grid = normal_grid_lineEdit.text().toDouble(&normal_ok);
	double radial_grid = radial_grid_lineEdit.text().toDouble(&radial_ok);
	double angular_grid = angular_grid_lineEdit.text().toDouble(&angular_ok);

	//check for errors
	if((normal_box.isChecked() and (not normal_ok or normal_grid <= 0)) or
	   (planar_box.isChecked() and (not radial_ok or radial_grid <= 0 or
			                        not angular_ok or angular_grid <= 0))) {

		QMessageBox::warning(this, "Input error:", "Grid sizes must be a real values > 0.");
		return;
	}

	//booleans for whether the normal/planar should be calculated are already set
	if(info.calculateNormalToSurface()) {

		//set grid value
		info.setNormalGridSize(normal_grid);

		if(normal_avg_button.isChecked()) {
			info.setNormalTimeAvg(true);
		}
		else if(normal_time_button.isChecked()) {

			info.setNormalFunctionOfTime(true);
			//get number of samples
			bool num_ok;
			info.setNormalSamplesPerPoint(normal_num_lineEdit.text().toUInt(&num_ok));

			if(not num_ok or info.normalSamplesPerPoint() < 1) {

				QMessageBox::warning(this, "Input error:", "Number of samples must be an integer > 0.");
				return;
			}
		}
	}

	if(info.calculateParallelToSurface()) {

		//try setting the grid values, if they don't work, give error msg and return

		info.setRadialGridSize(radial_grid);
		info.setAngularGridSize(angular_grid);

		if(planar_avg_button.isChecked()) {
			info.setPlanarTimeAvg(true);
		}
		else if(planar_time_button.isChecked()) {

			info.setPlanarFunctionOfTime(true);
			//get number of samples
			bool num_ok;
			info.setPlanarSamplesPerPoint(planar_num_lineEdit.text().toUInt(&num_ok));

			if(not num_ok or info.planarSamplesPerPoint() < 1) {

				QMessageBox::warning(this, "Input error:", "Number of samples must be an integer > 0.");
				return;
			}
		}
	}

	emit dataEntered(this); //pass pointer to this object so it can be used then deleted
}

//private slot that closes the window and clears all the window data.
//Emits: canceled() signal for main window to know to abort
void TrajectoryInfoWindow::cancel() {

	info.setDefaultValues(); //just in case, the object should go out of scope anyway
	emit canceled(this); //pass pointer to object so it can be deleted
}

//private slot that turns the read only state of the line edit on or off
//also clears any data if turning to read only
void TrajectoryInfoWindow::normalLineEditWritable(bool writable) {

	if(writable) {

		normal_num_label.setEnabled(true);
		normal_num_lineEdit.setEnabled(true);
	}
	else {

		normal_num_label.setEnabled(false);
		normal_num_lineEdit.clear();
		normal_num_lineEdit.setEnabled(false);
	}
}

//private slot that turns the read only state of the line edit on or off
//also clears any data if turning to read only
void TrajectoryInfoWindow::planarLineEditWritable(bool writable) {

	if(writable) {

		planar_num_label.setEnabled(true);
		planar_num_lineEdit.setEnabled(true);
	}
	else {

		planar_num_label.setEnabled(false);
		planar_num_lineEdit.clear();
		planar_num_lineEdit.setEnabled(false);
	}
}

//private slot that sets the boolean for whether to calculate the normal direction
void TrajectoryInfoWindow::calcNormal(bool normal_checked) {

	if(normal_checked)
		info.setCalcNormalToSurface(true);
	else
		info.setCalcNormalToSurface(false);
}

//private slot that sets the boolean for whether to calculate the normal direction
void TrajectoryInfoWindow::calcPlanar(bool planar_checked) {

	if(planar_checked)
		info.setCalcParallelToSurface(true);
	else
		info.setCalcParallelToSurface(false);
}

//Virtual slot from QWidget overwritten to make the 'close' button in the upper corner of the window do the
//same as cancel
void TrajectoryInfoWindow::closeEvent(QCloseEvent* event) {

	cancel();
	event->accept();
}

//Just assembles the list of trajectory options
void createTrajectoryLists() {

	traj_types.append(GROMACS_XTC);
	traj_types.append(NAMD_DCD);

	traj_titles.append(gromacs_title);
	traj_titles.append(namd_title);

	traj_opt_paths.append(gromacs_opt_path);
	traj_opt_paths.append(namd_opt_path);
}
