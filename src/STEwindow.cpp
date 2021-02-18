/*
 * STEwindow.cpp
 *
 *  Created on: Oct 1, 2012
 *      Author: bholland
 */

//for debugging only
#include <iostream>

//project includes
#include "STEwindow.h"
#include "KillBox.h"
#include "AtomicInfo.h"
#include "ScatteringTabFrame.h"
#include "Util.h"
#include "TrajectoryData.h"

//Qt includes
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
    #include <QtWidgets/QFileDialog>
    #include <QtWidgets/QScrollBar>
    #include <QtWidgets/QScrollArea>
    #include <QtWidgets/QDesktopWidget>
#else
    #include <QtGui/QMessageBox>
    #include <QtGui/QFileDialog>
    #include <QtGui/QScrollBar>
    #include <QtGui/QScrollArea>
    #include <QtGui/QDesktopWidget>
#endif

#include <QtCore/QFile>
#include <QtCore/QtAlgorithms>
#include <QtGui/QPen>
#include <QtGui/QPalette>


//C++ include
#include <cstdlib>
#include <cmath>
#include <cfloat>

using namespace std;
using namespace qutil;

const uint PLOT_FREQUENCY = 5;
const uint MAX_PLOTS = 20;
const int MIN_POINTS = 3;
const uint SYMBOL_WIDTH = 6; //in pixels i think
const int MAX_CURVES = 5;
const uint MAX_DIGITS = 5;
const uint DISPLAY_PRECISION = 3;
const uint NUMBER_OF_Q_LINE_EDITS = 6;

const int NUMBER_OF_DEFAULT_PLOTS = 6;

//plot margins - use percentages of the current domain and range
const double X_AXIS_MARGIN_FRACTION = 0;
const double Y_AXIS_MARGIN_FRACTION = 0.05;

//percentage for zoom menu actions
const double X_AXIS_ZOOM_MULTIPLIER = 2;
const double Y_AXIS_ZOOM_MULTIPLIER = 2;

//Fourier default values
const double XRAY_QMIN = 0.0;
const double XRAY_QMAX = 1.0;
const double XRAY_QSTEP = 0.001;
const double NEUT_QMIN = 0.0;
const double NEUT_QMAX = 0.4;
const double NEUT_QSTEP = 0.0025;
const double WATER_ED = 0.333;
const double WATER_NSLD = 6.38E-6;

STEwindow::STEwindow(QWidget* parent) : QMainWindow(parent) {

	//make this the parent for ExpParser
	exp.setParent(this);

	sim_open = false;
	sim_from_trajectory = false;
	form_factors_done = false;
	gro_file_open = false;
	scattering_tab_open = false;
	tab_area_maximized = false;
	converter = new QLocale(QLocale::English, QLocale::Canada);

	buildColourLibrary();
	buildSymbolLibrary();

	//for writing output / saving data
	out_writer.setParent(this);

	//initialize iterators that keep track of plot colour info
	xray_ff_colour_itr = plot_colours.begin();
	neutron_ff_colour_itr = plot_colours.begin();
	xray_ff_symbol_itr = plot_symbols.begin();
	neutron_ff_symbol_itr = plot_symbols.begin();

	//this passes the pointer through
	sim.setAtomInfo(info);

	//does most of the work that is in 'gui.h'
	setupUi(this);

	//make all the legends checkable, so the curves can be toggled on/off
    e_densities_legend.setDefaultItemMode(QwtLegendData::Checkable);
    neutron_SL_legend.setDefaultItemMode(QwtLegendData::Checkable);
    xray_formfactor_legend.setDefaultItemMode(QwtLegendData::Checkable);
    neutron_formfactor_legend.setDefaultItemMode(QwtLegendData::Checkable);
    num_densities_legend.setDefaultItemMode(QwtLegendData::Checkable);
    volume_legend.setDefaultItemMode(QwtLegendData::Checkable);

	//add the legends to the appropriate plots
	e_densities_plot->insertLegend(&e_densities_legend);
	neutron_SL_plot->insertLegend(&neutron_SL_legend);
	xray_formfactors_plot->insertLegend(&xray_formfactor_legend);
	neutron_formfactors_plot->insertLegend(&neutron_formfactor_legend);
	num_densities_plot->insertLegend(&num_densities_legend);
	vol_probs_plot->insertLegend(&volume_legend);

	//populate lineEdit containers and deletion containers
	buildLineEditContainers();
	buildDeleteActionsContainers();
	buildPlotContainers();
	buildScrollZoomerContainer();
	
	//need to initialize the margin containers
	initializeBoundingRects();

	//make any static connections for the window
	makeConnections();

	//center the window in the screen
	QRect frect = frameGeometry();
	frect.moveCenter(QDesktopWidget().availableGeometry().center());
	move(frect.topLeft());

	//For debugging purposes
	//connect(lineEdit_29, SIGNAL(selectionChanged()), this, SLOT(testButton()));

}

STEwindow::~STEwindow() {

	delete converter;
}

/*
 * Reads the simulation data from file and stores it in SimParser.  Does NOT
 * display, but triggers the signal to display if read in properly.
 */
void STEwindow::openSim() {

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open SIMtoEXP file:"), "", "SIMtoEXP: *.sim, *.SIM (*.sim, *.SIM);; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//returns null if 'cancel' is hit, so only continue if there is a path
	if(not fileName.isEmpty()) {

		//if there is already data open, need to delete first
		if(sim_open)
			removeSim();

		//trigger simOpened only if read was successful
		if(sim.readSimFile(fileName) == READ_SUCCESS) {
			emit simOpened();
		}
	}
}

/*
 * Displays simulation data in the simulation tab
 */
void STEwindow::displaySim() {

	//qDebug() << "displaySim() started from simOpened signal";

	//make the tabs usable
	showTab(num_densities_tab);
	showTab(e_densities_tab);
	showTab(neutron_SL_tab);

	// get curves - these are copies, but still pretty fast
	NormalDensities pds = sim.getElementalData();

	//check to make sure there is enough to plot!
	//qDebug() << "Number of particles: " << pds.size();
	if(pds[0].densities().size() >= MIN_POINTS) {

		//first deal with number densities
		uint size = pds.size();
		for(uint i = 0; i < size; i++) {

			//include the residue name if it exists
			QString curve_name, pName = pds[i].particleName();

			//only keep the string end up until the underscore
			//this gets rid of the MARTINI residue addition to the particle name
			QStringList split_name = pName.split('_');
			pName = split_name.last();

			if(not pds[i].residueName().isEmpty()) {

				curve_name = pName + ", " + pds[i].residueName();
			}
			else
				curve_name = pName;

			QwtPlotCurve* curve = new QwtPlotCurve(curve_name);

			curve->setSamples(pds[i].densities());
			curve->attach(num_densities_plot);
			curve->setPen(QPen (getColour(i*47))); //Make new QPen and give it a random colour, seed with a prime multiple
			
			//put curves into containers to keep track of them (data pointers can be extracted from the curve)
			sim_number_curves.append(curve);

            const QVariant itemInfo = num_densities_plot->itemToInfo(curve);
            QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( num_densities_legend.legendWidget( itemInfo ) );

            //only initially show either a maximum number of plots, or some frequency of them
			if(size <= MAX_PLOTS or i % PLOT_FREQUENCY == 0) {

                item->setChecked(true);
			}
			else {

				curve->hide();
                item->setChecked(false);
			}

		}

		//now deal with electron and neutron SL densities
		const AtomTypes& types = sim.getAtomTypes();
		uint num_atom_types = types.size();

		for(uint i = 0; i < num_atom_types; i++) {

			//get name to pass to plot curve
			QString atom_name(info.nameAndSymbol(types[i]).first.c_str());
			QwtPlotCurve* curve = new QwtPlotCurve(atom_name);
			curve->setSamples(sim.getTotalElectronDensity(types[i]));

			//set colour - use colour associated with the AtomType
			HSVColour cur_hsv = info.colour(types[i]);
			QColor cur_colour;
			cur_colour.setHsvF(cur_hsv.hue(), cur_hsv.saturation(), cur_hsv.value());

			curve->setPen(QPen (cur_colour));
			//attach curve to the plot - MUST DO THIS BEFORE CHECKING LEGEND ITEM, not sure why
			curve->attach(e_densities_plot);

            QVariant itemInfo = e_densities_plot->itemToInfo(curve);
            QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( e_densities_legend.legendWidget( itemInfo ) );

			item->setChecked(true); //make all of them visible


			//put curves into containers to keep track of them (data pointers can be extracted from the curve)
			sim_ed_curves.append(curve);

			curve = new QwtPlotCurve(atom_name);
			curve->setSamples(sim.getTotalNeutronDensity(types[i]));
			curve->attach(neutron_SL_plot);

			//use same colour for neutron SL plot
			curve->setPen(QPen (cur_colour));

            itemInfo = neutron_SL_plot->itemToInfo(curve);
            item = qobject_cast<QwtLegendLabel *>( neutron_SL_legend.legendWidget( itemInfo ) );

            item->setChecked(true); //make all of them visible

			//put curves into containers to keep track of them (data pointers can be extracted from the curve)
			sim_nsld_curves.append(curve);
		}

		//update the three plots
		updatePlot(tabIndex(num_densities_tab));
		updatePlot(tabIndex(e_densities_tab));
		updatePlot(tabIndex(neutron_SL_tab)); //THIS IS MESSED UP FOR SOME REASON

		//now that the bounds are up to date, set the margins for display
		QString left_margin_t = converter->toString(bounding_rects[tabIndex(num_densities_tab)].left());
		QString right_margin_t = converter->toString(bounding_rects[tabIndex(num_densities_tab)].right());

		//update the margin widgets
		left_margin_lineEdit->setText(left_margin_t);
		right_margin_lineEdit->setText(right_margin_t);

		//finally, focus on the number densities plot
		tabWidget->setCurrentIndex(tabIndex(num_densities_tab));

	} else {

		//outpout an error box and do not display anything
		QString message = "There must be >= ";
		message += (converter->toString(MIN_POINTS));
		message += " values to plot.";
		QMessageBox::warning(this, "Plotting error:", message);
		return;
	}
}

/*
 * Removes the simulation data; must be done prior to opening new file
 */
void STEwindow::removeSim() {

	//if no sim file open, simply return
	if(not sim_open)
		return;

	//go through all the curves, detach them and delete all their data
	int size = sim_number_curves.size();
	for(int i = 0; i < size; i++) {

		//detach them from the plot
		sim_number_curves[i]->detach();
    }

	size = sim_ed_curves.size();
	for(int i = 0; i < size; i++) {

		//detach them from the plot
		sim_ed_curves[i]->detach();
    }

	size = sim_nsld_curves.size();
	for(int i = 0; i < size; i++) {

		//detach them from the plot
		sim_nsld_curves[i]->detach();
	}

	qDeleteAll(sim_number_curves);
	qDeleteAll(sim_ed_curves);
	qDeleteAll(sim_nsld_curves);
	sim_number_curves.clear();
	sim_ed_curves.clear();
	sim_nsld_curves.clear();

	//do the same for any Form factor curves kicking around
	cleanFourierCurves();

	//remove data from SimParser
	sim.removeData();

	//update all the applicable plots
	updatePlot(tabIndex(num_densities_tab));
	updatePlot(tabIndex(e_densities_tab));
	updatePlot(tabIndex(neutron_SL_tab));
	updatePlot(tabIndex(xray_formfactors_tab));
	updatePlot(tabIndex(neutron_formfactors_tab));

	//disable the tabs if nothing else is using them
	hideTab(num_densities_tab); //no experimental data
	if(exp_ed_curves.empty())
		hideTab(e_densities_tab);
	if(exp_nsld_curves.empty())
		hideTab(neutron_SL_tab);
	if(xray_curves.empty())
		hideTab(xray_formfactors_tab);
	if(neutron_curves.empty())
		hideTab(neutron_formfactors_tab);

	//reset booleans
	form_factors_done = false;
	sim_open = false;
	sim_from_trajectory = false;
}

//Slot that opens window to get trajectory file info
void STEwindow::openTrajectoryFile() {

	setEnabled(false); //don't want to be able to access main window while other is open

	//first need to get information, use window
	TrajectoryInfoWindow* info_window = new TrajectoryInfoWindow(this);
	info_window->show();

	//use different slots to deal with return from window
	connect(info_window, SIGNAL(canceled(TrajectoryInfoWindow*)), this, SLOT(cancelTrajectoryOpen(TrajectoryInfoWindow*)));
	connect(info_window, SIGNAL(dataEntered(TrajectoryInfoWindow*)), this, SLOT(readTrajectoryFile(TrajectoryInfoWindow*)));
}

//Slot that deals with a canceled trajectory file open
void STEwindow::cancelTrajectoryOpen(TrajectoryInfoWindow* window) {

	//need to delete the trajectory window and re-enable main window
	delete window;
	setEnabled(true);
}

//Slot that reads trajectory files after the using the info from window
void STEwindow::readTrajectoryFile(TrajectoryInfoWindow* window) {

	TrajectoryInfo info = window->trajectoryInfo();

	//done with window, delete it
	delete window;

	//get the info and pass the work along to appropriate function... use switch
	switch(info.trajectoryType()) {

	case GROMACS_XTC:

		readGromacsXTCFile(info);
		break;

	case NAMD_DCD:

		//nothing for now
		break;

	case AMBER_CRD:

		//nothing for now
		break;
	}

	//re-enable the main window
	setEnabled(true);
}

//Slot that opens Gromacs XTC file, takes in user info and calculates density info. Must also have GRO file!
void STEwindow::readGromacsXTCFile(TrajectoryInfo traj_info) {

	//if there is anything sim data already open, get rid of it
	if(sim_open)
		removeSim();

	//need an object to hold all of the data
	TrajectoryData gro_data(&info);

	//set the axis that is normal to the bilayer
	gro_data.setNormalVector(traj_info.normalAxis());

	//Set any bin sizes for the grids
	gro_data.setNormalGridSize(traj_info.normalGridSize());
	gro_data.setPolarGridSizes(traj_info.radialGridSize(), traj_info.angularGridSize());

	//set the fraction of heavy water to be used
	gro_data.setHeavyWaterFraction(traj_info.heavyWaterFraction());
	gro_data.setSimGranularity(traj_info.simGranularity());

	//first try to open and read gro file
	FileCode read_result = gro_data.readGroFile(traj_info.optFilePath());

	//only continue if everything is ok, otherwise delete and exit
	if(read_result == FILE_ERROR) {

		//then something wrong, let object deal with error message... delete it and return
		gro_data.clearFrames();
		return;
	}

	//now try reading xtc file, again if something wrong, delete and return
	//also create and pass along a progress bar that will be attached to the main window
	QProgressDialog progress(this);
	progress.setWindowModality(Qt::WindowModal);
	progress.setMinimumDuration(0); //always show the progress bars

	read_result = gro_data.readXtcFile(traj_info.filePath(), progress);
	if(read_result == FILE_ERROR) {

		//then something wrong, let object deal with error message... delete it and return
		gro_data.clearFrames();
		return;
	}

	//data has been read, want to compute the particle densities from the data
	if(traj_info.normalTimeAvg()) {

		//Do time average over all frames and give the info to the SimManager, it will calculate other densities, e.g. elec, neutron ,etc
		sim.setSimDataFromTrajectory(gro_data.ensembleDensities());

		//set lipid area and number from the trajectory data
		sim.setLipidArea(gro_data.averageLipidArea());
		sim.setLipidNumber(gro_data.numLipids());

		//SimManager now has all the info, emit signal that finished (takes care of display)
		sim_from_trajectory = true;
		emit simOpened();

	} else if(traj_info.normalFunctionOfTime()) {

		//do function of time
	}
}


//Slot that calculate the fourier transform of the current loaded data and then displays
//the results in the FF tabs
void STEwindow::displayFourier() {

	//always clear the previous form factor calculation
	cleanFourierCurves();

	//first check to see if any data has been loaded
	if(not sim_open) {

		QMessageBox::warning(this, "SIMtoEXP error:", "Must load simulation data prior to calculating structure factors");
		return;
	}

	//retrieve the range of q values for xray and neutron
	double x_qmin, x_qmax, x_qstep, n_qmin, n_qmax, n_qstep;

	QList<bool> oks; //for the six q line edits
	for(uint i = 0; i < NUMBER_OF_Q_LINE_EDITS; i++)
		oks.append(true); //need to initialize

	x_qmin = xray_qlow_lineEdit->text().toDouble(&oks[0]);
	x_qmax = xray_qhigh_lineEdit->text().toDouble(&oks[1]);
	x_qstep = xray_qstep_lineEdit->text().toDouble(&oks[2]);
	n_qmin = neut_qlow_lineEdit->text().toDouble(&oks[3]);
	n_qmax = neut_qhigh_lineEdit->text().toDouble(&oks[4]);
	n_qstep = neut_qstep_lineEdit->text().toDouble(&oks[5]);

	//check all of the input, if any invalid show box and return
	int size = oks.size();
	for(int i = 0; i < size; i++) {

		if(not oks[i]) {

			QMessageBox::warning(this, "Input error:", "One of the q-interval values is invalid.");
			return;
		}
	}

	//if any qmax < qmin, or qstep > qmax - qmin, show appropriate message box and return
	if(x_qmax <= x_qmin or n_qmax <= n_qmin) {

		QMessageBox::warning(this, "Input error:", "The 'High' value must be greater than the 'Low' value");
		return;
	}
	else if(x_qstep > (x_qmax - x_qmin) or n_qstep > (n_qmax - n_qmin)) {

		QMessageBox::warning(this, "Input error:", "Step size too big; must fit between 'High' and 'Low'");
		return;
	}

	sim.setXrayQvalues(x_qmin, x_qmax, x_qstep);
	sim.setNeutronQvalues(n_qmin, n_qmax, n_qstep);
	sim.setSymmetrized(sym_checkBox->isChecked());

	//solvent scattering densities
	bool x_ok, n_ok;
	double x_sd = solv_ed_lineEdit->text().toDouble(&x_ok);
	double n_sd = solv_nsld_lineEdit->text().toDouble(&n_ok);

	//check the input values for errors
	if(not(x_ok and n_ok)) {

		QMessageBox::warning(this, "Input error:", "One of the solvent values is invalid.");
		return;
	}

	sim.setSolventScatteringDensities(x_sd, n_sd);

	//everything has been passed to sim, so calculate fourier
	sim.calcFourierTransforms(this);

	//get name to pass to plot curve
	QString x_title ("Simulated");
	QwtPlotCurve* curve = new QwtPlotCurve(x_title);
	curve->setSamples(sim.getXrayFormFactors());

	//attach curve to the plot
	curve->attach(xray_formfactors_plot);

    QVariant itemInfo = xray_formfactors_plot->itemToInfo(curve);
    QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( xray_formfactor_legend.legendWidget( itemInfo ) );

	curve->setLegendAttribute(QwtPlotCurve::LegendShowLine);

	item->setChecked(true); //make all of them visible

	//put curves into containers to keep track of them (data pointers can be extracted from the curve)
	sim_FF_xray_curve = curve;

	//get name to pass to plot curve
	QString n_title ("Simulated");
	curve = new QwtPlotCurve(n_title);
	curve->setSamples(sim.getNeutronFormFactors());

	//attach curve to the plot
	curve->attach(neutron_formfactors_plot);

    itemInfo = neutron_formfactors_plot->itemToInfo(curve);
    item = qobject_cast<QwtLegendLabel *>( neutron_formfactor_legend.legendWidget( itemInfo ) );

	curve->setLegendAttribute(QwtPlotCurve::LegendShowLine);

	item->setChecked(true); //make all of them visible

	//put curves into containers to keep track of them (data pointers can be extracted from the curve)
	sim_FF_neutron_curve = curve;

	//need to tell any open experimental plots about the simulation plot in order to scale them
	exp.scaleData(sim.getXrayFormFactors(), sim.getNeutronFormFactors());
	form_factors_done = true;

	//update the experimental curves with the new scaled data (if they exist)
	updateXrayCurves();
	updateNeutronCurves();

	//update the plots
	updatePlot(tabIndex(xray_formfactors_tab));
	updatePlot(tabIndex(neutron_formfactors_tab));

	//change the current view to the xray tab, if NOT on the neutron already
	if(currentTabIndex() != tabIndex(neutron_formfactors_tab))
		setCurrentTab(xray_formfactors_tab);

	//make both tabs visible if not already
	showTab(xray_formfactors_tab);
	showTab(neutron_formfactors_tab);
}

//SLOT for opening xray data from file; if file opens properly, then the GUI is update by
//calling displayXray()
void STEwindow::openXray() {

	//first check to make sure there is room, if not, then display message and return
	if(xray_curves.size() >= MAX_CURVES) {

		QMessageBox::warning(this, "SIMtoEXP error:", "A maximum of 5 experimental data sets can be viewed at a time.\nTo analyze new data, an old data set must be deleted");
		return;
	}

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open X-ray structure factor file:"), "", "X-ray structure factor: *.xsf, *.xff (*.xff *.xsf);; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//returns null if 'cancel' is hit, so only continue if there is a path
	if(not fileName.isEmpty()) {

		//trigger simOpened only if read was successful
		if(exp.readExpFile(fileName, XRAY) == READ_SUCCESS) {

			uint index = exp.getNumberXray()-1;

			//also connect the ExpData to the QLineEdit object for updates
			ExpData& expData = exp.getXrayDataObject(index);
			connect(&expData, SIGNAL(scalingFactorChanged(QString)), xray_scale[index], SLOT(setText(QString)));
			connect(&expData, SIGNAL(chiSquaredChanged(QString)), xray_chi[index], SLOT(setText(QString)));

			//change the status of the line edits
			xray_scale[index]->setReadOnly(false);
			qutil::setBackgroundColour(xray_scale[index], QColor(Qt::white)); //change to white
			qutil::setBackgroundColour(xray_chi[index], QColor(Qt::white)); //change to white

			emit xrayOpened(index); //pass the index of the recently added data
		}
	}
}

//Displays the 
void STEwindow::displayXray(uint index) {

	//get name to pass to plot curve
	QString x_title ("Sample ");
	x_title += QString::number(index + 1);

	QwtPlotCurve* curve = new QwtPlotCurve(x_title);
	curve->setSamples(exp.getXrayData(index));

	//set to plot symbols only
	curve->setStyle(QwtPlotCurve::NoCurve);
	QwtSymbol* symbol = new QwtSymbol(getNextXraySymbol());

	//set colour
	QColor cur_colour = getNextXrayColour();
	symbol->setPen(QPen(cur_colour));
	symbol->setBrush(QBrush(cur_colour));
	symbol->setSize(SYMBOL_WIDTH);

	curve->setSymbol(symbol);
	xray_symbols.append(symbol);

	//attach curve to the plot
	curve->attach(xray_formfactors_plot);

	//add the current plot to the legend
	curve->setLegendAttribute(QwtPlotCurve::LegendShowSymbol);

    const QVariant itemInfo = xray_formfactors_plot->itemToInfo(curve);
    QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( xray_formfactor_legend.legendWidget( itemInfo ) );

    item->setChecked(true); //make all of them visible

	//keep track of pointer
	xray_curves.append(curve);

	//if simulation form factors present, need to scale the experimental results
	if(form_factors_done) {

		exp.scaleData(sim.getXrayFormFactors(), XRAY);
		updateXrayCurves();
	}

	//change the visibility of the QAction in the menu so it can be deleted
	xray_deletes[index]->setVisible(true);

	//update the plot info
	updatePlot(tabIndex(xray_formfactors_tab));
	showTab(xray_formfactors_tab);

	//finally, change the current view to the xray tab
	setCurrentTab(xray_formfactors_tab);
}

//SLOT for opening neutron data from file; if file opens properly, then the GUI is update by
//calling displayNeutron()
void STEwindow::openNeutron() {

	//first check to make sure there is room, if not, then display message and return
	if(neutron_curves.size() >= MAX_CURVES) {

		QMessageBox::warning(this, "SIMtoEXP error:", "A maximum of 5 experimental data sets can be viewed at a time.\nTo analyze new data, an old data set must be deleted");
		return;
	}

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open neutron structure factor file:"), "", "Neutron structure factor: *.nsf, *.nff (*.nsf *.nff);; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//returns null if 'cancel' is hit, so only continue if there is a path
	if(not fileName.isEmpty()) {

		//trigger simOpened only if read was successful
		if(exp.readExpFile(fileName, NEUTRON) == READ_SUCCESS) {

			unsigned int index = exp.getNumberNeutron()-1;

			//also connect the ExpData to the QLineEdit object for updates
			ExpData& expData = exp.getNeutronDataObject(index);
			connect(&expData, SIGNAL(scalingFactorChanged(QString)), neutron_scale[index], SLOT(setText(QString)));
			connect(&expData, SIGNAL(chiSquaredChanged(QString)), neutron_chi[index], SLOT(setText(QString)));

			//change the status of the line edits
			neutron_scale[index]->setReadOnly(false);
			setBackgroundColour(neutron_scale[index], QColor(Qt::white)); //change to white
			setBackgroundColour(neutron_chi[index], QColor(Qt::white)); //change to white

			emit neutronOpened(index); //pass the index of the recently added data
		}
	}
}

//Displays the 
void STEwindow::displayNeutron(uint index) {

	//get name to pass to plot curve
	QString x_title ("Sample ");
	x_title += QString::number(index + 1);

	QwtPlotCurve* curve = new QwtPlotCurve(x_title);
	curve->setSamples(exp.getNeutronData(index));

	//set to plot symbols only
	curve->setStyle(QwtPlotCurve::NoCurve);
	QwtSymbol* symbol = new QwtSymbol(getNextNeutronSymbol());

	//set colour
	QColor cur_colour = getNextNeutronColour();
	symbol->setPen(QPen(cur_colour));
	symbol->setBrush(QBrush(cur_colour));
	symbol->setSize(SYMBOL_WIDTH);

	curve->setSymbol(symbol);
	neutron_symbols.append(symbol);

	//attach curve to the plot
	curve->attach(neutron_formfactors_plot);

	//add the current plot to the legend
	curve->setLegendAttribute(QwtPlotCurve::LegendShowSymbol);

    const QVariant itemInfo = neutron_formfactors_plot->itemToInfo(curve);
    QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( neutron_formfactor_legend.legendWidget( itemInfo ) );

	item->setChecked(true); //make all of them visible

	//track the pointer for possible deletion later
	neutron_curves.append(curve);

	//if simulation form factors present, need to scale the experimental results
	if(form_factors_done) {

		exp.scaleData(sim.getNeutronFormFactors(), NEUTRON);
		updateNeutronCurves();
	}

	//change the visibility of the QAction in the menu so it can be deleted
	neutron_deletes[index]->setVisible(true);

	//update the plot info
	updatePlot(tabIndex(neutron_formfactors_tab));
	showTab(neutron_formfactors_tab);

	//finally, change the current view to the xray tab
	setCurrentTab(neutron_formfactors_tab);
}

//SLOT for opening and reading an experimental electron density file
void STEwindow::openElecDensity() {

	//first need to remove any data if it is there
	if(exp.electronDensityLoaded())
		deleteElecDensity();

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open electron density file:"), "", "Electron density profile: *.edp (*.edp);; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//returns null if 'cancel' is hit, so only continue if there is a path
	if(not fileName.isEmpty()) {

		//trigger elecDensityOpened only if read was successful
		if(exp.readExpDensityFile(fileName, XRAY) == READ_SUCCESS)
			emit elecDensityOpened();
	}
}

//SLOT for displaying the previously loaded data for the experimental electron density
void STEwindow::displayElecDensity() {

	for(unsigned int i = 0; i < exp.getNumberElectronDensity(); i++) {

		//create curve and provide the data
		QwtPlotCurve* curve = new QwtPlotCurve(exp.getElectronDensityData(i).first);
		curve->setSamples(exp.getElectronDensityData(i).second);

		curve->setPen(QPen (getColour(i*7))); //use a "random" colour, multiply by a prime

		//attach curve to the plot
		curve->attach(e_densities_plot);

		//add the current plot to the legend
        const QVariant itemInfo = e_densities_plot->itemToInfo(curve);
        QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( e_densities_legend.legendWidget( itemInfo ) );

        item->setChecked(true); //make all of them visible

		//put curves into containers to keep track of them (data pointers can be extracted from the curve)
		exp_ed_curves.append(curve);
	}

	//replot, also set bounds and zoombase
	updatePlot(tabIndex(e_densities_tab));
	showTab(e_densities_tab);

	//finally, change the current view to the electron densities tab
	setCurrentTab(e_densities_tab);
}

//SLOT for opening and reading an experimental neutron scattering length density file
void STEwindow::openNeutDensity() {

	//first need to remove any data if it there
	if(exp.neutronSLDensityLoaded())
		deleteNeutDensity();

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open neutron scattering length file:"), "", "Neutron scattering length profile: *.nsld (*.nsld);; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//returns null if 'cancel' is hit, so only continue if there is a path
	if(not fileName.isEmpty()) {

		//trigger elecDensityOpened only if read was successful
		if(exp.readExpDensityFile(fileName, NEUTRON) == READ_SUCCESS)
			emit neutDensityOpened();
	}
}

//SLOT for displaying the previously loaded data for the experimental neutron scattering length density
void STEwindow::displayNeutDensity() {

	for(unsigned int i = 0; i < exp.getNumberNeutronSLDensity(); i++) {

		//create curve and provide the data
		QwtPlotCurve* curve = new QwtPlotCurve(exp.getNeutronSLDensityData(i).first);
		curve->setSamples(exp.getNeutronSLDensityData(i).second);

		curve->setPen(QPen (getColour(i*7)));

		//attach curve to the plot
		curve->attach(neutron_SL_plot);

		//add the current plot to the legend
        const QVariant itemInfo = neutron_SL_plot->itemToInfo(curve);
        QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( neutron_SL_legend.legendWidget( itemInfo ) );

		item->setChecked(true); //make all of them visible

		//put curves into containers to keep track of them (data pointers can be extracted from the curve)
		exp_nsld_curves.append(curve);
	}

	//replot, also get bounds and set zoombase
	updatePlot(tabIndex(neutron_SL_tab));
	showTab(neutron_SL_tab);

	//finally, change the current view to the electron densities tab
	setCurrentTab(neutron_SL_tab);
}

//SLOT for opening and reading the component file
void STEwindow::openComponents() {

	//check to see if sim file even opened, otherwise do nothing
	if(not sim_open) {

		QMessageBox::warning(this, "SIMtoEXP error:", "You must open simulation data before loading components.");
		return;
	}

	//first need to remove any data if it there
	if(sim.componentFileLoaded())
		sim.removeComponents();

	bool file_ok = false;
	QString fileName;

	//do loop to either get a file that exists, or cancel the input
	while(not file_ok) {

		//use static method from file dialog to get simulation file name
		fileName = QFileDialog::getOpenFileName(this, tr("Open component file:"), "", "Component: *.cmp (*.cmp);; All files (*)");

		if(QFile::exists(fileName) or fileName.isNull())
			file_ok = true;
		else {
			QMessageBox::warning(this, "File selection error:", "The file selected for input does not exist.\nPlease select a different file or cancel input.");
		}
	}

	//returns null if 'cancel' is hit, so only continue if there is a path
	if(not fileName.isEmpty()) {

		//trigger elecDensityOpened only if read was successful
		if(sim.readComponentFile(fileName) == READ_SUCCESS)
			emit componentsOpened();
	}
}

//Displays the components on the three plots: 'Electron density', 'Neutron SL' and 'Number Density'
void STEwindow::displayComponents() {

	QList<Component> components = sim.getComponents();

	//first deal with number densities - want to hide all of the invididual atom curves and remove them from
	//the legend
	for(int i = 0; i < sim_number_curves.size(); i++) {

		sim_number_curves[i]->hide();
	}

	//then create new plots for the components
	for(int i = 0; i < components.size(); i++) {

		//for total only, make black, but don't make a curve for number densities
		QColor cur_colour;
		QwtPlotCurve* curve;
		if(components[i].name() == "Total")
			cur_colour = QColor (Qt::black);
		else
			cur_colour = getColour(i*47);

		QPen cur_pen(cur_colour);

		if(components[i].name() != "Total") {

			curve = new QwtPlotCurve(components[i].name());
			curve->setSamples(components[i].getNumberDensities());
			curve->attach(num_densities_plot);
			curve->setPen(cur_pen); //Make new QPen and give it a random colour, seed with a prime multiple

			//deal with legend
            const QVariant itemInfo = num_densities_plot->itemToInfo(curve);
            QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( num_densities_legend.legendWidget( itemInfo ) );

			item->setChecked(true); //make all of them visible

			//put curves into containers to keep track of them (data pointers can be extracted from the curve)
			comp_number_curves.append(curve);
		}

		//now do electron densities
		curve = new QwtPlotCurve(components[i].name());
		curve->setSamples(components[i].getElectronDensities());
		curve->setLegendAttribute(QwtPlotCurve::LegendShowLine); //must go prior to attachment, kinda stupid
		curve->attach(e_densities_plot);
		cur_pen.setStyle(Qt::DashLine);
		curve->setPen(cur_pen); //Make new QPen and give it a random colour, seed with a prime multiple

		//deal with legend
        QVariant itemInfo = e_densities_plot->itemToInfo(curve);
        QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( e_densities_legend.legendWidget( itemInfo ) );

		item->setChecked(true); //make all of them visible

		//put curves into containers to keep track of them (data pointers can be extracted from the curve)
		comp_ed_curves.append(curve);

		//and neutron SL
		curve = new QwtPlotCurve(components[i].name());
		curve->setSamples(components[i].getNeutronSLDensities());
		curve->setLegendAttribute(QwtPlotCurve::LegendShowLine); //must go prior to attachment, kinda stupid
		curve->attach(neutron_SL_plot);
		curve->setPen(cur_pen); //Make new QPen and give it a random colour, seed with a prime multiple

		//deal with legend
        itemInfo = neutron_SL_plot->itemToInfo(curve);
        item = qobject_cast<QwtLegendLabel *>( neutron_SL_legend.legendWidget( itemInfo ) );

		item->setChecked(true); //make all of them visible

		//put curves into containers to keep track of them (data pointers can be extracted from the curve)
		comp_nsld_curves.append(curve);
	}

	//update the three plots
	updatePlot(tabIndex(num_densities_tab));
	updatePlot(tabIndex(e_densities_tab));
	updatePlot(tabIndex(neutron_SL_tab));

	//focus on the number densities plot
	setCurrentTab(num_densities_tab);

	//also need to display the values in the lineEdits; just loop through and add
	//can only display a maximum of 10 components in the lineEdits
	uint max_index = id_edits.size();
	if(max_index > sim.numberComponents())
		max_index = sim.numberComponents();

	for(uint i = 0; i < max_index; i++) {

		id_edits[i]->setText(components[i].name());
		natoms_edits[i]->setText(QString::number(components[i].numberOfAtoms()));

		QString apg = QString::number(components[i].atomsPerGroup());
		apg.truncate(MAX_DIGITS);
		atoms_pgroup_edits[i]->setText(apg);

		ngroups_edits[i]->setText(QString::number(components[i].numberOfGroups()));
		nelec_pgroup_edits[i]->setText(QString::number(components[i].electronsPerGroup()));

		QString nsld;
		nsld.setNum(components[i].neutronSLPerGroup(),'E', DISPLAY_PRECISION);
		nsld_pgroup_edits[i]->setText(nsld);

		volume_edits[i]->setText(QString::number(components[i].volume()));

		setBackgroundColour(id_edits[i], QColor(Qt::white));
		setBackgroundColour(natoms_edits[i], QColor(Qt::white));
		setBackgroundColour(atoms_pgroup_edits[i], QColor(Qt::white));
		setBackgroundColour(ngroups_edits[i], QColor(Qt::white));
		setBackgroundColour(nelec_pgroup_edits[i], QColor(Qt::white));
		setBackgroundColour(nsld_pgroup_edits[i], QColor(Qt::white));
		setBackgroundColour(volume_edits[i], QColor(Qt::white));
	}
}

//SLOT for displaying the probability densities of the components
void STEwindow::displayVolumes() {

	QList<Component> components = sim.getComponents();

	//then create new plots for the components
	for(int i = 0; i < components.size(); i++) {

		//for total only, make black, but don't make a curve for number densities
		QColor cur_colour;
		QwtPlotCurve* curve;
		if(components[i].name() == "Total")
			cur_colour = QColor (Qt::black);
		else
			cur_colour = getColour(i*47);

		QPen cur_pen(cur_colour);

		curve = new QwtPlotCurve(components[i].name());

		if(components[i].name() != "Total")
			curve->setSamples(components[i].getProbabilities());
		else
			curve->setSamples(sim.getTotalComponentProbs());

		curve->attach(vol_probs_plot);
		curve->setPen(cur_pen); //Make new QPen and give it a random colour, seed with a prime multiple

		//deal with legend
        const QVariant itemInfo = vol_probs_plot->itemToInfo(curve);
        QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( volume_legend.legendWidget( itemInfo ) );

		item->setChecked(true); //make all of them checked

		//put curves into containers to keep track of them (data pointers can be extracted from the curve)
		volume_curves.append(curve);
	}

	//replot and set the bounds
	updatePlot(tabIndex(vol_probs_tab));
	showTab(vol_probs_tab);

	//finally, focus on the number densities plot
	setCurrentTab(vol_probs_tab);

	//also need to update the line edits
	uint max_index = id_edits.size();
	if(max_index > sim.numberComponents())
		max_index = sim.numberComponents();

	for(uint i = 0; i < max_index; i++)
		volume_edits[i]->setText(QString::number(components[i].volume()));

	//also update the RMS value line edit
	vol_rms_lineEdit->setText(QString::number(sim.volumeRMS()));
}

//private function that removes all of the component curves and component data from the program
void STEwindow::removeComponents() {

	//if no component file open, simply return
	if(not sim.componentFileLoaded())
		return;

	//go through all the curves, detach them and delete all their data
	for(int i = 0; i < comp_number_curves.size(); i++) {

		//detach them from the plot
		comp_number_curves[i]->detach();

		//remove them from the legend
	}

	for(int i = 0; i < comp_ed_curves.size(); i++) {

		//detach them from the plot
		comp_ed_curves[i]->detach();
	}

	for(int i = 0; i < comp_nsld_curves.size(); i++) {

		//detach them from the plot
		comp_nsld_curves[i]->detach();
	}

	//reset the lineEdits
	uint num_components = id_edits.size();
	for(uint i = 0; i < num_components; i++) {

		id_edits[i]->setText("");
		natoms_edits[i]->setText("0.0");
		atoms_pgroup_edits[i]->setText("0.0");
		ngroups_edits[i]->setText("0.0");
		nelec_pgroup_edits[i]->setText("0.0");
		nsld_pgroup_edits[i]->setText("0.0");
		volume_edits[i]->setText("0.0");

		setBackgroundColour(id_edits[i], QColor(Qt::lightGray));
		setBackgroundColour(natoms_edits[i], QColor(Qt::lightGray));
		setBackgroundColour(atoms_pgroup_edits[i], QColor(Qt::lightGray));
		setBackgroundColour(ngroups_edits[i], QColor(Qt::lightGray));
		setBackgroundColour(nelec_pgroup_edits[i], QColor(Qt::lightGray));
		setBackgroundColour(nsld_pgroup_edits[i], QColor(Qt::lightGray));
		setBackgroundColour(volume_edits[i], QColor(Qt::lightGray));
	}

	//reset the rms lineEdit
	vol_rms_lineEdit->setText("");

	qDeleteAll(comp_number_curves);
	qDeleteAll(comp_ed_curves);
	qDeleteAll(comp_nsld_curves);
	qDeleteAll(volume_curves);
	comp_number_curves.clear();
	comp_ed_curves.clear();
	comp_nsld_curves.clear();
	volume_curves.clear();

	//remove data from SimParser
	sim.removeComponents();

	//update the three plots
	updatePlot(tabIndex(num_densities_tab));
	updatePlot(tabIndex(e_densities_tab));
	updatePlot(tabIndex(neutron_SL_tab));

	//hide the tab
	hideTab(vol_probs_tab);

	//force a replot anyway so nothing weird happens
	vol_probs_plot->replot();
}

//SLOTs for save/printing to file the currently loaded data.  If nothing loaded, simply tell the user and return.
void STEwindow::writeElementalNumberDensities() {

	//first check for data
	if(sim_open)
		out_writer.writeNumberDensities(sim.getElementalData(), sim.numLipids(), sim.lipidArea());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}

void STEwindow::writeParticleNumberDensities() {

	//first check for data
	if(sim_open)
		out_writer.writeNumberDensities(sim.getParticleData(), sim.numLipids(), sim.lipidArea());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}

void STEwindow::writeComponentNumberDensities() {

	//first check for data
	if(sim.componentFileLoaded())
		out_writer.writeComponentNumberDensities(sim.getComponents());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}


void STEwindow::writeElementalElectronDensities() {

	//first check for data
	if(sim_open)
		out_writer.writeElectronDensities(sim.getElementalData());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}

void STEwindow::writeComponentElectronDensities() {

	//first check for data
	if(sim.componentFileLoaded())
		out_writer.writeComponentElectronDensities(sim.getComponents());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}

void STEwindow::writeElementalNeutronDensities() {

	//first check for data
	if(sim_open)
		out_writer.writeNeutronDensities(sim.getElementalData());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}

void STEwindow::writeComponentNeutronDensities() {

	//first check for data
	if(sim.componentFileLoaded())
		out_writer.writeComponentNeutronDensities(sim.getComponents());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}

void STEwindow::writeXrayFactors() {

	//first check for data
	if(form_factors_done)
		out_writer.writeXrayStructureFactors(sim.getXrayFormFactors());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");

}

void STEwindow::writeNeutronFactors() {

	//first check for data
	if(form_factors_done)
		out_writer.writeNeutronStructureFactors(sim.getNeutronFormFactors());

	//if not provide an error message
	else
		QMessageBox::warning(this, "Export error:", "There is no data loaded to save.");
}
void STEwindow::writeVolumetricFit() {}
void STEwindow::writeEverything() {}

//SLOT for deleting the current simulation data in the program
void STEwindow::deleteSim() {

	removeSim();
	removeComponents();

	if(tabHidden(currentTab()))
		goToFirstShownTab();

	//need to clear the margin line edits
	left_margin_lineEdit->setText("");
	right_margin_lineEdit->setText("");
}

//NEXT 10 functions: slots that pass along work to the functions that delete the curves
void STEwindow::deleteXraySample1() {deleteXrayData(0);}
void STEwindow::deleteXraySample2() {deleteXrayData(1);}
void STEwindow::deleteXraySample3() {deleteXrayData(2);}
void STEwindow::deleteXraySample4() {deleteXrayData(3);}
void STEwindow::deleteXraySample5() {deleteXrayData(4);}
void STEwindow::deleteNeutronSample1() {deleteNeutronData(0);}
void STEwindow::deleteNeutronSample2() {deleteNeutronData(1);}
void STEwindow::deleteNeutronSample3() {deleteNeutronData(2);}
void STEwindow::deleteNeutronSample4() {deleteNeutronData(3);}
void STEwindow::deleteNeutronSample5() {deleteNeutronData(4);}

//Deletes the curve provided by the index. No error checking, as there should be no way to pass
//the wrong index (i.e. out of bounds)
void STEwindow::deleteXrayData(uint index) {

	//index of largest curve
	int max_index = xray_curves.size() - 1;

	//first detach curve from plot
	xray_curves[index]->detach();

	//actually need to disconnect all ExpData objects from the lineEdits, and then reconnect below (somewhat annoying)
	uint num_xrays = exp.getNumberXray();
	for(unsigned int i = 0; i < num_xrays; i++) {

		ExpData& expData = exp.getXrayDataObject(i);
		disconnect(&expData, SIGNAL(scalingFactorChanged(QString)), xray_scale[i], SLOT(setText(QString)));
		disconnect(&expData, SIGNAL(chiSquaredChanged(QString)), xray_chi[i], SLOT(setText(QString)));
	}
	//delete the data from the ExpParser
	exp.removeXrayData(index);

	//delete the curve object
	delete xray_curves[index];

	//delete the curve pointer from the list
	xray_curves.removeAt(index);

	//rename the curves left in the list
	int xray_size = xray_curves.size();
	for(int i = 0; i < xray_size; i++) {

		//get name to pass to plot curve
		QString title ("Sample ");
		title += QString::number(i+1);
		xray_curves[i]->setTitle(title);
	}

	//hide the largest action button
	xray_deletes[max_index]->setVisible(false);

	//set the connections back up
	uint num_xray = exp.getNumberXray();
	for(uint i = 0; i < num_xray; i++) {

		ExpData& expData = exp.getXrayDataObject(i);
		connect(&expData, SIGNAL(scalingFactorChanged(QString)), xray_scale[i], SLOT(setText(QString)));
		connect(&expData, SIGNAL(chiSquaredChanged(QString)), xray_chi[i], SLOT(setText(QString)));
	}

	//shuffle all the lineEdit text past 'index' down one
	for(int i = index; i < max_index; i++) {

		xray_scale[i]->setText(xray_scale[i+1]->text());
		xray_chi[i]->setText(xray_chi[i+1]->text());
	}

	//remove the use of the last sample lineEdit
	xray_scale[max_index]->setReadOnly(true);
	setBackgroundColour(xray_scale[max_index], QColor (Qt::lightGray));
	xray_scale[max_index]->setText("1.0");
	xray_chi[max_index]->setReadOnly(true);
	setBackgroundColour(xray_chi[max_index], QColor (Qt::lightGray));
	xray_chi[max_index]->setText("0.0");

	//update the plot
	updatePlot(tabIndex(xray_formfactors_tab));

	//if nothing left, hide the tab and go to another one if on it
	if((sim_number_curves.empty() or not sim.formFactorsCalculated()) and xray_curves.empty()) {

		hideTab(xray_formfactors_tab);
		if(tabHidden(currentTab()))
			goToFirstShownTab();
	}
}

void STEwindow::deleteNeutronData(uint index) {

	//index of largest curve
	int max_index = neutron_curves.size() - 1;

	//first detach curve from plot
	neutron_curves[index]->detach();

	//actually need to disconnect all ExpData objects from the lineEdits, and then reconnect below (somewhat annoying)
	uint num_neutron = exp.getNumberNeutron();
	for(uint i = 0; i < num_neutron; i++) {

		ExpData& expData = exp.getNeutronDataObject(i);
		disconnect(&expData, SIGNAL(scalingFactorChanged(QString)), neutron_scale[i], SLOT(setText(QString)));
		disconnect(&expData, SIGNAL(chiSquaredChanged(QString)), neutron_chi[i], SLOT(setText(QString)));
	}
	//delete the data from the ExpParser
	exp.removeNeutronData(index);

	//delete the curve object
	delete neutron_curves[index];

	//delete the curve pointer from the list
	neutron_curves.removeAt(index);

	//rename the curves left in the list
	int neut_size = neutron_curves.size();
	for(int i = 0; i < neut_size; i++) {

		//get name to pass to plot curve
		QString title ("Sample ");
		title += QString::number(i+1);
		neutron_curves[i]->setTitle(title);
	}

	//hide the largest action button
	neutron_deletes[max_index]->setVisible(false);

	//set the connections back up
	num_neutron = exp.getNumberNeutron();
	for(uint i = 0; i < num_neutron; i++) {

		ExpData& expData = exp.getNeutronDataObject(i);
		connect(&expData, SIGNAL(scalingFactorChanged(QString)), neutron_scale[i], SLOT(setText(QString)));
		connect(&expData, SIGNAL(chiSquaredChanged(QString)), neutron_chi[i], SLOT(setText(QString)));
	}

	//shuffle all the lineEdit text past 'index' down one
	for(int i = index; i < max_index; i++) {

		neutron_scale[i]->setText(neutron_scale[i+1]->text());
		neutron_chi[i]->setText(neutron_chi[i+1]->text());
	}

	//remove the use of the last sample lineEdit
	neutron_scale[max_index]->setReadOnly(true);
	setBackgroundColour(neutron_scale[max_index], QColor (Qt::lightGray));
	neutron_scale[max_index]->setText("1.0");
	neutron_chi[max_index]->setReadOnly(true);
	setBackgroundColour(neutron_chi[max_index], QColor (Qt::lightGray));
	neutron_chi[max_index]->setText("0.0");

	//update the plot info
	updatePlot(tabIndex(neutron_formfactors_tab));

	//if nothing left, hide the tab and go to another one if on that tab
	if((sim_number_curves.empty() or not sim.formFactorsCalculated()) and neutron_curves.empty()) {

		hideTab(neutron_formfactors_tab);
		if(tabHidden(currentTab()))
			goToFirstShownTab();
	}
}

//Removes all of the experimental electron density data loaded
void STEwindow::deleteElecDensity() {

	//remove the curves from the plot and the legend
	int num_ed_curves = exp_ed_curves.size();
	for(int i = 0; i < num_ed_curves; i++) {

		exp_ed_curves[i]->detach();
	}

	//delete all of the curve objects and clear the container
	qDeleteAll(exp_ed_curves);
	exp_ed_curves.clear();

	//delete all of the stored data
	exp.removeElectronDensityData();

	//update the plot bounds and the plot, including zoom base
	updatePlot(tabIndex(e_densities_tab));

	//if nothing left, hide the tab and go to another one if on that tab
	if(sim_ed_curves.empty()) {

		hideTab(e_densities_tab);
		if(tabHidden(currentTab()))
			goToFirstShownTab();
	}
}

//Removes all of the experimental neutron scattering length density data loaded
void STEwindow::deleteNeutDensity() {

	//remove the curves from the plot and the legend
	int num_nsld_curves = exp_nsld_curves.size();
	for(int i = 0; i < num_nsld_curves; i++) {

		exp_nsld_curves[i]->detach();
	}

	//delete all of the curve objects and clear the container
	qDeleteAll(exp_nsld_curves);
	exp_nsld_curves.clear();

	//delete all of the stored data
	exp.removeNeutronSLDensityData();

	//update the plot bounds and the plot, including zoom base
	updatePlot(tabIndex(neutron_SL_tab));

	//if nothing left, hide the tab and go to another one if on that tab
	if(sim_nsld_curves.empty()) {

		hideTab(neutron_SL_tab);
		if(tabHidden(currentTab()))
			goToFirstShownTab();
	}
}

//SLOT that delete the components if loaded
void STEwindow::deleteComponents() {

	removeComponents();

	//also need to return the number density plots; make them visible and put back in the legend
	int num_curves = sim_number_curves.size();
	for(int i = 0; i < num_curves; i++) {

		sim_number_curves[i]->show();
        const QVariant itemInfo = num_densities_plot->itemToInfo(sim_number_curves[i]);
        QwtLegendLabel *item = qobject_cast<QwtLegendLabel *>( num_densities_legend.legendWidget( itemInfo ) );
\
        if(i % PLOT_FREQUENCY == 0) {

			item->setChecked(true);
		}
		else {
			sim_number_curves[i]->hide();
            item->setChecked(false);
        }
	}
}

//For debugging - really to help determine which button is named what!
void STEwindow::testButton() {

	qDebug() << "Connected button was just pressed!";
}

//Helper function that returns a random colour - only meant to be used
//in the number densities.  The AtomTypes should have set colours as per
//AtomicInfo.
QColor STEwindow::getColour(int seed) {

	double c1, c2, c3;
	srand(seed);

	c1 = (double)rand()/RAND_MAX;
	c2 = (double)rand()/RAND_MAX;
	c3 = (double)rand()/RAND_MAX;

	RGB colour(c1,c2,c3); //use to convert to integer

	return QColor(colour.getIntRed(),colour.getIntGreen(),colour.getIntBlue());
}

//Updates only the xray curves
void STEwindow::updateXrayCurves() {

	//Go through the curve pointer list, and update each curve with the new data
	for(int i = 0; i < xray_curves.size(); i++)
		xray_curves[i]->setSamples(exp.getXrayData(i));
}

//Updates only the neutron curves
void STEwindow::updateNeutronCurves() {

	int size = neutron_curves.size();
	for(int i = 0; i < size; i++)
		neutron_curves[i]->setSamples(exp.getNeutronData(i));
}

//Just forces a replot.
void STEwindow::replotXrayFormFactor() {

	xray_formfactors_plot->replot();
}

//Updates the info on the GUI for the scaling, and forces a replot.
void STEwindow::replotNeutronFormFactor() {

	neutron_formfactors_plot->replot();
}

//helper function for cleaning out all form factor curves.  If not done yet,
//don't do anything
void STEwindow::cleanFourierCurves() {

	if(form_factors_done) {

		//eliminate from the SimParser
		sim.clearFourier();

		//detach both curves
		sim_FF_xray_curve->detach();
		sim_FF_neutron_curve->detach();

		//delete them
		delete sim_FF_xray_curve;
		delete sim_FF_neutron_curve;

		form_factors_done = false;
	}
}

//Next 2 - Slots that scales each experimental form factor data set using the values provided in the GUI
void STEwindow::scaleXrayUser() {

	//first get all of the scaling factors
	Doubles scale_factors;

	uint num_xray = exp.getNumberXray();
	for(uint i = 0; i < num_xray; i++)
		scale_factors.push_back(xray_scale[i]->text().toDouble());

//	//if sim data exists, pass it along, otherwise don't bother
//	if(form_factors_done)
//		exp.scaleXrayData(scale_factors, &(sim.getXrayFormFactors()));
//	else
//		exp.scaleXrayData(scale_factors);

	exp.scaleXrayData(scale_factors);

	//update everything
	updateXrayCurves();
	updatePlot(tabIndex(xray_formfactors_tab));
}

void STEwindow::scaleNeutronUser() {

	//first get all of the scaling factors
	Doubles scale_factors;

	uint num_neut = exp.getNumberNeutron();
	for(uint i = 0; i < num_neut; i++)
		scale_factors.push_back(neutron_scale[i]->text().toDouble());

	//if sim data exists, pass it along, otherwise don't bother
//	if(form_factors_done)
//		exp.scaleNeutronData(scale_factors, &(sim.getNeutronFormFactors()));
//	else
//		exp.scaleNeutronData(scale_factors);

	exp.scaleNeutronData(scale_factors);

	updateNeutronCurves();
	updatePlot(tabIndex(neutron_formfactors_tab));
}

//return the next available colour in the colour container and increments for next call
QColor STEwindow::getNextXrayColour() {

	RGB rgb = *xray_ff_colour_itr;

	//increment for next time, cycle if necessary
	xray_ff_colour_itr++;
	if(xray_ff_colour_itr == plot_colours.end())
		xray_ff_colour_itr = plot_colours.begin();

	return QColor(rgb.getIntRed(), rgb.getIntGreen(), rgb.getIntBlue());
}

//return the next available colour in the colour container
QColor STEwindow::getNextNeutronColour() {

	RGB rgb = *neutron_ff_colour_itr;

	//increment for next time
	neutron_ff_colour_itr++;
	if(neutron_ff_colour_itr == plot_colours.end())
		neutron_ff_colour_itr = plot_colours.begin();

	return QColor(rgb.getIntRed(), rgb.getIntGreen(), rgb.getIntBlue());
}

//return the next available symbol in the container, and increment for next call
QwtSymbol::Style STEwindow::getNextXraySymbol() {

	QwtSymbol::Style style = *xray_ff_symbol_itr;

	xray_ff_symbol_itr++;
	if(xray_ff_symbol_itr == plot_symbols.end())
		xray_ff_symbol_itr = plot_symbols.begin();

	return style;
}

//return the next available symbol in the container, and increment for next call
QwtSymbol::Style STEwindow::getNextNeutronSymbol() {

	QwtSymbol::Style style = *neutron_ff_symbol_itr;

	neutron_ff_symbol_itr++;
	if(neutron_ff_symbol_itr == plot_symbols.end())
		neutron_ff_symbol_itr = plot_symbols.begin();

	return style;
}

//helper function that creates a vector of colours for use with plotting the data.  Only
//five colours available since only five exp samples can be plotted at a time
void STEwindow::buildColourLibrary() {

	plot_colours.push_back(RGB (1, 0, 0));
	plot_colours.push_back(RGB (0.25, 0.75, 0.25));
	plot_colours.push_back(RGB (0, 0, 1));
	plot_colours.push_back(RGB (1, 0, 1));
	plot_colours.push_back(RGB (0, 1, 1));
}

//helper function that creates a vector of symbol styles for use with plotting. Again, only five
//available due to the limit of five plots at once
void STEwindow::buildSymbolLibrary() {

	plot_symbols.push_back(QwtSymbol::Ellipse);
	plot_symbols.push_back(QwtSymbol::Rect);
	plot_symbols.push_back(QwtSymbol::Diamond);
	plot_symbols.push_back(QwtSymbol::UTriangle);
	plot_symbols.push_back(QwtSymbol::DTriangle);
}

void STEwindow::buildDeleteActionsContainers() {

	xray_deletes.push_back(actionDeleteXray1);
	xray_deletes.push_back(actionDeleteXray2);
	xray_deletes.push_back(actionDeleteXray3);
	xray_deletes.push_back(actionDeleteXray4);
	xray_deletes.push_back(actionDeleteXray5);

	neutron_deletes.push_back(actionDeleteNeutron1);
	neutron_deletes.push_back(actionDeleteNeutron2);
	neutron_deletes.push_back(actionDeleteNeutron3);
	neutron_deletes.push_back(actionDeleteNeutron4);
	neutron_deletes.push_back(actionDeleteNeutron5);
}

//helper function that puts specific gui QLineEdit objects into specific container, for making connections later
//NOTE: all of these QWidgets could have been created and stored in container, would be much nice, but
//this is how QtDesigner creates the objects, so don't want to mess with gui.h too much (out of laziness). 
void STEwindow::buildLineEditContainers() {

	xray_scale.push_back(xray_scale_sam1_lineEdit);
	xray_scale.push_back(xray_scale_sam2_lineEdit);
	xray_scale.push_back(xray_scale_sam3_lineEdit);
	xray_scale.push_back(xray_scale_sam4_lineEdit);
	xray_scale.push_back(xray_scale_sam5_lineEdit);

	xray_chi.push_back(xray_chi_sam1_lineEdit);
	xray_chi.push_back(xray_chi_sam2_lineEdit);
	xray_chi.push_back(xray_chi_sam3_lineEdit);
	xray_chi.push_back(xray_chi_sam4_lineEdit);
	xray_chi.push_back(xray_chi_sam5_lineEdit);

	neutron_scale.push_back(neut_scale_sam1_lineEdit);
	neutron_scale.push_back(neut_scale_sam2_lineEdit);
	neutron_scale.push_back(neut_scale_sam3_lineEdit);
	neutron_scale.push_back(neut_scale_sam4_lineEdit);
	neutron_scale.push_back(neut_scale_sam5_lineEdit);

	neutron_chi.push_back(neut_chi_sam1_lineEdit);
	neutron_chi.push_back(neut_chi_sam2_lineEdit);
	neutron_chi.push_back(neut_chi_sam3_lineEdit);
	neutron_chi.push_back(neut_chi_sam4_lineEdit);
	neutron_chi.push_back(neut_chi_sam5_lineEdit);

	//set all of these lineEdit to grey background and uneditable until samples are opened
	int i;
	int size = xray_scale.size();
	for(i = 0; i < size; i++) {

		xray_scale[i]->setReadOnly(true);
		setBackgroundColour(xray_scale[i], QColor (Qt::lightGray));
	}

	size = xray_chi.size();
	for(i = 0; i < size; i++) {

		xray_chi[i]->setReadOnly(true);
		setBackgroundColour(xray_chi[i], QColor (Qt::lightGray));
	}

	size = neutron_scale.size();
	for(i = 0; i < size; i++) {

		neutron_scale[i]->setReadOnly(true);
		setBackgroundColour(neutron_scale[i], QColor (Qt::lightGray));
	}

	size = xray_chi.size();
	for(i = 0; i < size; i++) {

		neutron_chi[i]->setReadOnly(true);
		setBackgroundColour(neutron_chi[i], QColor (Qt::lightGray));
	}

	//set component background colours to grey until they have something to display
	size = id_edits.size();
	for(i = 0; i < size; i++)
		setBackgroundColour(id_edits[i], QColor (Qt::lightGray));

	size = natoms_edits.size();
	for(i = 0; i < size; i++)
		setBackgroundColour(natoms_edits[i], QColor (Qt::lightGray));

	size = atoms_pgroup_edits.size();
	for(i = 0; i < size; i++)
		setBackgroundColour(atoms_pgroup_edits[i], QColor (Qt::lightGray));

	size = nelec_pgroup_edits.size();
	for(i = 0; i < size; i++)
		setBackgroundColour(nelec_pgroup_edits[i], QColor (Qt::lightGray));

	size = nsld_pgroup_edits.size();
	for(i = 0; i < size; i++)
		setBackgroundColour(nsld_pgroup_edits[i], QColor (Qt::lightGray));

	size = ngroups_edits.size();
	for(int i = 0; i < size; i++)
		setBackgroundColour(ngroups_edits[i], QColor (Qt::lightGray));

	size = volume_edits.size();
	for(int i = 0; i < size; i++)
		setBackgroundColour(volume_edits[i], QColor (Qt::lightGray));
}

//sets the margins using the data currently available for each type of data
QRectF STEwindow::setNumDensityBoundingRect() {

	//just get from SimManager
	bounding_rects[tabIndex(num_densities_tab)] = sim.numDensitiesRect();

	//set all plot boundaries to the rectangle with margins added
	QRectF rect_with_margins = addMargins(bounding_rects[tabIndex(num_densities_tab)]);

	num_densities_plot->setAxisScale(QwtPlot::xBottom, rect_with_margins.left(), rect_with_margins.right());
	num_densities_plot->setAxisScale(QwtPlot::yLeft, rect_with_margins.bottom(), rect_with_margins.top());

	return rect_with_margins;
}

//Assumes Simulation bounds are up to date!
QRectF STEwindow::setElectronDensityBoundingRect() {

	const QRectF& sim_rect = sim.electronDensityRect();
	const QRectF& exp_rect = exp.electronDensityRect();

	//if no experimental data loaded, don't compare
	if(not exp.electronDensityLoaded()) {

		bounding_rects[tabIndex(e_densities_tab)] = sim_rect;

	} else if(not sim.dataLoaded()) {

		bounding_rects[tabIndex(e_densities_tab)] = exp_rect;

	} else {

		double left, right, top, bottom;

		//determine bounds
		if(sim_rect.left() < exp_rect.left())
			left = sim_rect.left();
		else
			left = exp_rect.left();

		if(sim_rect.right() > exp_rect.right())
			right = sim_rect.right();
		else
			right = exp_rect.right();

		if(sim_rect.bottom() < exp_rect.bottom())
			bottom = sim_rect.bottom();
		else
			bottom = exp_rect.bottom();

		if(sim_rect.top() > exp_rect.top())
			top = sim_rect.top();
		else
			top = exp_rect.top();

		//set bounds on rect
		bounding_rects[tabIndex(e_densities_tab)].setCoords(left, top, right, bottom);
	}

	//set all plot boundaries to the rectangle with margins added
	QRectF rect_with_margins = addMargins(bounding_rects[tabIndex(e_densities_tab)]);

	e_densities_plot->setAxisScale(QwtPlot::xBottom, rect_with_margins.left(), rect_with_margins.right());
	e_densities_plot->setAxisScale(QwtPlot::yLeft, rect_with_margins.bottom(), rect_with_margins.top());

	return rect_with_margins;
}

//Assumes Simulation margins are up to date!
QRectF STEwindow::setNeutronSLDensityBoundingRect() {

	const QRectF& sim_rect = sim.neutronSLDensityRect();
	const QRectF& exp_rect = exp.neutronSLDensityRect();

	if(not exp.neutronSLDensityLoaded()) {

		bounding_rects[tabIndex(neutron_SL_tab)] = sim_rect;

	} else if(not sim.dataLoaded()) {

		bounding_rects[tabIndex(e_densities_tab)] = exp_rect;
	}

	else {

		double left, right, top, bottom;

		//determine the bounds
		if(sim_rect.left() < exp_rect.left())
			left = sim_rect.left();
		else
			left = exp_rect.left();

		if(sim_rect.right() > exp_rect.right())
			right = sim_rect.right();
		else
			right = exp_rect.right();

		if(sim_rect.bottom() < exp_rect.bottom())
			bottom = sim_rect.bottom();
		else
			bottom = exp_rect.bottom();

		if(sim_rect.top() > exp_rect.top())
			top = sim_rect.top();
		else
			top = exp_rect.top();

		//set the rect
		bounding_rects[tabIndex(neutron_SL_tab)].setCoords(left, top, right, bottom);
	}

	//set all plot boundaries to the rectangle with margins added
	QRectF rect_with_margins = addMargins(bounding_rects[tabIndex(neutron_SL_tab)]);

	neutron_SL_plot->setAxisScale(QwtPlot::xBottom, rect_with_margins.left(), rect_with_margins.right());
	neutron_SL_plot->setAxisScale(QwtPlot::yLeft, rect_with_margins.bottom(), rect_with_margins.top());

	return rect_with_margins;
}

//set the xray bounds using both experimental and simulation data
QRectF STEwindow::setXrayFormFactorBoundingRect() {

	//to keep from calculating over and over again below
	const QRectF& sim_rect = sim.xrayFormFactorRect();
	const QRectF& exp_rect = exp.xrayFormFactorRect();

	//just make a default unit box - this should probably never happen anyway
	if(exp.xrayFormFactorEmpty() and not sim.formFactorsCalculated()) {

		bounding_rects[tabIndex(xray_formfactors_tab)].setCoords(0, 1, 1, 0);

	} else if(exp.xrayFormFactorEmpty()) {

		//then just set to the calculated simulation values
		bounding_rects[tabIndex(xray_formfactors_tab)] = sim_rect;

	} else if(not sim.formFactorsCalculated()) {

		//then just set to the experimental values
		bounding_rects[tabIndex(xray_formfactors_tab)] = exp_rect;

	} else {

		double left, right, top, bottom;

		//determine the bounds
		if(sim_rect.left() < exp_rect.left())
			left = sim_rect.left();
		else
			left = exp_rect.left();

		if(sim_rect.right() > exp_rect.right())
			right = sim_rect.right();
		else
			right = exp_rect.right();

		if(sim_rect.bottom() < exp_rect.bottom())
			bottom = sim_rect.bottom();
		else
			bottom = exp_rect.bottom();

		if(sim_rect.top() > exp_rect.top())
			top = sim_rect.top();
		else
			top = exp_rect.top();

		//set the rect
		bounding_rects[tabIndex(xray_formfactors_tab)].setCoords(left, top, right, bottom);
	}

	//set all plot boundaries to the rectangle with margins added
	QRectF rect_with_margins = addMargins(bounding_rects[tabIndex(xray_formfactors_tab)]);

	xray_formfactors_plot->setAxisScale(QwtPlot::xBottom, rect_with_margins.left(), rect_with_margins.right());
	xray_formfactors_plot->setAxisScale(QwtPlot::yLeft, rect_with_margins.bottom(), rect_with_margins.top());

	return rect_with_margins;
}

//set the neutron bounds using both experimental and simulation data
QRectF STEwindow::setNeutronFormFactorBoundingRect() {

	//to keep from calculating over and over again below
	const QRectF& sim_rect = sim.neutronFormFactorRect();
	const QRectF& exp_rect = exp.neutronFormFactorRect();

	if(exp.neutronFormFactorEmpty() and not sim.formFactorsCalculated()) {

		bounding_rects[tabIndex(neutron_formfactors_tab)].setCoords(0, 1, 1, 0);

	} else if(exp.neutronFormFactorEmpty()) {

		//then just set to the calculated simulation values
		bounding_rects[tabIndex(neutron_formfactors_tab)] = sim_rect;

	} else if(not sim.formFactorsCalculated()) {

		//then just set to the experimental values
		bounding_rects[tabIndex(neutron_formfactors_tab)] = exp_rect;

	} else {

		double left, right, top, bottom;

		//determine the bounds
		if(sim_rect.left() < exp_rect.left())
			left = sim_rect.left();
		else
			left = exp_rect.left();

		if(sim_rect.right() > exp_rect.right())
			right = sim_rect.right();
		else
			right = exp_rect.right();

		if(sim_rect.bottom() < exp_rect.bottom())
			bottom = sim_rect.bottom();
		else
			bottom = exp_rect.bottom();

		if(sim_rect.top() > exp_rect.top())
			top = sim_rect.top();
		else
			top = exp_rect.top();

		//set the rect
		bounding_rects[tabIndex(neutron_formfactors_tab)].setCoords(left, top, right, bottom);
	}

	//set all plot boundaries to the rectangle with margins added
	QRectF rect_with_margins = addMargins(bounding_rects[tabIndex(neutron_formfactors_tab)]);

	neutron_formfactors_plot->setAxisScale(QwtPlot::xBottom, rect_with_margins.left(), rect_with_margins.right());
	neutron_formfactors_plot->setAxisScale(QwtPlot::yLeft, rect_with_margins.bottom(), rect_with_margins.top());

	return rect_with_margins;
}

//Assumes Simulation bounds are up to date!
QRectF STEwindow::setVolumeBoundingRect() {

	//go through all the curves, get the bounding rectangles and calculate the bounding box
	double left = DBL_MAX, right = DBL_MIN, top = DBL_MIN, bottom = DBL_MAX;
	for(int i = 0; i < volume_curves.size(); i++) {

		//Rect is calculated, so set to a temp variable
		QRectF cur_rect = volume_curves[i]->boundingRect();

		//Qwt has a bug, for the curves, top = bottom and bottom = top
		if(cur_rect.left() < left)
			left = cur_rect.left();

		if(cur_rect.right() > right)
			right = cur_rect.right();

		if(cur_rect.bottom() > top)
			top = cur_rect.bottom();

		if(cur_rect.top() < bottom)
			bottom = cur_rect.top();
	}

	bounding_rects[tabIndex(vol_probs_tab)] = QRectF (QPointF(left,top), QPointF(right,bottom));

	//set plot axis to the bounding box plus margins
	QRectF rect_with_margins = addMargins(bounding_rects[tabIndex(vol_probs_tab)]);

	vol_probs_plot->setAxisScale(QwtPlot::xBottom, rect_with_margins.left(), rect_with_margins.right());
	vol_probs_plot->setAxisScale(QwtPlot::yLeft, rect_with_margins.bottom(), rect_with_margins.top());

	return rect_with_margins;
}

//SLOT for calculating component volumes and displaying the result in the vol_probs_plot
void STEwindow::calcComponentVolumes() {

	//first check to see if components loaded, if not provide warning
	if(not sim.componentFileLoaded()) {
		QMessageBox::warning(0, "SIMtoEXP error:", "Components must be loaded to calculate volumes.");
		return;
	}

	//if they have been loaded, but already calculated, then do nothing
	if(not sim.volumeCalculated()) {

		sim.solveVolumeEquations();
		displayVolumes();
	}
}

//private function that puts the plots and legends into a map container using the tab index as the key
void STEwindow::buildPlotContainers() {

	plots.insert(tabIndex(e_densities_tab), e_densities_plot);
	plots.insert(tabIndex(neutron_SL_tab), neutron_SL_plot);
	plots.insert(tabIndex(xray_formfactors_tab), xray_formfactors_plot);
	plots.insert(tabIndex(neutron_formfactors_tab), neutron_formfactors_plot);
	plots.insert(tabIndex(num_densities_tab), num_densities_plot);
	plots.insert(tabIndex(vol_probs_tab), vol_probs_plot);
}

//private function for populating the list of ScrollZoomer objects that are attached to each plot
void STEwindow::buildScrollZoomerContainer() {

	scroll_zoomers.insert(tabIndex(e_densities_tab), e_densities_scroll);
	scroll_zoomers.insert(tabIndex(neutron_SL_tab), neutron_SL_scroll);
	scroll_zoomers.insert(tabIndex(xray_formfactors_tab), xray_formfactors_scroll);
	scroll_zoomers.insert(tabIndex(neutron_formfactors_tab), neutron_formfactors_scroll);
	scroll_zoomers.insert(tabIndex(num_densities_tab), num_densities_scroll);
	scroll_zoomers.insert(tabIndex(vol_probs_tab), vol_probs_scroll);
}

//private function that takes care of initializing all of the rectangles used
//to set the bounding / zooming boxes
void STEwindow::initializeBoundingRects() {
      
	//initialize with a rectangle with zero size/area
	bounding_rects.insert(tabIndex(e_densities_tab), QRectF(QPointF(0,0), QPointF(0,0)));
	bounding_rects.insert(tabIndex(neutron_SL_tab), QRectF(QPointF(0,0), QPointF(0,0)));
	bounding_rects.insert(tabIndex(xray_formfactors_tab), QRectF(QPointF(0,0), QPointF(0,0)));
	bounding_rects.insert(tabIndex(neutron_formfactors_tab), QRectF(QPointF(0,0), QPointF(0,0)));
	bounding_rects.insert(tabIndex(num_densities_tab), QRectF(QPointF(0,0), QPointF(0,0)));
	bounding_rects.insert(tabIndex(vol_probs_tab), QRectF(QPointF(0,0), QPointF(0,0)));
}

//SLOT: sets the zoom back to the original bounding box with margins
void STEwindow::setZoomToBounds() {

	int cur_tab = currentTabIndex();
	scroll_zoomers[cur_tab]->zoom(scroll_zoomers[cur_tab]->zoomBase());
}

//Private function that returns a new rectangle equivalent to the passed rectangle with the margins added
QRectF STEwindow::addMargins(const QRectF& rect) {

	QRectF with_margins;

	//QRectF has width and height functions, but for some reason seem to be buggy (probably due to the severe challenge
	//of calculating them).  Anyway, need to calculate here
	double width = rect.right() - rect.left();
	double height = rect.top() - rect.bottom();

	with_margins.setLeft(rect.left() - X_AXIS_MARGIN_FRACTION * width * 0.5);
	with_margins.setRight(rect.right() + X_AXIS_MARGIN_FRACTION * width * 0.5);
	with_margins.setBottom(rect.bottom() - (Y_AXIS_MARGIN_FRACTION * height * 0.5));
	with_margins.setTop(rect.top() + (Y_AXIS_MARGIN_FRACTION * height * 0.5));

	return with_margins;
}

//SLOT: resets the Fourier section line edits to default values defined above
void STEwindow::setFourierDefaults() {

    xray_qlow_lineEdit->setText(QString::number(XRAY_QMIN, 'f', 1));
    xray_qhigh_lineEdit->setText(QString::number(XRAY_QMAX, 'f', 1));
    xray_qstep_lineEdit->setText(QString::number(XRAY_QSTEP));
    neut_qlow_lineEdit->setText(QString::number(NEUT_QMIN, 'f', 1));
    neut_qhigh_lineEdit->setText(QString::number(NEUT_QMAX));
    neut_qstep_lineEdit->setText(QString::number(NEUT_QSTEP));
    solv_ed_lineEdit->setText(QString::number(WATER_ED));
    solv_nsld_lineEdit->setText(QString::number(WATER_NSLD, 'E', 2));
}

//private functions that passes the work of determining the bounding rectangle along to the appropriate function
QRectF STEwindow::setBoundingRect(int index) {

	//have to use if/else instead of switch, since the tab indexes are technically not constant
	if(index == tabIndex(num_densities_tab))
		return setNumDensityBoundingRect();

	else if(index == tabIndex(e_densities_tab))
		return setElectronDensityBoundingRect();

	else if(index == tabIndex(neutron_SL_tab))
		return setNeutronSLDensityBoundingRect();

	else if(index == tabIndex(xray_formfactors_tab))
		return setXrayFormFactorBoundingRect();

	else if(index == tabIndex(neutron_formfactors_tab))
		return setNeutronFormFactorBoundingRect();

	else if(index == tabIndex(vol_probs_tab))
		return setVolumeBoundingRect();

	//a default rect
	else
		return QRectF();
}

//private function that updates the plot boundaries (including margins), replots the curves and sets the zoom
//to the correct base value.
void STEwindow::updatePlot(int index) {

	QRectF rect = setBoundingRect(index);
	plots[index]->replot();
	scroll_zoomers[index]->setZoomBase(rect);
}

//SLOT - for toggling whether the tabs should take up the whole program window or just the upper portion
void STEwindow::togglePlotArea() {

	//accomplish this by hiding the bottom layout if it is shown, and vice-versa - very simple :)
	if(not tab_area_maximized) {

	    margins->hide();
	    simFormFactors->hide();
	    expFormFactors->hide();
	    components->hide();

	    //toggle boolean
	    tab_area_maximized = true;

	} else {

	    margins->show();
	    simFormFactors->show();
	    expFormFactors->show();
	    components->show();

	    //toggle boolean
		tab_area_maximized = false;
	}
}

//SLOT - for displaying a new tab with the scattering density values that can be viewed/edited
void STEwindow::displayScatteringDensities() {

	//check to see if a tab is already open, if so do nothing
	if(scattering_tab_open)
		return;

	//create a new tab and add it to the tab widget, also make it the current tab
    Tab* scat_tab = new Tab(Scattering);
    tabWidget->addTab(scat_tab, QApplication::translate("MainWindow", "Scattering Densities", 0 ));//QApplication::UnicodeUTF8));
    scat_tab->setTabIndex(tabWidget->indexOf(scat_tab));
    tabWidget->setCurrentIndex(tabWidget->indexOf(scat_tab));

    //add tab to list
    tabs.append(scat_tab);

    //create a connection from the "close" button to the slot that does the work
    connect(tabWidget->closeButton(tabWidget->indexOf(scat_tab)), SIGNAL(clicked(bool)), scat_tab, SLOT(closeTab()));
    connect(scat_tab, SIGNAL(needToClose(int)), this, SLOT(closeTab(int)));

    //Use a separate class for dealing with tab contents and their i/o - this object must be a child of the
    //tab to ensure it is deleted properly when closed (along with the rest of the dynamically
    //created objects below)
    ScatteringTabFrame* scat_frame = new ScatteringTabFrame(&info, scat_tab);

    //put the frame in a scroll area for when the window is reduced in size
    QScrollArea* scat_scroll = new QScrollArea(scat_tab);
    scat_scroll->setWidget(scat_frame);

    //put in a grid layout to fill the tab properly (although it doesn't expand to fill it, not sure why)
    QGridLayout* scat_grid = new QGridLayout(scat_tab);
    scat_grid->addWidget(scat_scroll);

    scat_scroll->show();

    //set open bool
    scattering_tab_open = true;
}

//SLOT - closes the tab for which the 'close' button was pressed and deletes the page widget
void STEwindow::closeTab(int index) {

	//get tab type and change bool
	TabName cur_name = ((Tab*) tabWidget->widget(index))->getTabName();
	resetTabBoolean(cur_name);

	//remove the appropriate tab and delete the object
	tabWidget->removeTab(index);

	int list_index = index - NUMBER_OF_DEFAULT_PLOTS;
	tabs.removeAt(list_index);

	//if any dynamically created tabs remaining, need to renumber them
	if(not tabs.empty()) {

		int size = tabs.size();
		for(int i = 0; i < size; i++) {

			int index = tabIndex(tabs[i]); //get it
			tabs[i]->setTabIndex(index); //set it
		}
	}
}

//private function that provides a switch for turning off the appropriate QTab boolean
void STEwindow::resetTabBoolean(TabName name) {

	switch(name) {

	case Scattering:

		scattering_tab_open = false;
		break;

	case AtomicFF:

		atomicff_tab_open = false;
		break;

	case SimSpecifics:

		sim_specifics_tab_open = false;
		break;

	case Unknown:
		break;
	}
}

//private functions for displaying or making tabs unavailable
void STEwindow::showTab(QWidget* tab) {

	//if not already enabled, do so
	if(not tabWidget->isTabEnabled(tabIndex(tab)))
        tabWidget->setTabEnabled(tabIndex(tab), true);
}

void STEwindow::hideTab(QWidget* tab) {

	//if already enabled, turn off
	if(tabWidget->isTabEnabled(tabIndex(tab)))
        tabWidget->setTabEnabled(tabIndex(tab), false);
}

//private functions for whether the tab is enabled / disabled
bool STEwindow::tabShown(QWidget* tab) {

	return tabWidget->isTabEnabled(tabIndex(tab));
}

bool STEwindow::tabHidden(QWidget* tab) {

	return not tabShown(tab);
}

//simple function that returns the index of the tab in the tab widget
int STEwindow::tabIndex(QWidget* tab) {

	return tabWidget->indexOf(tab);
}

//just return the current tab
int STEwindow::currentTabIndex() {

	return tabWidget->currentIndex();
}

QWidget* STEwindow::currentTab() {

	return tabWidget->currentWidget();
}

//simple function that sets the current tab
void STEwindow::setCurrentTab(QWidget* tab) {

	tabWidget->setCurrentIndex(tabIndex(tab));
}

//function that finds the first enabled tab and makes it the current tab.  If nothing enabled
//just sets to index = 0
void STEwindow::goToFirstShownTab() {

	//unfortunately, QTabWidget does not give a list of tabs, so need to go through them manually
	if(tabShown(num_densities_tab))
		setCurrentTab(num_densities_tab);

	else if(tabShown(e_densities_tab))
		setCurrentTab(e_densities_tab);

	else if(tabShown(neutron_SL_tab))
		setCurrentTab(neutron_SL_tab);

	else if(tabShown(xray_formfactors_tab))
		setCurrentTab(xray_formfactors_tab);

	else if(tabShown(neutron_formfactors_tab))
		setCurrentTab(neutron_formfactors_tab);

	else if(tabShown(vol_probs_tab))
		setCurrentTab(vol_probs_tab);

	//if there are extra tabs, just go to the first one, otherwise set to zero
	else if(not tabs.empty())
		setCurrentTab(tabs[0]);

	else
		setCurrentTab(num_densities_tab);

}

//private slot that sets the boolean 'sim_open' to true
void STEwindow::setSimOpenTrue() {sim_open = true;}

#if QT_VERSION >= 0x050000
void STEwindow::showEDensities(const QVariant &itemInfo, bool on, int)
{
    QwtPlotItem *plotItem = e_densities_plot->infoToItem( itemInfo );
    if ( plotItem )
        plotItem->setVisible( on );
        plotItem->plot()->replot();
}

void STEwindow::showNeutronSLDensities(const QVariant &itemInfo, bool on, int)
{
    QwtPlotItem *plotItem = neutron_SL_plot->infoToItem( itemInfo );
    if ( plotItem )
        plotItem->setVisible( on );
        plotItem->plot()->replot();
}

void STEwindow::showXRayFF(const QVariant &itemInfo, bool on, int)
{
    QwtPlotItem *plotItem = xray_formfactors_plot->infoToItem( itemInfo );
    if ( plotItem )
        plotItem->setVisible( on );
        plotItem->plot()->replot();
}

void STEwindow::showNeutronFF(const QVariant &itemInfo, bool on, int)
{
    QwtPlotItem *plotItem = neutron_formfactors_plot->infoToItem( itemInfo );
    if ( plotItem )
        plotItem->setVisible( on );
        plotItem->plot()->replot();
}

void STEwindow::showNumDensities(const QVariant &itemInfo, bool on, int)
{
    QwtPlotItem *plotItem = num_densities_plot->infoToItem( itemInfo );
    if ( plotItem )
        plotItem->setVisible( on );
        plotItem->plot()->replot();
}

void STEwindow::showVolProbs(const QVariant &itemInfo, bool on, int)
{
    QwtPlotItem *plotItem = vol_probs_plot->infoToItem( itemInfo );
    if ( plotItem )
        plotItem->setVisible( on );
        plotItem->plot()->replot();
}
#else
//SLOT for setting visibility of plot items
void STEwindow::setVisibility(QwtPlotItem* item, bool value) {

    item->setVisible(value);
    item->plot()->replot(); //replot the plot associated with the item
}
#endif

//private function for making all static connections that are created with the main window
void STEwindow::makeConnections() {

	//'open' options -> functions that read in data
	connect(actionSIMtoEXP_format, SIGNAL(triggered(bool)), this, SLOT(openSim()));
	connect(actionTrajFiles, SIGNAL(triggered(bool)), this, SLOT(openTrajectoryFile()));

	//function that reads in data ->  function that displays data
	connect(this, SIGNAL(simOpened()), this, SLOT(displaySim()));
	connect(this, SIGNAL(simOpened()), this, SLOT(setSimOpenTrue()));

	//'Fourier transform' button -> perform transform and display
	connect(fourier_button, SIGNAL(clicked(bool)), this, SLOT(displayFourier()));
	connect(defaults_button, SIGNAL(clicked(bool)), this, SLOT(setFourierDefaults()));

	//open experimental action -> function that reads in data
	connect(actionX_ray, SIGNAL(triggered(bool)), this, SLOT(openXray()));
	connect(this, SIGNAL(xrayOpened(unsigned int)), this, SLOT(displayXray(unsigned int)));

	connect(actionNeutron, SIGNAL(triggered(bool)), this, SLOT(openNeutron()));
	connect(this, SIGNAL(neutronOpened(unsigned int)), this, SLOT(displayNeutron(unsigned int)));

	//open experimental scattering densities -> function that reads in data
	connect(actionElectron_density_profile, SIGNAL(triggered(bool)), this, SLOT(openElecDensity()));
	connect(this, SIGNAL(elecDensityOpened()), this, SLOT(displayElecDensity()));
	connect(actionNeutron_Scattering_Length_Density, SIGNAL(triggered(bool)), this, SLOT(openNeutDensity()));
	connect(this, SIGNAL(neutDensityOpened()), this, SLOT(displayNeutDensity()));

	//open component file
	connect(actionComponent_file, SIGNAL(triggered(bool)), this, SLOT(openComponents()));
	connect(this, SIGNAL(componentsOpened()), this, SLOT(displayComponents()));

	//for updating the scaling info after scaling is finished
	connect(&exp, SIGNAL(xrayDataScaled()), this, SLOT(replotXrayFormFactor()));
	connect(&exp, SIGNAL(neutronDataScaled()), this, SLOT(replotNeutronFormFactor()));

	//scales the form factor plots manually as per the GUI input
	connect(xray_scale_button, SIGNAL(clicked(bool)), this, SLOT(scaleXrayUser()));
	connect(neut_scale_button, SIGNAL(clicked(bool)), this, SLOT(scaleNeutronUser()));

	//connect volume button to calculation
	connect(calc_volume_button, SIGNAL(clicked(bool)), this, SLOT(calcComponentVolumes()));

    //for checking/unchecking curves on all the plots
#if QT_VERSION >= 0x050000
    connect(&e_densities_legend, SIGNAL(checked( const QVariant &, bool, int )), this, SLOT(showEDensities(const QVariant &,bool,int)));
    connect(&neutron_SL_legend, SIGNAL(checked( const QVariant &, bool, int )), this, SLOT(showNeutronSLDensities(const QVariant &,bool,int)));
    connect(&xray_formfactor_legend, SIGNAL(checked( const QVariant &, bool, int )), this, SLOT(showXRayFF(const QVariant &,bool,int)));
    connect(&neutron_formfactor_legend, SIGNAL(checked( const QVariant &, bool, int )), this, SLOT(showNeutronFF(const QVariant &,bool,int)));
    connect(&num_densities_legend, SIGNAL(checked( const QVariant &, bool, int )), this, SLOT(showNumDensities(const QVariant &,bool,int)));
    connect(&volume_legend, SIGNAL(checked( const QVariant &, bool, int )), this, SLOT(showVolProbs(const QVariant &,bool,int)));
#else
    connect(e_densities_plot, SIGNAL(legendChecked(QwtPlotItem*,bool)), this, SLOT(setVisibility(QwtPlotItem*,bool)));
    connect(neutron_SL_plot, SIGNAL(legendChecked(QwtPlotItem*,bool)), this, SLOT(setVisibility(QwtPlotItem*,bool)));
    connect(xray_formfactors_plot, SIGNAL(legendChecked(QwtPlotItem*,bool)), this, SLOT(setVisibility(QwtPlotItem*,bool)));
    connect(neutron_formfactors_plot, SIGNAL(legendChecked(QwtPlotItem*,bool)), this, SLOT(setVisibility(QwtPlotItem*,bool)));
    connect(num_densities_plot, SIGNAL(legendChecked(QwtPlotItem*,bool)), this, SLOT(setVisibility(QwtPlotItem*,bool)));
    connect(vol_probs_plot, SIGNAL(legendChecked(QwtPlotItem*,bool)), this, SLOT(setVisibility(QwtPlotItem*,bool)));
#endif

	//for exporting/saving data
	connect(actionElementalNumb_density, SIGNAL(triggered(bool)), this, SLOT(writeElementalNumberDensities()));
	connect(actionParticleNumb_density, SIGNAL(triggered(bool)), this, SLOT(writeParticleNumberDensities()));
	connect(actionComponentNumb_density, SIGNAL(triggered(bool)), this, SLOT(writeComponentNumberDensities()));
	connect(actionElementalElectron_density, SIGNAL(triggered(bool)), this, SLOT(writeElementalElectronDensities()));
	connect(actionComponentElectron_density, SIGNAL(triggered(bool)), this, SLOT(writeComponentElectronDensities()));
	connect(actionElementalNeutron_SL_density, SIGNAL(triggered(bool)), this, SLOT(writeElementalNeutronDensities()));
	connect(actionComponentNeutron_SL_density, SIGNAL(triggered(bool)), this, SLOT(writeComponentNeutronDensities()));
	connect(action_X_ray_form_factors, SIGNAL(triggered(bool)), this, SLOT(writeXrayFactors()));
	connect(actionNeutron_form_factors, SIGNAL(triggered(bool)), this, SLOT(writeNeutronFactors()));
	connect(actionVolumetric_fit, SIGNAL(triggered(bool)), this, SLOT(writeVolumetricFit()));
	connect(actionEverything, SIGNAL(triggered(bool)), this, SLOT(writeEverything()));

	//bunch of connections for deleting curves/data
	connect(actionDeleteSim, SIGNAL(triggered(bool)), this, SLOT(deleteSim()));
	connect(actionDeleteED, SIGNAL(triggered(bool)), this, SLOT(deleteElecDensity()));
	connect(actionDeleteNSLD, SIGNAL(triggered(bool)), this, SLOT(deleteNeutDensity()));
	connect(actionDeleteXray1, SIGNAL(triggered(bool)), this, SLOT(deleteXraySample1()));
	connect(actionDeleteXray2, SIGNAL(triggered(bool)), this, SLOT(deleteXraySample2()));
	connect(actionDeleteXray3, SIGNAL(triggered(bool)), this, SLOT(deleteXraySample3()));
	connect(actionDeleteXray4, SIGNAL(triggered(bool)), this, SLOT(deleteXraySample4()));
	connect(actionDeleteXray5, SIGNAL(triggered(bool)), this, SLOT(deleteXraySample5()));
	connect(actionDeleteNeutron1, SIGNAL(triggered(bool)), this, SLOT(deleteNeutronSample1()));
	connect(actionDeleteNeutron2, SIGNAL(triggered(bool)), this, SLOT(deleteNeutronSample2()));
	connect(actionDeleteNeutron3, SIGNAL(triggered(bool)), this, SLOT(deleteNeutronSample3()));
	connect(actionDeleteNeutron4, SIGNAL(triggered(bool)), this, SLOT(deleteNeutronSample4()));
	connect(actionDeleteNeutron5, SIGNAL(triggered(bool)), this, SLOT(deleteNeutronSample5()));
	connect(actionDeleteComponents, SIGNAL(triggered(bool)), this, SLOT(deleteComponents()));

	//for plot actions
	connect(actionFitToWindow, SIGNAL(triggered(bool)), this, SLOT(setZoomToBounds()));
	connect(actionTogglePlot, SIGNAL(triggered(bool)), this, SLOT(togglePlotArea()));

	//for displaying tools
	connect(actionScattering_lengths, SIGNAL(triggered(bool)), this, SLOT(displayScatteringDensities()));
}
