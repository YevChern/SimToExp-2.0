/*
 * SimParser.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: bholland
 */ 

#include "SimManager.h"
#include "LinearSolver.h"
#include "Util.h"

#include <cmath>
#include <iostream>
#include <limits>

#include <QtCore/QIODevice>
#include <QtCore/QtDebug>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QMessageBox>
#endif

using namespace std;
using namespace qutil;

//Default Constructor
SimManager::SimManager() {

	sim_data_loaded = false;
	sim_data_from_traj = false;
	atomTypesSet = false;
	compFileLoaded = false;
	form_factors_calculated = false;
	volCalculated = false;
	traj_num_lipids = 0;
	traj_lipid_area = 0;

	//set the pointer to NULL for now just to initialize
	xray_fourier_FF = NULL;
	neut_fourier_FF = NULL;
	info = NULL;

	//initialize the bounds
	num_densities_rect = QRectF (QPointF(0,0), QPointF(0,0));
	ed_rect = QRectF (QPointF(0,0), QPointF(0,0));
	nsld_rect = QRectF (QPointF(0,0), QPointF(0,0));
	xrayFF_rect = QRectF (QPointF(0,0), QPointF(0,0));
	neutFF_rect = QRectF (QPointF(0,0), QPointF(0,0));

	//initialize volume rms
	volume_rms = 0;
}

//Destructor
SimManager::~SimManager() {

	//created dynamically in the constructor
	if(xray_fourier_FF != 0)
		delete xray_fourier_FF;
	if(neut_fourier_FF != 0)
		delete neut_fourier_FF;
}

//setter for the info, must be done right away for the FormFactor objects
//not the best way to do this
void SimManager::setAtomInfo(AtomicInfo& info) {

	this->info = &info;

	xray_fourier_FF = new FourierFormFactor(XRAY, info);
	neut_fourier_FF = new FourierFormFactor(NEUTRON, info);
}

//setter for whether the fourier transform should be symmetrized or not.  Sets both to the same value.
void SimManager::setSymmetrized(bool symmetrized) {

	xray_fourier_FF->setSymmetrized(symmetrized);
	neut_fourier_FF->setSymmetrized(symmetrized);
}

//setters for the 'q' values
void SimManager::setXrayQvalues(double max, double min, double step) {

	xray_fourier_FF->setQvalues(max,min,step);
}

void SimManager::setNeutronQvalues(double max, double min, double step) {

	neut_fourier_FF->setQvalues(max,min,step);
}

//set the solvent scattering densities for both FourierFormFactor objects
void SimManager::setSolventScatteringDensities(double electron_d, double nsld) {

	xray_fourier_FF->setSolventScatteringDensity(electron_d);
	neut_fourier_FF->setSolventScatteringDensity(nsld);
}

//Return the elemental data read in.  If nothing read yet, returns an empty vector
NormalDensities SimManager::getElementalData() const {

	//return the appropriate list of NormalDensity objects depending on how data was procured
	if(sim_data_from_traj)
		return ensemble_frames[elemental_index].normalDensities();

	return reader.getParticleDensities();
}

//Return the particle data read in.  If nothing read yet, returns an empty vector
NormalDensities SimManager::getParticleData() const {

	//return the appropriate list of NormalDensity objects depending on how data was procured
	if(sim_data_from_traj)
		return ensemble_frames[particle_index].normalDensities();

	return reader.getParticleDensities();
}

//just returns whether or not the simulation data is currently loaded
bool SimManager::dataLoaded() const {return sim_data_loaded;}
bool SimManager::simDataIsFromTrajectory() const {return sim_data_from_traj;}

//setters for trajectory data - should be set after trajectory is read in
void SimManager::setLipidArea(double area) {

	traj_lipid_area = area;
}

void SimManager::setLipidNumber(uint number) {

	traj_num_lipids = number;
}

//getters for same values
double SimManager::lipidArea() const {

	if(sim_data_from_traj)
		return traj_lipid_area;

	return reader.lipidArea();
}

uint SimManager::numLipids() const {

	if(sim_data_from_traj)
		return traj_num_lipids;

	return reader.numLipids();
}


//return the components read in by 'readComponentFile'.  If nothing read yet, returns an empty QList
Components& SimManager::getComponents() {

	return reader.getComponents();
}

const Components& SimManager::getComponents() const {

	return reader.getComponents();
}

//Boolean for whether the form factors have been calculated for the current data
bool SimManager::formFactorsCalculated() const {return form_factors_calculated;}

//Returns the number of components currently loaded - "Total" is NOT included
uint SimManager::numberComponents() const {

	if(not compFileLoaded)
		return 0;

	uint num_components = reader.getComponents().size() - 1;
	return num_components;
}

//Return a vector list of all of the AtomTypes seen in the sim file (and thus in the maps).
//Returned in the order they appear in the maps
const AtomTypes& SimManager::getAtomTypes() {

	//if container not yet filled, do it
	if(not atomTypesSet) {

		AtomTypes types = total_elec_densities.keys();

		int size = types.size();
		for(int i = 0; i < size; ++i) {
			atom_types.append(types[i]);
		}

		atomTypesSet = true;
	}

	return atom_types;
}

//Returns the number of atom types read in; if 'getAtomTypes' not called previously, this
//should always return zero
uint SimManager::numberAtomTypes() const {

	return atom_types.size();
}

//Next two return total density functions for the given atom type
Points& SimManager::getTotalElectronDensity(AtomType type) {

	return *total_elec_densities.find(type);
}

const Points& SimManager::getTotalElectronDensity(AtomType type) const {

	return *total_elec_densities.find(type);
}

Points& SimManager::getTotalNeutronDensity(AtomType type) {

	return *total_neut_densities.find(type);
}

const Points& SimManager::getTotalNeutronDensity(AtomType type) const {

	return *total_neut_densities.find(type);
}

//Return the calculated form factor functions
Points& SimManager::getXrayFormFactors() {

	return xray_fourier_FF->getFormFactors();
}

const Points& SimManager::getXrayFormFactors() const {

	return xray_fourier_FF->getFormFactors();
}

Points& SimManager::getNeutronFormFactors() {

	return neut_fourier_FF->getFormFactors();
}

const Points& SimManager::getNeutronFormFactors() const {

	return neut_fourier_FF->getFormFactors();
}

bool SimManager::volumeCalculated() const {return volCalculated;}
double SimManager::volumeRMS() const {return volume_rms;}

//Opens and reads the data in the file given by 'path'.
FileCode SimManager::readSimFile(const QString& path) {

	//pass job onto reader
	FileCode code = reader.readSimFile(path);

	//if everything okay, do calculations
	if(code == READ_SUCCESS) {

		//now that everything has been read in properly, go through sim data and calculate electron
		//and neutron SL densities, and also set up atomic info
		NormalDensity* density = reader.getFirstParticleDensity();
		for( ; not reader.atEnd(); density = reader.getNextParticleDensity()) {

			setNormalDensityAtomInfo(*density); //set up the atomic info
			density->calcAllDensities(false);
		}

		//finally, go through and calculate the TOTAL density functions for each AtomType and store them
		calcTotalDensities();
		sim_data_loaded = true;
	}

	return code;
}

/*
 * Sets the information read from the trajectory file; this is for the ensemble calculation
 *
 * MUST follow the order:
 *
 * 1) elemental frame
 * 2) particle frame
 */
void SimManager::setSimDataFromTrajectory(const DensityFrames& frame) {

	//first just set the given frame to stored frame
	ensemble_frames = frame;

	DensityFrame& elem_frame = ensemble_frames[elemental_index];
	DensityFrame& particle_frame = ensemble_frames[particle_index];

	//adjust the frames so that every bin is represented by all particles
	alignDensityFrame(elem_frame);
	alignDensityFrame(particle_frame);

	//finally, go through and calculate the TOTAL density functions for each AtomType and store them
	sim_data_from_traj = true;
	calcTrajDataTotalDensities();
	sim_data_loaded = true;
}

/*
 * Private function that determines the z-axis limits of the given DensityFrame object, and trims
 * off any NormalDensity object that goes beyond the limits. Basically needed to ensure there are
 * no bins that some particles were never sampled in.
 */
void SimManager::alignDensityFrame(DensityFrame& frame) {

	//first go through all of the calculated densities and ensure they have the
	//same number of bins and bin values, trim as necessary
	float min_bin = -numeric_limits<float>::max();
	float max_bin = numeric_limits<float>::max();

	//loop through once to set the min and max
	NormalDensity* density = frame.firstNormalDensity();
	for(; not frame.normalAtEnd(); density = frame.nextNormalDensity()) {

		float cur_min = density->firstPoint().x();
		float cur_max = density->lastPoint().x();

		//want the BIGGEST minimum
		if(cur_min > min_bin)
			min_bin = cur_min;

		//want the SMALLEST maximum
		if(cur_max < max_bin)
			max_bin = cur_max;
	}

	//go through sim data and calculate electron and neutron SL densities, while trimming if necessary
	for(density = frame.firstNormalDensity(); not frame.normalAtEnd(); density = frame.nextNormalDensity()) {

		//trim densities if necessary
		density->trimIndependentAxis(min_bin, max_bin);
		density->calcAllDensities(true);
	}
}



FileCode SimManager::readComponentFile(const QString& path) {

	//pass job onto reader
	FileCode code = reader.readComponentFile(path);

	//if everything okay, set boolean
	if(code == READ_SUCCESS)
		compFileLoaded = true;

	return code;
}

//Clears the manager of all data read from previous file
void SimManager::removeData() {

	//should go through and clear all memory of each
	reader.clearSimData();
	total_elec_densities.clear();
	total_neut_densities.clear();
	clearFourier();

	//if data has been loaded from trajectory, delete necessary traj objects
	ensemble_frames.clear(); //TODO: change to clearAll() when polar introduced
	time_function_frames.clear();

	//reset atomTypes
	atom_types.clear();
	atomTypesSet = false;
	sim_data_loaded = false;
	sim_data_from_traj = false;
}

//only clears the FourierFormFactor objects; for calculating from the same data
void SimManager::clearFourier() {

	xray_fourier_FF->clear();
	neut_fourier_FF->clear();

	form_factors_calculated = false;
}

//Calculates the total density functions for each AtomType and stores them in maps for easy access later
void SimManager::calcTotalDensities() {

	//go through the container of NormalDensity objects
	NormalDensity* density = reader.getFirstParticleDensity();

    //create empty vector for total_elec_densities
    //we have no total density yet, so set all density values to 0
    Points dummy_electrons = density->electronDensities();
    for (int i=0; i<dummy_electrons.size(); ++i) {
        dummy_electrons[i].setY(0.0);
    }

    //create empty vector for total_neut_densities
    //we have no total density yet, so set all density values to 0
    Points dummy_neutrons = density->neutronDensities();
    for (int i=0; i<dummy_neutrons.size(); ++i) {
        dummy_neutrons[i].setY(0.0);
    }

	//add totals to the density map container
    total_elec_densities.insert(SYSTEM_TOTAL, dummy_electrons);
    total_neut_densities.insert(SYSTEM_TOTAL, dummy_neutrons);

	for(; not reader.atEnd(); density = reader.getNextParticleDensity()) {

		insertDensities(density);
	}
}

//Calculates the total density functions for each AtomType and stores them in maps for easy access later
//This is different from above in that it uses the data from the trajectory read in, not the 'sim' file
void SimManager::calcTrajDataTotalDensities() {

	//go through the container of NormalDensity objects
	DensityFrame& elem_frame = ensemble_frames[elemental_index];
	NormalDensity* density = elem_frame.firstNormalDensity();

    //create empty vector for total_elec_densities
    //we have no total density yet, so set all density values to 0
    Points dummy_electrons = density->electronDensities();
    for (int i=0; i<dummy_electrons.size(); ++i) {
        dummy_electrons[i].setY(0.0);
    }

    //create empty vector for total_neut_densities
    //we have no total density yet, so set all density values to 0
    Points dummy_neutrons = density->neutronDensities();
    for (int i=0; i<dummy_neutrons.size(); ++i) {
        dummy_neutrons[i].setY(0.0);
    }

	//add totals to the density map container
    total_elec_densities.insert(SYSTEM_TOTAL, dummy_electrons);
    total_neut_densities.insert(SYSTEM_TOTAL, dummy_neutrons);

	for(; not elem_frame.normalAtEnd(); density = elem_frame.nextNormalDensity()) {

		insertDensities(density);
	}
}

/*
 * Helper function that does the work of actually inserting the densities once they are obtained
 * from the appropriate source
 */
void SimManager::insertDensities(NormalDensity* density) {

	TotalDensities::iterator e_itr, n_itr, e_total_itr, n_total_itr;
	AtomType type = density->atomType();
	e_itr = total_elec_densities.find(type);
	n_itr = total_neut_densities.find(type);

	const Points& cur_elec = density->electronDensities();
	const Points& cur_neut = density->neutronDensities();

	int cur_size = cur_elec.size();

	//if the AtomType not yet created, then do it here
	//if not in one container, then it shouldn't be in the other
	if(e_itr == total_elec_densities.end()) {

		total_elec_densities.insert(type, cur_elec);
		total_neut_densities.insert(type, cur_neut);

	} else {

		//otherwise go through and add to it

		//all vectors should be the same size, write out an error if not true
		//TODO: add error check

		for(int j = 0; j < cur_size; ++j) {

			double cur_ed = e_itr->at(j).y() + cur_elec[j].y();
			(*e_itr)[j].setY(cur_ed);

			double cur_nsld = n_itr->at(j).y() + cur_neut[j].y();
			(*n_itr)[j].setY(cur_nsld);
		}
	}

	//either way, add to the totals as well
	e_total_itr = total_elec_densities.find(SYSTEM_TOTAL);
	n_total_itr = total_neut_densities.find(SYSTEM_TOTAL);

	//check sizes to make sure there is no funny business
	if(cur_size != e_total_itr->size()) {

		cerr << "SimManager::insertDensities - Atom density and total density containers must be the same size!" << endl;
		cerr << "SimManager::insertDensities - Atom density size: " << cur_size << endl;
		cerr << "SimManager::insertDensities - Total density size: " << e_total_itr->size() << endl;
	}

	//otherwise go ahead and add
	for(int j = 0; j < cur_size; ++j) {

		double cur_ed = e_total_itr->at(j).y() + cur_elec[j].y();
		(*e_total_itr)[j].setY(cur_ed);

		double cur_nsld = n_total_itr->at(j).y() + cur_neut[j].y();
		(*n_total_itr)[j].setY(cur_nsld);
	}
}

//helper function for setting the AtomicInfo of the NormalDensity
void SimManager::setNormalDensityAtomInfo(NormalDensity& pd) {

	//TODO: see if the first part of this is necessary anymore with the new '/' notation, don't think it is
	//want to look for this as a marker for the AtomType in the SIM file
	char dash = '-';

	//set the three pieces of info for the atom
	string std_name = pd.particleName().toStdString();
	string particle_name = std_name.substr(0, std_name.find_first_of(dash));
	AtomType element = info->atomType(particle_name);

	//if not found, output an error box message and return
	if(element == UNKNOWN) {

		//first check to see if it is a MARTINI particle - the first two characters are used to describe the acyl chains
		string res_name = pd.residueName().toStdString();
		string mart_name = res_name.substr(0, 2) + "_" + std_name;
		element = info->atomType(mart_name);

		//if still unknown, create an error message box
		if(element == UNKNOWN) {

			QString msg = "The particle '";
			msg += pd.particleName();
			msg += "' is not recognized by SIMtoEXP, please see manual for all accepted particles.";
			QMessageBox::critical(0, "SIM File Input Error:", msg);
			return;
		}

		//if it is a MARTINI particle, change the name in the ParticleLabel
		pd.setParticleName(QString(mart_name.c_str()));
	}

	double num_electrons = info->numberElectrons(element);
	double neutron_sl = info->neutronSL(element);

	//pass to the ParticleDensity object
	pd.setAtomInfo(element, num_electrons, neutron_sl);
	return; //for debugger
}

//Assumes all data has been provided and calculates the transforms for both Xray and Neutron
void SimManager::calcFourierTransforms(QWidget* progress_parent) {

	//if the electron map is empty, return an error message for debug
	if(total_elec_densities.empty()) {

		cerr << "SimManager::calcFourierTransforms: Electric densities is empty, nothing to calculate." << endl;
		return;
	}

	//need to fill in the last bit of info first
	xray_fourier_FF->createAtomicFormFactors(getAtomTypes());
	neut_fourier_FF->createAtomicFormFactors(getAtomTypes());

	Points points = *total_elec_densities.begin();
	xray_fourier_FF->setRealSpaceValues(points);
	neut_fourier_FF->setRealSpaceValues(points);

	//create a progress dialog box to show the progress of the calculation
	QProgressDialog progress(progress_parent);
	progress.setWindowModality(Qt::WindowModal);
	progress.setLabelText("Calculating form factors...");

	//get the number of Qvalues for each situation
	uint num_xray_qvalues = xray_fourier_FF->getNumberQValues();
	uint num_total_qvalues = num_xray_qvalues + neut_fourier_FF->getNumberQValues();

	progress.setMaximum(num_total_qvalues);
	progress.setMinimumDuration(0); //always show the progress bars

	//everything should be set, now just calculate
	xray_fourier_FF->calcFormFactors(getElementalData(), progress, 0);
	neut_fourier_FF->calcFormFactors(getElementalData(), progress, num_xray_qvalues);

	//finalize the progress bar
	progress.setValue(num_total_qvalues);

	//set booleans
	form_factors_calculated = true;
}

//Return the bounding rectangle for the number density data by searching through the ParticleDensity objects, unless components
//have been loaded, then only goes through the components
const QRectF& SimManager::numDensitiesRect() {

	//if no components, then go through particles
	if(not compFileLoaded) {

		//need to go through all curves to get a proper bounding rectangle
		NormalDensities pds;
		if(sim_data_from_traj)
			pds = ensemble_frames[elemental_index].normalDensities();
		else
			pds = reader.getParticleDensities();

		//loop through, get the rectangle for each curve and compare to the current bounding box
		int size = pds.size();
		for(int i = 0; i < size; ++i) {

			QRectF cur_rect = getBoundingRect(pds[i].densities());

			//need to initialize the rect as (0,0) may not even be in final rectangle
			if(i == 0)
				num_densities_rect = cur_rect;

			else {

				if(num_densities_rect.left() > cur_rect.left())
					num_densities_rect.setLeft(cur_rect.left());

				if(num_densities_rect.right() < cur_rect.right())
					num_densities_rect.setRight(cur_rect.right());

				if(num_densities_rect.bottom() > cur_rect.bottom())
					num_densities_rect.setBottom(cur_rect.bottom());

				if(num_densities_rect.top() < cur_rect.top())
					num_densities_rect.setTop(cur_rect.top());
			}
		}
	}

	//otherwise go through components
	else {

		//need to go through all curves to get a proper bounding rectangle
		Components comps = reader.getComponents();

		//loop through, get the rectangle for each curve and compare to the current bounding box
		int comps_size = comps.size();
		for(int i = 0; i < comps_size; ++i) {

			QRectF cur_rect = getBoundingRect(comps[i].getNumberDensities());

			//need to initialize the rect as (0,0) may not even be in final rectangle
			if(i == 0)
				num_densities_rect = cur_rect;

			else {

				if(num_densities_rect.left() > cur_rect.left())
					num_densities_rect.setLeft(cur_rect.left());

				if(num_densities_rect.right() < cur_rect.right())
					num_densities_rect.setRight(cur_rect.right());

				if(num_densities_rect.bottom() > cur_rect.bottom())
					num_densities_rect.setBottom(cur_rect.bottom());

				if(num_densities_rect.top() < cur_rect.top())
					num_densities_rect.setTop(cur_rect.top());
			}
		}
	}

	return num_densities_rect;
}

//Return the bounding rectangle for the electron density data by searching through the curves
const QRectF& SimManager::electronDensityRect() {

	//need to go through all curves to get a proper bounding rectangle
	QMap<AtomType,Points>::iterator itr, itr_end = total_elec_densities.end();

	//loop through, get the rectangle for each curve and compare to the current bounding box
	for(itr = total_elec_densities.begin(); itr != itr_end; ++itr) {

		QRectF cur_rect = getBoundingRect(itr.value());

		//need to initialize the rect as (0,0) may not even be in final rectangle
		if(itr == total_elec_densities.begin())
			ed_rect = cur_rect;

		else {

			if(ed_rect.left() > cur_rect.left())
				ed_rect.setLeft(cur_rect.left());

			if(ed_rect.right() < cur_rect.right())
				ed_rect.setRight(cur_rect.right());

			if(ed_rect.bottom() > cur_rect.bottom())
				ed_rect.setBottom(cur_rect.bottom());

			if(ed_rect.top() < cur_rect.top())
				ed_rect.setTop(cur_rect.top());
		}
	}

	//if components are loaded, also need to go through them
	if(compFileLoaded) {

		Components comps = reader.getComponents();

		//loop through, get the rectangle for each curve and compare to the current bounding box
		int comps_size = comps.size();
		for(int i = 0; i < comps_size; ++i) {

			QRectF cur_rect = getBoundingRect(comps[i].getElectronDensities());

			if(ed_rect.left() > cur_rect.left())
				ed_rect.setLeft(cur_rect.left());

			if(ed_rect.right() < cur_rect.right())
				ed_rect.setRight(cur_rect.right());

			if(ed_rect.bottom() > cur_rect.bottom())
				ed_rect.setBottom(cur_rect.bottom());

			if(ed_rect.top() < cur_rect.top())
				ed_rect.setTop(cur_rect.top());
		}
	}

	return ed_rect;
}

//Return the bounding rectangle for the neutron scattering length density data by searching through the curves
const QRectF& SimManager::neutronSLDensityRect() {

	//need to go through all curves to get a proper bounding rectangle
	QMap<AtomType,Points>::iterator itr, itr_end = total_neut_densities.end();

	//loop through, get the rectangle for each curve and compare to the current bounding box
	for(itr = total_neut_densities.begin(); itr != itr_end; ++itr) {

		QRectF cur_rect = getBoundingRect(itr.value());

		//need to initialize the rect as (0,0) may not even be in final rectangle
		if(itr == total_neut_densities.begin())
			nsld_rect = cur_rect;

		else {

			if(nsld_rect.left() > cur_rect.left())
				nsld_rect.setLeft(cur_rect.left());

			if(nsld_rect.right() < cur_rect.right())
				nsld_rect.setRight(cur_rect.right());

			if(nsld_rect.bottom() > cur_rect.bottom())
				nsld_rect.setBottom(cur_rect.bottom());

			if(nsld_rect.top() < cur_rect.top())
				nsld_rect.setTop(cur_rect.top());
		}
	}

	//if components are loaded, also need to go through them
	if(compFileLoaded) {

		Components comps = reader.getComponents();

		//loop through, get the rectangle for each curve and compare to the current bounding box
		int comps_size = comps.size();
		for(int i = 0; i < comps_size; ++i) {

			QRectF cur_rect = getBoundingRect(comps[i].getNeutronSLDensities());

			if(nsld_rect.left() > cur_rect.left())
				nsld_rect.setLeft(cur_rect.left());

			if(nsld_rect.right() < cur_rect.right())
				nsld_rect.setRight(cur_rect.right());

			if(nsld_rect.bottom() > cur_rect.bottom())
				nsld_rect.setBottom(cur_rect.bottom());

			if(nsld_rect.top() < cur_rect.top())
				nsld_rect.setTop(cur_rect.top());
		}
	}

	return nsld_rect;
}

//Return the bounding rectangle for the xray form factor data
QRectF SimManager::xrayFormFactorRect() const {

	return getBoundingRect(xray_fourier_FF->getFormFactors());
}

//Return the bounding rectangle for the neutron form factor data
QRectF SimManager::neutronFormFactorRect() const {

	return getBoundingRect(neut_fourier_FF->getFormFactors());
}

//simple boolean for whether a component file is currently loaded
bool SimManager::componentFileLoaded() const {return compFileLoaded;}

//Clears the container of all Component objects
void SimManager::removeComponents() {

	reader.getComponents().clear();
	compFileLoaded = false;
	volCalculated = false;
}

//Returns true if the given particle has been loaded into the reader
bool SimManager::particleLoaded(const QString& name) const {

	return reader.particleLoaded(name);
}

//As advertised, populates the component derived matrix and vector that are use to calculate the volumes 
void SimManager::populateVolumeEquations() {

	//number of components
	int num_cmps = reader.getComponents().size()-1;

	//initialize the members with zeros
	cmp_matrix = QVector<Doubles> (num_cmps);
	for(int i = 0; i < num_cmps; ++i)
		cmp_matrix[i] = Doubles(num_cmps);

	cmp_vector = Doubles(num_cmps);

	double matrix_elem, vector_elem;

	Points nd1, nd2; //number densities from the components

	//matrix will be symmetric, so only calculate half
	for(int i = 0; i < num_cmps; ++i) {

		nd2 = reader.getComponents()[i].getNumberDensities();
		vector_elem = 0;

		int nd2_size = nd2.size();
		for(int k = 0; k < nd2_size; ++k)
			vector_elem += nd2[k].y();

		for(int j = i; j < num_cmps; ++j) {

			//reset sums for each loop
			matrix_elem = 0;

			nd1 = reader.getComponents()[j].getNumberDensities();


			//now need to loop through and get a sum of the products
			int nd1_size = nd1.size();
			for(int k = 0; k < nd1_size; ++k)
				matrix_elem += nd1[k].y() * nd2[k].y();

			//update the members
			cmp_matrix[i][j] = matrix_elem;
			cmp_matrix[j][i] = matrix_elem; //since symmetric
		}

		cmp_vector[i] = vector_elem;
	}
}

//Calculates the component volumes using a set of linear equations from SIMtoEXP paper (Eq. 5)
void SimManager::solveVolumeEquations() {

	//first need to populate the matrix and the vector
	populateVolumeEquations();
	Doubles volumes = lsolver::gauss_elim(cmp_matrix, cmp_vector);

	//go through the components and add the volumes (Total should not be added yet)
	int comp_size = reader.getComponents().size()-1;
	for(int i = 0; i < comp_size; ++i)
		reader.getComponents()[i].setVolume(volumes[i]);

	//with all volumes set, calculate the RMS for this solution
	calcVolumeRMS(volumes);

	//set boolean
	volCalculated = true;
}

//private helper function for calculating RMS of volume solution. Uses Eq. 6 from
//SIMtoEXP paper
void SimManager::calcVolumeRMS(const Doubles& volumes) {

	double sum_squares = 0;

	uint num_points = 0;
	if(numberComponents() > 0)
		num_points = reader.getComponents()[0].getNumberDensities().size();

	for(uint i = 0; i < num_points; ++i) {

		double prob_sum = 0;

		int vol_size = volumes.size();
		for(int j = 0; j < vol_size; ++j) {

			double prob = reader.getComponents()[j].getProbabilities().at(i).y();
			prob_sum += prob - 1;
		}

		sum_squares += prob_sum * prob_sum;
	}

	volume_rms = sqrt(sum_squares) / sqrt(num_points - volumes.size());
}

//Convenience function that go through and creates a plot for the TOTAL probability density for all components.
//Assumes that volumes are set!!  Will probably just get zero otherwise
QVector<QPointF> SimManager::getTotalComponentProbs() {

	Points total_probs;
	QList<Points> cmp_probs;

	//first go through and collect the functions
	int comp_size = reader.getComponents().size()-1;
	for(int i = 0; i < comp_size; ++i)
		cmp_probs.append(reader.getComponents()[i].getProbabilities());

	//should put an error check to ensure not empty
	int cmp_prob_size = cmp_probs[0].size();
	for(int i = 0; i < cmp_prob_size; ++i) {

		double cur_prob = 0;

		int cmp_probs_size = cmp_probs.size();
		for(int j = 0; j < cmp_probs_size; ++j)
			cur_prob += cmp_probs[j][i].y();

		total_probs.append(QPointF (cmp_probs[0][i].x(), cur_prob));
	}

	return total_probs;
}
