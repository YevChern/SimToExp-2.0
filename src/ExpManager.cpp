
#include "ExpManager.h"
#include "KillBox.h"
#include "Util.h"

#include <cfloat>

#include <QtCore/QDebug>

using namespace std;

//Default constructor
ExpManager::ExpManager(QWidget* parent) : QWidget(parent) {

	//initialize booleans
	ed_data_loaded = false;
	nsld_data_loaded = false;
}

//Destructor
ExpManager::~ExpManager() {}

//Getters for data
ExpData& ExpManager::getXrayDataObject(uint index) {

	//pass job along to reader
	return *(reader.getExpDataObject(index, XRAY));
}

const ExpData& ExpManager::getXrayDataObject(uint index) const {

	//pass job along to reader
	return *(reader.getExpDataObject(index, XRAY));
}

ExpData& ExpManager::getNeutronDataObject(uint index) {

	//pass job along to reader
	return *(reader.getExpDataObject(index, NEUTRON));
}

const ExpData& ExpManager::getNeutronDataObject(uint index) const {

	//pass job along to reader
	return *(reader.getExpDataObject(index, NEUTRON));
}

Points& ExpManager::getXrayData(uint index) {

	return reader.getExpDataObject(index,XRAY)->getData();
}

const Points& ExpManager::getXrayData(uint index) const {

	return reader.getExpDataObject(index,XRAY)->getData();
}

Points& ExpManager::getNeutronData(uint index) {

	return reader.getExpDataObject(index,NEUTRON)->getData();
}

const Points& ExpManager::getNeutronData(uint index) const {

	return reader.getExpDataObject(index,NEUTRON)->getData();
}

QStringPointsPair& ExpManager::getElectronDensityData(uint index) {

	return reader.getElectronDensity(index);
}

const QStringPointsPair& ExpManager::getElectronDensityData(uint index) const {

	return reader.getElectronDensity(index);
}

QStringPointsPair& ExpManager::getNeutronSLDensityData(uint index) {

	return reader.getNeutronSLDensity(index);
}

const QStringPointsPair& ExpManager::getNeutronSLDensityData(uint index) const {

	return reader.getNeutronSLDensity(index);
}

//Return the scaling factor for sample 'index' and ScatteringType 'type'
double ExpManager::getScalingFactor(uint index, ScatteringType type) const {

	if(type == XRAY)
		return reader.xray_data[index]->getScalingFactor();

	else if(type == NEUTRON)
		return reader.neutron_data[index]->getScalingFactor();

	//otherwise return -1 as an error
	return -1;
}

//Return the chi squared value for sample 'index' and ScatteringType 'type'
double ExpManager::getChiSquared(uint index, ScatteringType type) const {

	if(type == XRAY)
		return reader.xray_data[index]->getChiSquared();

	else if(type == NEUTRON)
		return reader.neutron_data[index]->getChiSquared();

	//otherwise return -1 as error
	return -1;
}

//getters that return the number of open data sets
uint ExpManager::getNumberXray() const {return reader.numberDataSets(XRAY);}
uint ExpManager::getNumberNeutron() const {return reader.numberDataSets(NEUTRON);}

//getters that return the number of components stored
uint ExpManager::getNumberElectronDensity() const {return reader.numberElectronDensitySets();}
uint ExpManager::getNumberNeutronSLDensity() const {return reader.numberNeutronSLDensitySets();}

//getters for booleans
bool ExpManager::electronDensityLoaded() const {return ed_data_loaded;}
bool ExpManager::neutronSLDensityLoaded() const {return nsld_data_loaded;}

bool ExpManager::xrayFormFactorEmpty() const {
	return reader.formFactorEmpty(XRAY);
}

bool ExpManager::neutronFormFactorEmpty() const {
	return reader.formFactorEmpty(NEUTRON);
}

//Deletes the data set at the given index
void ExpManager::removeXrayData(uint index) {

	//pass along to reader
	reader.removeData(index, XRAY);
}

void ExpManager::removeNeutronData(uint index) {

	//pass along to reader
	reader.removeData(index, NEUTRON);
}

//Clears all of the data from the loaded experimental scattering density containers
void ExpManager::removeElectronDensityData() {

	reader.removeElectronDensityData();
	ed_data_loaded = false;
}

void ExpManager::removeNeutronSLDensityData() {

	reader.removeNeutronSLDensityData();
	nsld_data_loaded = false;
}

//MUST be in following format: 3 columns of values starting on the line after their titles.
//All comment lines ('#') are ignored. ScatteringType determines which QList the data goes in
FileCode ExpManager::readExpFile(const QString& path, ScatteringType type) {

	//pass work along to reader
	return reader.readExpFile(path, type);
}

//Reads in a data file of experimental electron or neutron scattering length densities, depending on
//ScatteringType, although both file formats should be the same.
FileCode ExpManager::readExpDensityFile(const QString& path, ScatteringType type) {

	FileCode code;
	code = reader.readExpDensityFile(path, type);

	if(code == READ_SUCCESS) {

		//set booleans if successful
		if(type == XRAY) {

			ed_data_loaded = true;
		}
		else if(type == NEUTRON) {

			nsld_data_loaded = true;
		}
	}

	return code;
}

//Scales the data automatically by calculating the scaling factor using Eq. 2 of SIMtoEXP paper.  Emits both 'dataScaled'
//signal after finished
void ExpManager::scaleData(const Points& sim_xray, const Points& sim_neutron) {

	//go through each list and scale; if nothing in list, then nothing is done
	for(int i = 0; i < reader.xray_data.size(); i++) {

		//first reset in case scaling has been done manually
		if(reader.xray_data[i]->getScalingFactor() != 1.0)
			reader.xray_data[i]->resetData();

		reader.xray_data[i]->calcScalingFactor(sim_xray);
	}

	for(int i = 0; i < reader.neutron_data.size(); i++) {

		//first reset in case scaling has been done manually
		if(reader.neutron_data[i]->getScalingFactor() != 1.0)
			reader.neutron_data[i]->resetData();

		reader.neutron_data[i]->calcScalingFactor(sim_neutron);
	}

	emit xrayDataScaled();
	emit neutronDataScaled();
}

//Scales the data automatically by calculating the scaling factor using Eq. 2 of SIMtoEXP paper.  Emits both 'dataScaled'
//signal after finished
void ExpManager::scaleData(const Points& sim, ScatteringType type) {

	//go through each list and scale; if nothing in list, then nothing is done
	if(type == XRAY) {

		for(int i = 0; i < reader.xray_data.size(); i++) {

			reader.xray_data[i]->calcScalingFactor(sim);
		}

		emit xrayDataScaled();
	}

	else if(type == NEUTRON) {

		for(int i = 0; i < reader.neutron_data.size(); i++) {

			reader.neutron_data[i]->calcScalingFactor(sim);
		}

		emit neutronDataScaled();
	}
}

//Scales the data using the scaling factors provided.  Emits 'dataScaled' signal when finished
void ExpManager::scaleXrayData(const Doubles& factors) {

	//vectors must be the same size, otherwise kill the application with a QKillBox (in lieu of a Seg fault)
	if(factors.size() != reader.xray_data.size())
		KillBox (parentWidget(), "Vectors do not match: must have same number of scaling factors as plots. Function: ExpParser::scaleXrayData()");


	//go through each list and scale; if nothing in list, then nothing is done
	int xray_size = reader.xray_data.size();
	for(int i = 0; i < xray_size; ++i) {

		//do nothing if scaling factor hasn't changed
		if(factors[i] != reader.xray_data[i]->getScalingFactor()) {

			reader.xray_data[i]->scaleData(factors[i]);
		}
	}

	emit xrayDataScaled();
}

//Scales the data using the scaling factors provided.  Emits 'dataScaled' signal when finished
void ExpManager::scaleNeutronData(const Doubles& factors) {

	//vectors must be the same size, otherwise kill the application with a QKillBox (in lieu of a Seg fault)
	if(factors.size() != reader.neutron_data.size())
		KillBox (parentWidget(), "Vectors do not match: must have same number of scaling factors as plots. Function: ExpParser::scaleNeutronData()");


	//go through each list and scale; if nothing in list, then nothing is done
	int neutron_size = reader.neutron_data.size();
	for(int i = 0; i < neutron_size; ++i) {

		//do nothing if scaling factor is one
		if(factors[i] != reader.neutron_data[i]->getScalingFactor()) {

			reader.neutron_data[i]->scaleData(factors[i]);
		}
	}

	emit neutronDataScaled();
}

//Return the bounding rectangle for the electron density data by searching through the curves
const QRectF& ExpManager::electronDensityRect() {

	//loop through, get the rectangle for each curve and compare to the current bounding box
	for(uint i = 0; i < reader.numberElectronDensitySets(); i++) {

		QRectF cur_rect = qutil::getBoundingRect(reader.getElectronDensity(i).second);

		//need to initialize the rect as (0,0) may not even be in final rectangle
		if(i == 0)
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

	return ed_rect;
}

//Return the bounding rectangle for the neutron scattering length density data by searching through the curves
const QRectF& ExpManager::neutronSLDensityRect() {

	//loop through, get the rectangle for each curve and compare to the current bounding box
	for(uint i = 0; i < reader.numberNeutronSLDensitySets(); i++) {

		QRectF cur_rect = qutil::getBoundingRect(reader.getNeutronSLDensity(i).second);

		//need to initialize the rect as (0,0) may not even be in final rectangle
		if(i == 0)
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

	return nsld_rect;
}

//Return the bounding rectangle for the xray form factor data
const QRectF& ExpManager::xrayFormFactorRect() {

	//loop through, get the rectangle for each curve and compare to the current bounding box
	for(uint i = 0; i < reader.numberDataSets(XRAY); i++) {

		QRectF cur_rect = qutil::getBoundingRect(reader.getExpDataObject(i, XRAY)->getData());

		//need to initialize the rect since (0,0) may not even be in final rectangle
		if(i == 0)
			xrayFF_rect = cur_rect;

		else {

			if(xrayFF_rect.left() > cur_rect.left())
				xrayFF_rect.setLeft(cur_rect.left());

			if(xrayFF_rect.right() < cur_rect.right())
				xrayFF_rect.setRight(cur_rect.right());

			if(xrayFF_rect.bottom() > cur_rect.bottom())
				xrayFF_rect.setBottom(cur_rect.bottom());

			if(xrayFF_rect.top() < cur_rect.top())
				xrayFF_rect.setTop(cur_rect.top());
		}
	}

	return xrayFF_rect;
}

//Return the bounding rectangle for the neutron form factor data
const QRectF& ExpManager::neutronFormFactorRect() {

	//loop through, get the rectangle for each curve and compare to the current bounding box
	for(uint i = 0; i < reader.numberDataSets(NEUTRON); i++) {

		QRectF cur_rect = qutil::getBoundingRect(reader.getExpDataObject(i, NEUTRON)->getData());

		//need to initialize the rect since (0,0) may not even be in final rectangle
		if(i == 0)
			neutFF_rect = cur_rect;

		else {

			if(neutFF_rect.left() > cur_rect.left())
				neutFF_rect.setLeft(cur_rect.left());

			if(neutFF_rect.right() < cur_rect.right())
				neutFF_rect.setRight(cur_rect.right());

			if(neutFF_rect.bottom() > cur_rect.bottom())
				neutFF_rect.setBottom(cur_rect.bottom());

			if(neutFF_rect.top() < cur_rect.top())
				neutFF_rect.setTop(cur_rect.top());
		}
	}

	return neutFF_rect;
}

//Resets all of the ExpData objects with respect to their scaling.
//Ex: the simulation data has been deleted, so reset all scaling
void ExpManager::resetAllData() {

	int num_xrays = reader.xray_data.size();
	int num_neutron = reader.neutron_data.size();

	for(int i = 0; i < num_xrays; ++i)
		reader.xray_data[i]->resetData();

	for(int i = 0; i < num_neutron; ++i)
		reader.neutron_data[i]->resetData();

	//emit signal so plot can be updated
	emit xrayDataScaled();
	emit neutronDataScaled();
}
