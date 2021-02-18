/*
 * ExpReader.cpp
 *
 *  Created on: Nov 15, 2012
 *      Author: bholland
 */

#include "ExpReader.h"

#include <QtCore/QFile>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QMessageBox>
#endif
#include <QtCore/QTextStream>
#include <QtCore/QStringList>

const unsigned int Z_COLUMN = 0;

//Default constructor
ExpReader::ExpReader() {}

//destructor
ExpReader::~ExpReader() {

	//use qt macro for deleting objects in QContainers
	qDeleteAll(xray_data);
	qDeleteAll(neutron_data);
}

//Return the ExpData object pointer at 'index' for the given ScatteringType
ExpData* ExpReader::getExpDataObject(uint index, ScatteringType type) {

	if(type == XRAY)
		return xray_data[index];

	else if(type == NEUTRON)
		return neutron_data[index];

	else
		return NULL;
}

//Return the ExpData object pointer at 'index' for the given ScatteringType, const version
const ExpData* ExpReader::getExpDataObject(uint index, ScatteringType type) const {

	if(type == XRAY)
		return xray_data[index];

	else if(type == NEUTRON)
		return neutron_data[index];

	else
		return NULL;
}

//Return the electron density at the given 'index'
QStringPointsPair& ExpReader::getElectronDensity(uint index) {

	return elec_densities[index];
}

const QStringPointsPair& ExpReader::getElectronDensity(uint index) const {

	return elec_densities[index];
}

//Return the neutron density at the given 'index'
QStringPointsPair& ExpReader::getNeutronSLDensity(uint index) {

	return neut_densities[index];
}

const QStringPointsPair& ExpReader::getNeutronSLDensity(uint index) const {

	return neut_densities[index];
}

//returns true if the appropriate container is empty of experimental data
bool ExpReader::formFactorEmpty(ScatteringType type) const {

	if(type == XRAY)
		return xray_data.empty();

	else if(type == NEUTRON)
		return neutron_data.empty();

	//if there are other ScatteringTypes, they must be empty since they don't have containers
	return true;
}

//Returns the number of data sets for the ScatteringType
uint ExpReader::numberDataSets(ScatteringType type) const {

	if(type == XRAY)
		return xray_data.size();

	else if(type == NEUTRON)
		return neutron_data.size();

	//otherwise return 0, since there is no container for 'type'
	return 0;
}

//Return the number of electron density sets loaded
uint ExpReader::numberElectronDensitySets() const {return elec_densities.size();}
uint ExpReader::numberNeutronSLDensitySets() const {return neut_densities.size();}

//Does the work of reading the experimental data file - MUST be in following format: 3 columns of
//values starting on the line after their titles.  All comment lines ('#') are ignored. ScatteringType
//determines which QList the data goes in
FileCode ExpReader::readExpFile(const QString& path, ScatteringType type) {

	//open the file
	QFile file(path);

	//try opening file, if doesn't work, output an error box
	if(not file.open(QIODevice::ReadOnly | QIODevice::Text)) {

		QMessageBox::critical(0, "Unknown file error:", "Cannot open file :-(");
		return FILE_ERROR;
	}

	//Read lines until "Data"
	QTextStream input(&file);
	bool pastLabels = false;

	ExpData* data = new ExpData; //for storing data

	while(not input.atEnd()) {

		//get the next line
		QString line = input.readLine();

		//check the first character, if it's a '#', then skip to next line
		if(line[0] == '#')
			continue;

		//skip the first line after the comments
		if(not pastLabels) {
			pastLabels = true;
			continue;
		}

		//otherwise start to parse and store the data
		QStringList list = line.split(QRegExp("\\s"), QString::SkipEmptyParts); //split the line using any 'non-word' character

		//if the list is not 2 or 3 columns wide, there is a problem
		if(list.length() < 2 or list.length() > 3) {

			//show message, clean up and return an error
			QMessageBox::critical(0, "Input error:", "One or more of the input values is invalid.");
			delete data;
			return FILE_ERROR;
		}

		//'q' should be in first column, form factor in the next, then uncertainty
		bool q_ok, FF_ok, err_ok = true; //for error check
		double q = list[0].toDouble(&q_ok);
		double FF = list[1].toDouble(&FF_ok);
		double err = 0;

		//if there is error included, use it
		if(list.length() == 3)
			err = list[2].toDouble(&err_ok);

		//if any problems with the data, show message and exit input
		if(not(q_ok and FF_ok and err_ok)) {

			//show message, clean up and return an error
			QMessageBox::critical(0, "Input error:", "One or more of the input values is invalid.");
			delete data;
			return FILE_ERROR;
		}

		data->addPoint(q,FF,err);
	}

	if(type == XRAY)
		xray_data.append(data);

	else if(type == NEUTRON)
		neutron_data.append(data);

	//clean up
	file.close();

	return READ_SUCCESS;
}

//Reads in a data file of experimental electron or neutron scattering length densities, depending on
//ScatteringType, although both file formats should be the same.
FileCode ExpReader::readExpDensityFile(const QString& path, ScatteringType type) {

	QStringPointsPairs* cur_container;

	if(type == XRAY)
		cur_container = &elec_densities;
	else if(type == NEUTRON)
		cur_container = &neut_densities;

	//open the file
	QFile file(path);

	//try opening file, if doesn't work, output an error box
	if(not file.open(QIODevice::ReadOnly | QIODevice::Text)) {

		QMessageBox::critical(0, "Unknown file error:", "Cannot open file :-(");
		return FILE_ERROR;
	}

	//Read lines until "Data"
	QTextStream input(&file);

	//now read data into map container
	bool first_row = true;
	while(not input.atEnd()) {

		QString line = input.readLine();

		//check the first character, if it's a '#', then skip to next line
		if(line[0] == '#')
			continue;

		QStringList list = line.split(QRegExp("\\s"), QString::SkipEmptyParts); //split the line using any 'non-word' character

		if(first_row) {

			//skip 'z(A)' in first column
			for(int i =  1; i < list.length(); i++)
				cur_container->append(QPair<QString,QVector<QPointF> >(list[i],QVector<QPointF>())); //initialize each pair in the container

			first_row = false;

		} else {

			//now go through and for each line, add the values to the vectors
			for(int i =  1; i < list.length(); i++) {

				double x = list[Z_COLUMN].toDouble();
				double y = list[i].toDouble();

				//add the point to the appropriate pair
				(*cur_container)[i-1].second.push_back(QPointF (x,y));
			}
		}
	}

	//clean up
	file.close();

	return READ_SUCCESS;
}

//Remove the data at 'index' for the given scattering type
void ExpReader::removeData(int index, ScatteringType type) {

	if(type == XRAY) {

		delete xray_data[index]; //first delete the object
		xray_data.removeAt(index); //then remove pointer
	}
	else if (type == NEUTRON) {

		delete neutron_data[index]; //first delete the object
		neutron_data.removeAt(index); //then remove pointer
	}
}

//Clears all electron density data
void ExpReader::removeElectronDensityData() {

	//simply clear the container
	elec_densities.clear();
}

//Clears all neutron density data
void ExpReader::removeNeutronSLDensityData() {

	neut_densities.clear();
}
