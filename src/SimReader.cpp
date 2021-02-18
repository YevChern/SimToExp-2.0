/*
 * SimReader.cpp
 *
 *  Created on: Nov 14, 2012
 *      Author: bholland
 */

#include "SimReader.h"

#include <cmath>
#include <iostream>

#include <QtCore/QFile>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QMessageBox>
#endif
#include <QtCore/QTextStream>
#include <QtCore/QtDebug>

using namespace std;

//simulation input file parameters
const QString DATA = "data";
const QString AREA = "area";
const QString NLIP = "Nlip";
const QString NORM = "norm";

const unsigned int Z_COLUMN = 0; //column that the 'z' axis is located in the sim file

//Default constructor
SimReader::SimReader() {

	//initialize
	area = 1.0;
	nlip = 1;
	normalize_volume = false;
	thickness = 0;
}

//Destructor
SimReader::~SimReader() {}

//getters
double SimReader::lipidArea() const {

	return area;
}

uint SimReader::numLipids() const {

	return nlip;
}

//return the PD object if there, a default constructed NormalDensity object otherwise (see QMap documentation)
NormalDensity& SimReader::getParticleDensity(const QString& name) {

	//get the value and check to see if it was constructed by QMap
	NormalDensity& density = sim_data[name];

	//if it was constructed by QMap (i.e. not in the container), output message to cerr
	if(density.densitiesAtEnd())
		cerr << "Warning: name not found in map, default NormalDensity created and added with name '" << name.toStdString() << "'" << endl;

	return density;
}

const NormalDensity& SimReader::getParticleDensity(const QString& name) const {

	//get the value and check to see if it was constructed by QMap
	const NormalDensity& density = sim_data[name];

	//if it was constructed by QMap (i.e. not in the container), output message to cerr
	if(density.isEmpty())
		cerr << "Warning: name likely not found in map, default NormalDensity created and added with name '" << name.toStdString() << "'" << endl;

	return density;
}

//Returns a QList of the NormalDensity objects stored in the map, but ordered according to
//order read in
NormalDensities SimReader::getParticleDensities() const {

	return sim_data.values();
}

//Removes all data from the NormalDensity map and reset the other attributes
void SimReader::clearSimData() {

	sim_data.clear();
	area = 0;
	nlip = 0;
	normalize_volume = false;
	thickness = 0;
	first_line.clear();
	components.clear();
}

//Return the first NormalDensity in the map and set the iterator to the beginning
NormalDensity* SimReader::getFirstParticleDensity() {

	sim_data_itr = sim_data.begin();

	//if at end, return an error through cerr for debugging
	if(sim_data_itr == sim_data.end())
		cerr << "SimReader::getFirstParticleDensity: NormalDensity map is empty, no particle to retrieve" << endl;

	//set result to a variable to pass back, then increment the iterator for next time
	//Need to increment here or can't tell if at end
	return &sim_data_itr.value();
}

//Return the next NormalDensity in the map; does not initialize the iterator to beginning!
NormalDensity* SimReader::getNextParticleDensity() {

	//first update
	++sim_data_itr;

	//if at end, return a dummy object, should not happen in a properly written for loop
	if(sim_data_itr == sim_data.end())
		return &dummy_density;

	//Iterator is already at the "next" density, so just return the current one and increment
	return &sim_data_itr.value();
}

//Boolean that checks if there are any NormalDensity objects in the map
bool SimReader::isEmpty() const {return sim_data.empty();}

//Checks to see if iterator at end of container - should be used in conjunction with the iterative methods
bool SimReader::atEnd() const {

	if(sim_data_itr == sim_data.end())
		return true;
	return false;
}

//Checks to see if the given string is in the map
bool SimReader::particleLoaded(const QString& name) const {

	if(sim_data.find(name) == sim_data.end())
		return false;
	return true;
}

//getters for the components
Components& SimReader::getComponents() {

	return components;
}

const Components& SimReader::getComponents() const {

	return components;
}

//Opens and reads the data in the file given by 'path'. If any error occurs,
//simply returns (giving hopefully useful error message).
FileCode SimReader::readSimFile(const QString& path) {

	//open the file
	QFile file(path);

	//try opening file, if doesn't work, output an error box
	if(not file.open(QIODevice::ReadOnly | QIODevice::Text)) {

		QMessageBox::critical(0, "Unknown file error:", "Cannot open file.");
		return FILE_ERROR;
	}

	//Read lines until "Data"
	QTextStream input(&file);
	while(not input.atEnd()) {

		//boolean for checking conversions
		bool ok = true;

		//get the next line
		QString line = input.readLine();

		//check the first character, if it's a '#', then skip to next line
		if(line[0] == '#')
			continue;

		//check to see if at 'data'
		if(line.startsWith(DATA, Qt::CaseInsensitive))
			break;

		//otherwise, must be one of input parameters
		if(line.startsWith(AREA, Qt::CaseInsensitive)) {

			line.remove(0, AREA.length()+1);
			area = line.toDouble(&ok);

		} else if(line.startsWith(NLIP, Qt::CaseInsensitive)) {

			line.remove(0, NLIP.length()+1);
			nlip = line.toUInt(&ok);

		} else if(line.startsWith(NORM, Qt::CaseInsensitive)) {

			line.remove(0, NORM.length()+1);
			if(line.toLower() == "yes")
				normalize_volume = true;
			else if(line.toLower() == "no")
				normalize_volume = false;
			else
				ok = false;

		} else {

			//if not one of the parameters, there is an error in the file
			QMessageBox::critical(0, "Unknown parameter:", "At least one parameter in the simulation file is not recognized by SIMtoEXP.");
			qDebug() << "Sim file input error: current string = " << line;
			return FILE_ERROR;
		}

		//if conversion was unsuccessful, return message box and file error
		if(not ok) {
			QMessageBox::critical(0, "Input error:", "One or more parameters has an invalid value.");
			return FILE_ERROR;
		}
	}

	//now read data into container
	bool first_row = true;
	while(not input.atEnd()) {

		QString line = input.readLine();
		QStringList list = line.split(QRegExp("\\s"), QString::SkipEmptyParts); //split the line using any 'non-word' character

		if(first_row) {

			first_line = list;

			//skip 'z(A)' in first column, so start at i = 1
			uint length = list.length();
			for(uint i =  1; i < length; ++i) {

				QStringList cur_names = list.at(i).split('/');

				//should only be one or two, if more there is an error
				if(cur_names.length() < 1 or cur_names.length() > 2) {

					QString msg = "There is an error with in the particle/residue name combination: ";
					msg += list.at(i);
					msg += ".";
					QMessageBox::critical(0, "Input error:", msg);
					return FILE_ERROR;
				}

				QString pName = cur_names[0]; //even if no residue, still the particle name
				QString rName = "";
				if(cur_names.length() == 2)
					rName = cur_names[1];

				//make a PD object with the full particle name from the input file, also keep track of position
				ParticleLabel label(rName, pName);
				sim_data.insert(list.at(i), NormalDensity(label));
			}

			first_row = false;

		} else {

			//now go through and for each line, add the values to the vectors
			int length = list.length();
			for(int i =  1; i < length; ++i) {

				const QString& cur_name = first_line[i];
				bool x_ok, y_ok;
				double x = list[Z_COLUMN].toDouble(&x_ok);
				double y = list[i].toDouble(&y_ok);

				//if either conversion unsuccessful or a probability value is less than zero, show box and return
				if(not(x_ok and y_ok) or y_ok < 0) {

					QString msg = "A density value for particle ";
					msg += sim_data.value(cur_name).particleName();
					msg += " is invalid.";
					QMessageBox::critical(0, "Input error:", msg);
					return FILE_ERROR;
				}

				sim_data[cur_name].addDensity(x,y);
			}
		}
	}

	//now all read in, normalize the data if necessary
	if(normalize_volume)
		normalizeData();

	//clean up
	file.close();

	//return successful
	return READ_SUCCESS;
}

FileCode SimReader::readComponentFile(const QString& path) {

	//open the file
	QFile file(path);

	//try opening file, if doesn't work, output an error box
	if(not file.open(QIODevice::ReadOnly | QIODevice::Text)) {

		QMessageBox::critical(0, "Unknown file error:", "Cannot open file :-(");
		return FILE_ERROR;
	}

	//Have one component for 'Total' values
	Component total_comp("Total", 1);

	//Read lines until "Data"
	QTextStream input(&file);

	while(not input.atEnd()) {

		QString line = input.readLine();

		//check the first character, if it's a '#', then skip to next line
		if(line[0] == '#')
			continue;

		QStringList list = line.split(QRegExp("\\s"), QString::SkipEmptyParts); //split the line using any 'non-word' character

		//If list is not at least three long, there is an error
		if(list.size() < 3) {

			QMessageBox::critical(0, "Input error:", "There must be a component name, followed by the number of groups, followed by at least one atom.");
			components.clear();
			return FILE_ERROR;
		}

		//First is the component name, then the number of groups
		bool ok; //for conversion check
		QString name = list[0];
		double num_groups = list[1].toDouble(&ok);

		//Error check for group number entry
		if(not ok or num_groups <= 0) {

			QMessageBox::critical(0, "Input error:", "'Number of groups' must be a number greater than zero");
			components.clear();
			return FILE_ERROR;
		}

		//Create a Component object
		Component cur_comp(name, num_groups);

		//add the atoms to the component
		int length = list.size();
		for(int i = 2; i < length; i++) {

			//if name not found, there is an error in the file
			NormalDensity& pd = getParticleDensity(list[i]);
			//object should have values in it from sim file
			if(pd.isEmpty()) {

				QString msg = "The particle id ";
				msg += list[i];
				msg += " has not been loaded or has not data. Possibly a typo?";
				QMessageBox::critical(0, "Input error:", msg);

				components.clear();
				return FILE_ERROR;
			}

			//otherwise, add it to components
			cur_comp.addParticle(pd);
			total_comp.addParticle(pd);
		}

		//add the component to the list and calculate densities
		cur_comp.calcDensities();
		components.append(cur_comp);
	}

	total_comp.calcDensities();
	components.append(total_comp);

	//clean up
	file.close();

	//return successful
	return READ_SUCCESS;
}

//private function that goes through and normalizes the data from a histogram style input to a set of probability values
void SimReader::normalizeData() {

	//if nothing read yet, do nothing
	if(sim_data.empty())
		return;

	//otherwise go through the ParticleDensity objects one by one and make the appropriate adjustments
	NormalDensityMap::iterator itr, itr_end;
	itr_end = sim_data.end();
	for(itr = sim_data.begin(); itr != itr_end; ++itr) {

		Points points = itr.value().densities();

		//calculate the bin thickness
		if(thickness == 0) {

			if(points.size() >= 2)
				thickness =  fabs(points[0].x() - points[1].x());

			else {

				QMessageBox::critical(0, "Input error:", "There must be at least two histogram bins.");
				return;
			}
		}

		int length = points.size() - 1;
		for(int j = 0; j < length; ++j) {

			double prob = points[j].y() / (nlip * area * thickness);
			points[j].setY(prob);
		}
	}
}
