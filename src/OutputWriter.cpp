/*
 * OutputWriter.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: bholland
 */

#include "OutputWriter.h"

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QFileDialog>
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QFileDialog>
    #include <QtGui/QMessageBox>
#endif

#include <QtCore/QVector>

#include <cmath>

using namespace std;

const uint OUTPUT_PRECISION = 6;
const uint Z_PRECISION = 5;

const QString DEFAULT_FILENAME = "simtoexp_untitled";
const QString NUMBER_DENSITY_FILE_SUFFIX = ".sim";
const QString ELECTRON_DENSITY_FILE_SUFFIX = ".edp";
const QString NEUTRON_DENSITY_FILE_SUFFIX = ".nsld";
const QString XRAY_STRUCTURE_FACTOR_SUFFIX = ".xsf";
const QString NEUTRON_STRUCTURE_FACTOR_SUFFIX = ".nsf";
const QString VOLUMETRIC_FIT_FILE_SUFFIX = ".cvf";

//Default constructor
OutputWriter::OutputWriter(QWidget* parent) {

	this->parent = parent;
	setHeader();

	out_stream.precision(OUTPUT_PRECISION);
}

//Destructor
OutputWriter::~OutputWriter() {}

//Setter for the parent window/widget
void OutputWriter::setParent(QWidget* parent) {

	this->parent = parent;
}

/*
 * Write the currently loaded number densities into the SIMtoEXP format, that can then be reloaded in the future.
 * i.e. save the loaded data, especially if calculated from a trajectory.
 */
void OutputWriter::writeNumberDensities(const NormalDensities& densities, uint num_lipids, double area) {

	QString file_path = getFilePath(NUMBER_DENSITY);

	//if the returned string is good, open the output stream
	if(not (file_path.isNull() or file_path.isEmpty())) {

		//append the suffix if it isn't already there
		if(not file_path.endsWith(NUMBER_DENSITY_FILE_SUFFIX, Qt::CaseInsensitive))
			file_path.append(NUMBER_DENSITY_FILE_SUFFIX);

		out_stream.open(file_path.toStdString());

		//first output header and commented area
		out_stream << header.toStdString() << endl;
		out_stream << "# Any comments can go here, as long as they start with '#'" << endl;
		out_stream << "#" << endl;
		out_stream << "# Automatically generated data file, change values at your own peril!" << endl;
		out_stream << "#" << endl;
		out_stream << "# Average total area of lipids (Ang^2)" << endl;
		out_stream << "area " << area << endl;
		out_stream << "#" << endl;
		out_stream << "#" << endl;
		out_stream << "# Number of lipid molecules" << endl;
		out_stream << "Nlip " << num_lipids << endl;
		out_stream << "#" << endl;
		out_stream << "#" << endl;
		out_stream << "# Area per lipid (Ang^2) - calculated assuming a bilayer!" << endl;
		out_stream << "# " << 2*area/num_lipids << endl;
		out_stream << "#" << endl;
		out_stream << "#" << endl;
		out_stream << "#***************************************************************************************" << endl;
		out_stream << "Data" << endl;

		//get all of the densities and put them in a container
		QVector<Points> data;
		for(auto density : densities) {

			data.append(density.densities());
		}

		writeDensityData(data, densities);
	}
}

void OutputWriter::writeComponentNumberDensities(const Components& components) {

	//pass the work along to private function
	writeComponents(NUMBER_DENSITY, components);
}

void OutputWriter::writeElectronDensities(const NormalDensities& densities) {

	QString file_path = getFilePath(ELECTRON_DENSITY);

	//if the returned string is good, open the output stream
	if(not (file_path.isNull() or file_path.isEmpty())) {

		//append the suffix if it isn't already there
		if(not file_path.endsWith(ELECTRON_DENSITY_FILE_SUFFIX, Qt::CaseInsensitive))
			file_path.append(ELECTRON_DENSITY_FILE_SUFFIX);

		out_stream.open(file_path.toStdString());

		//first output header and commented area
		out_stream << header.toStdString() << endl;
		out_stream << "# Electron Density Profile (EDP)" << endl;
		out_stream << "#***************************************************************************************" << endl;

		//get all of the densities and put them in a container
		QVector<Points> data;
		for(auto density : densities) {

			data.append(density.electronDensities());
		}

		writeDensityData(data, densities);
	}
}

void OutputWriter::writeComponentElectronDensities(const Components& components) {

	//pass the work along to private function
	writeComponents(ELECTRON_DENSITY, components);
}

void OutputWriter::writeNeutronDensities(const NormalDensities& densities) {

	QString file_path = getFilePath(NEUTRON_DENSITY);

	//if the returned string is good, open the output stream
	if(not (file_path.isNull() or file_path.isEmpty())) {

		//append the suffix if it isn't already there
		if(not file_path.endsWith(NEUTRON_DENSITY_FILE_SUFFIX, Qt::CaseInsensitive))
			file_path.append(NEUTRON_DENSITY_FILE_SUFFIX);

		out_stream.open(file_path.toStdString());

		//first output header and commented area
		out_stream << header.toStdString() << endl;
		out_stream << "# Neutron Scattering Length Distribution (NSLD)" << endl;
		out_stream << "#***************************************************************************************" << endl;

		//get all of the densities and put them in a container
		QVector<Points> data;
		for(auto density : densities) {

			data.append(density.neutronDensities());
		}

		writeDensityData(data, densities);
	}
}

void OutputWriter::writeComponentNeutronDensities(const Components& components) {

	//pass the work along to private function
	writeComponents(NEUTRON_DENSITY, components);
}

void OutputWriter::writeXrayStructureFactors(const Points& data) {

	QString file_path = getFilePath(XRAY_STRUCTURE_FACTOR);

	//if the returned string is good, open the output stream
	if(not (file_path.isNull() or file_path.isEmpty())) {

		//append the suffix if it isn't already there
		if(not file_path.endsWith(XRAY_STRUCTURE_FACTOR_SUFFIX, Qt::CaseInsensitive))
			file_path.append(XRAY_STRUCTURE_FACTOR_SUFFIX);

		out_stream.open(file_path.toStdString());

		//first output header and commented area
		out_stream << header.toStdString() << endl;
		out_stream << "# Xray Structure Factor calculated from simulation (XSF)" << endl;
		out_stream << "#***************************************************************************************" << endl;

		writeStructureFactorData(data);
	}
}

void OutputWriter::writeNeutronStructureFactors(const Points& data) {

	QString file_path = getFilePath(NEUTRON_STRUCTURE_FACTOR);

	//if the returned string is good, open the output stream
	if(not (file_path.isNull() or file_path.isEmpty())) {

		//append the suffix if it isn't already there
		if(not file_path.endsWith(NEUTRON_STRUCTURE_FACTOR_SUFFIX, Qt::CaseInsensitive))
			file_path.append(NEUTRON_STRUCTURE_FACTOR_SUFFIX);

		out_stream.open(file_path.toStdString());

		//first output header and commented area
		out_stream << header.toStdString() << endl;
		out_stream << "# Neutron Structure Factor calculated from simulation (NSF)" << endl;
		out_stream << "#***************************************************************************************" << endl;

		writeStructureFactorData(data);
	}
}

/*
 * Private function that dose the work of writing the data - all output formats are essentially the same
 * Just columns of data for each atom/particle type as a function of z
 */
void OutputWriter::writeDensityData(const QVector<Points>& data, const NormalDensities& densities) {

	//then the data, first line is the particle names.
	//Note: total space needed for z is: sign + ceiling(log(x)) + precision + period
	double order_mag = ceil(log10(abs(data[0][0].x())));
	uint z_column_size = Z_PRECISION + 2 + order_mag;

	//Note: total space needed for densities is precision + "E-?" + period
	uint density_column_size = OUTPUT_PRECISION + 6;

	//output a row for the names / column titles
	QString z_title = "z(A)";
	out_stream << addSpaces(z_title, z_column_size).toStdString();
	for(const auto& density : densities) {

		//if there is an underscore in the particle name, this makes it a MARTINI particle
		//only want to keep the last portion of the name though, res part is redundant
		QString particle_name = density.particleName();
		QStringList pName_list = particle_name.split('_');
		particle_name = pName_list.last() + "/" + density.residueName();

		out_stream << addSpaces(particle_name, density_column_size).toStdString();
	}
	out_stream << endl;

	writeData(data, z_column_size);

	//writing done, close stream
	out_stream.close();
}

/*
 * Private function that writes out all exported component data, as components use a different
 * class then NormalDensity.
 */
void OutputWriter::writeComponentData(const QVector<Points>& data, const Components& components) {

	//then the data, first line is the particle names.
	//Note: total space needed for z is: sign + ceiling(log(x)) + precision + period
	double order_mag = ceil(log10(abs(data[0][0].x())));
	uint z_column_size = Z_PRECISION + 2 + order_mag;

	//Note: total space needed for densities is precision + "E-?" + period
	uint column_size = OUTPUT_PRECISION + 6;

	//output a row for the names / column titles
	QString z_title = "z(A)";
	out_stream << addSpaces(z_title, z_column_size).toStdString();
	for(const auto& component : components) {

		//Components aren't necessarily part of any residue (could be a mixture), so only
		//output the name - if they can't be loaded back, not going to support this anyway
		out_stream << addSpaces(component.name(), column_size).toStdString();
	}
	out_stream << endl;

	writeData(data, z_column_size);

	//writing done, close stream
	out_stream.close();
}

/*
 * Writes out the actual numerical data; intended to be used after a row of titles has been
 * printed out. DOES NOT close the file, leaves this up to the calling function.
 */
void OutputWriter::writeData(const QVector<Points>& data, double z_column_size) {

	//now output the actual data

	//now go through this vector 'length' times for each row
	uint num_rows = data[0].size();
	uint num_columns = data.size();

	for(uint row = 0; row < num_rows; ++row) {

		//output the 'z' column
		QString value;
		value.setNum(data[0][row].x(), 'f', Z_PRECISION);
		out_stream << addSpaces(value, z_column_size).toStdString();

		for(uint column = 0; column < num_columns; ++column) {

			QString value;
			value.setNum(data[column][row].y(), 'E', OUTPUT_PRECISION);
			value.append("  ");
			out_stream << value.toStdString();
		}

		out_stream << endl;
	}
}

/*
 * Private function that dose the work of writing the structure factor data - both output formats (xray and neutron) are the same
 * Just columns of data for each atom/particle type as a function of z
 */
void OutputWriter::writeStructureFactorData(const Points& data) {

	//then the data, first line is the particle names.
	//Note: total space needed for z is: sign + ceiling(log(x)) + precision + period
	double order_mag = ceil(log10(abs(data[0].x())));
	uint z_column_size = Z_PRECISION + 2 + order_mag;

	//output a row for the names / column titles
	QString z_title = "z(A)";
	out_stream << addSpaces(z_title, z_column_size).toStdString();
	out_stream << "S(q)" << endl;

	//now output the actual data

	//now go through this vector 'length' times for each row
	uint num_rows = data.size();

	for(uint row = 0; row < num_rows; ++row) {

		//output the 'z' column
		QString value;
		value.setNum(data[row].x(), 'f', Z_PRECISION);
		out_stream << addSpaces(value, z_column_size).toStdString();

		//then the value
		value.setNum(data[row].y(), 'E', OUTPUT_PRECISION);
		out_stream << value.toStdString() << endl;
	}

	//writing done, close stream
	out_stream.close();
}

/*
 * Private function that uses a dialog window to choose the filename/path to save to.
 */
QString OutputWriter::getFilePath(OutputType type) const {

	//The filter used will depend on the output type
	QString filter, default_name;

	switch(type) {

	case NUMBER_DENSITY:

		filter = "SIMtoEXP: *.sim, *.SIM (*.sim, *.SIM);; All files (*)";
		default_name = DEFAULT_FILENAME + NUMBER_DENSITY_FILE_SUFFIX;
		break;

	case ELECTRON_DENSITY:

		filter = "Electron density profile: *.edp, *.EDP (*.edp, *.EDP);; All files (*)";
		default_name = DEFAULT_FILENAME + ELECTRON_DENSITY_FILE_SUFFIX;
		break;

	case NEUTRON_DENSITY:

		filter = "Neutron Scattering Length density: *.nsld, *.NSLD (*.nsld, *.NSLD);; All files (*)";
		default_name = DEFAULT_FILENAME + NEUTRON_DENSITY_FILE_SUFFIX;
		break;

	case XRAY_STRUCTURE_FACTOR:

		filter = "X-ray structure factor: *.xsf, *.XSF (*.xsf, *.XSF);; All files (*)";
		default_name = DEFAULT_FILENAME + XRAY_STRUCTURE_FACTOR_SUFFIX;
		break;

	case NEUTRON_STRUCTURE_FACTOR:

		filter = "Neutron structure factor: *.nsf, *.NSF (*.nsf, *.NSF);; All files (*)";
		default_name = DEFAULT_FILENAME + NEUTRON_STRUCTURE_FACTOR_SUFFIX;
		break;

	case VOLUMETRIC_FIT:

		filter = "Component volumetric fit: *.cvf, *.CVF (*.cvf, *.CVF);; All files (*)";
		default_name = DEFAULT_FILENAME + VOLUMETRIC_FIT_FILE_SUFFIX;
		break;
	}

	//has built-in error/overwrite checking
	return QFileDialog::getSaveFileName(parent, "Select a path for saving:", default_name, filter);
}

/*
 * Does the work of writing all component output depending on the output type chosen
 */
void OutputWriter::writeComponents(OutputType type, const Components& components) {

	if(type == NUMBER_DENSITY) {

		QString file_path = getFilePath(NUMBER_DENSITY);

		//if the returned string is good, open the output stream
		if(not (file_path.isNull() or file_path.isEmpty())) {

			//append the suffix if it isn't already there
			if(not file_path.endsWith(NUMBER_DENSITY_FILE_SUFFIX, Qt::CaseInsensitive))
				file_path.append(NUMBER_DENSITY_FILE_SUFFIX);

			out_stream.open(file_path.toStdString());

			//first output header and commented area
			out_stream << header.toStdString() << endl;
			out_stream << "# Component Number Density" << endl;
			out_stream << "#***************************************************************************************" << endl;

			//get all of the densities and put them in a container
			QVector<Points> data;
			for(const auto& component : components) {

				data.append(component.getNumberDensities());
			}

			writeComponentData(data, components);
		}

	} else if(type == ELECTRON_DENSITY) {

		QString file_path = getFilePath(ELECTRON_DENSITY);

		//if the returned string is good, open the output stream
		if(not (file_path.isNull() or file_path.isEmpty())) {

			//append the suffix if it isn't already there
			if(not file_path.endsWith(ELECTRON_DENSITY_FILE_SUFFIX, Qt::CaseInsensitive))
				file_path.append(ELECTRON_DENSITY_FILE_SUFFIX);

			out_stream.open(file_path.toStdString());

			//first output header and commented area
			out_stream << header.toStdString() << endl;
			out_stream << "# Component Electron Density Profile (EDP)" << endl;
			out_stream << "#***************************************************************************************" << endl;

			//get all of the densities and put them in a container
			QVector<Points> data;
			for(const auto& component : components) {

				data.append(component.getElectronDensities());
			}

			writeComponentData(data, components);
		}

	} else if(type == NEUTRON_DENSITY) {

		QString file_path = getFilePath(NEUTRON_DENSITY);

		//if the returned string is good, open the output stream
		if(not (file_path.isNull() or file_path.isEmpty())) {

			//append the suffix if it isn't already there
			if(not file_path.endsWith(NEUTRON_DENSITY_FILE_SUFFIX, Qt::CaseInsensitive))
				file_path.append(NEUTRON_DENSITY_FILE_SUFFIX);

			out_stream.open(file_path.toStdString());

			//first output header and commented area
			out_stream << header.toStdString() << endl;
			out_stream << "# Component Neutron Scattering Length Profile (NSLP)" << endl;
			out_stream << "#***************************************************************************************" << endl;

			//get all of the densities and put them in a container
			QVector<Points> data;
			for(const auto& component : components) {

				data.append(component.getNeutronSLDensities());
			}

			writeComponentData(data, components);
		}
	}
}

/*
 *
 */
QString OutputWriter::addSpaces(const QString& input, uint total_length) const {

	QString padded = input;

	//add necessary spaces to even out columns
	while(padded.length() < total_length) {

		padded.append(' ');
	}

	//add two more spaces for delimiter
	padded.append("  ");

	return padded;
}

/*
 * Private function for defining the header that goes at the top of all text files
 */
void OutputWriter::setHeader() {

	header = "\
#***************************************************************************************\n\
# SIMtoEXP: software for the density analysis of simulated profiles\n\
#\n\
# v2.0 has switched to using Qt as a GUI library over Tcl/Tk/BLT\n\
#\n\
# original code and analysis routines by Norbert Kucerka at NRC Canada\n\
# 'Qt' code, trajectory reader and MARTINI added by:\n\
#\n\
# Bryan W. Holland\n\
# Centre for Molecular Simulations\n\
# University of Calgary\n\
# Calgary, Alberta, Canada\n\
#\n\
#***************************************************************************************\n\
#\n\
# If you use this program, you MUST cite:\n\
#\n\
# N. Kucerka, J. Katsaras and J.F. Nagle (2010) Comparing Membrane\n\
# Simulations to Scattering Experiments: Introducing the SIMtoEXP\n\
# Software. J. Membr. Biol., 235(1): 43-50\n\
#***************************************************************************************\n\
#                               ****\n\
#                   Norbert.Kucerka@nrc-cnrc.gc.ca\n\
#                     bryan.holland@ucalgary.ca\n\
#                               ****\n\
#***************************************************************************************";
}
