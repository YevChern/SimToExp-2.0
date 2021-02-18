/*
 * GroCoordFile.cpp
 *
 *  Created on: Feb 25, 2013
 *      Author: bholland
 */

#include "TrajectoryData.h"

#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QMessageBox>
#endif
#include <QtCore/QtDebug>

#include <iostream>
#include <cmath>

using namespace std;

//Default constructor
TrajectoryData::TrajectoryData(QObject* parent) : QObject(parent) {

	info = NULL;
	gro_file_loaded = false;
	norm = Vector3D();

	num_particles = 0;
	num_lipids = 0;
	first_step = 0;
	start_time = 0;
	time_step = 0;
	heavy_water_fraction = 0;
	num_frames = 0;

	norm_axis = Z_AXIS; // default to z, but should be user input
	granularity = SimGranularity::UNKNOWN_GRANULARITY;

	bilayer_thickness = 0;
}

//constructor that sets the AtomicInfo pointer
TrajectoryData::TrajectoryData(AtomicInfo* info, QObject* parent) : QObject(parent) {

	//set info pointer for all objects that need it
	this->info = info;
	norm_grid.setAtomicInfo(info);
	//TODO: also set AtomicInfo pointer for polar grid!!

	gro_file_loaded = false;
	norm = Vector3D();

	num_particles = 0;
	num_lipids = 0;
	first_step = 0;
	start_time = 0;
	time_step = 0;
	heavy_water_fraction = 0;
	num_frames = 0;

	norm_axis = Axis::Z_AXIS;
	granularity = SimGranularity::UNKNOWN_GRANULARITY;

	bilayer_thickness = 0;
}

//Destructor
TrajectoryData::~TrajectoryData() {}

/*
 * Return the calculated average total area of the lipids
 */
double TrajectoryData::averageLipidArea() const {

	double area = 0;

	//just the area along the NON-normal axes
	if(norm_axis == X_AXIS) {

		area = avg_box.y() * avg_box.z();
	}
	else if(norm_axis == Y_AXIS) {

		area = avg_box.x() * avg_box.z();
	}
	else if(norm_axis == Z_AXIS) {

		area = avg_box.x() * avg_box.y();
	}

	return area;
}

/*
 * Return the number of lipids in the system
 */
uint TrajectoryData::numLipids() const {

	return num_lipids;
}

/*
 * Set which axis corresponds to the normal vector. This is very specific: MUST be an axis (i.e. x, y or z).
 * Does not include an error check as the GUI will only allow
 * one of the three axes as options, and it is up to me not to screw up!
 */
void TrajectoryData::setNormalVector(Axis axis) {

	switch(axis) {

	case X_AXIS:

		norm.set(1, 0, 0);
		break;

	case Y_AXIS:

		norm.set(0, 1, 0);
		break;

	case Z_AXIS:

		norm.set(0, 0, 1);
		break;

	default:
		break;
	}

	//also set the axis in the object for reference
	this->norm_axis = axis;
}

//Set the grid size for the NormalGrid
void TrajectoryData::setNormalGridSize(double grid_size) {

	norm_grid.setBinSize(grid_size);
}

//Set the grid sizes for the PolarGrid
void TrajectoryData::setPolarGridSizes(double r_size, double th_size) {

	//TODO: finish this function once PolarGrid is ready
}

//Set the granularity of the trajectory as given by the user
void TrajectoryData::setSimGranularity(SimGranularity granularity) {

	this->granularity = granularity;
}


//Set the fraction of heavy water to be used for this particular run
void TrajectoryData::setHeavyWaterFraction(double fraction) {

	heavy_water_fraction = fraction;
}

//Reads in the 'gro' coordinate file, ignoring any velocity data.
FileCode TrajectoryData::readGroFile(const QString& path) {

	//open the file
	QFile file(path);

	//try opening file, if doesn't work, output an error box
	if(not file.open(QIODevice::ReadOnly | QIODevice::Text)) {

		QMessageBox::critical(0, "Unknown file error:", "Cannot open file :-(");
		return FILE_ERROR;
	}

	//Set file as a QTextStream object - easy to use and already using Qt everywhere
	QTextStream input(&file);

	//first line is just name of system
	QString next_line = input.readLine();
	system_.setName(next_line);

	//get number of particles to set for loop length
	bool nparts_ok;
	uint num_particles = (input.readLine()).toUInt(&nparts_ok);

	//if any problems with the data, show message and exit input
	if(not nparts_ok) {

		//show message, clean up and return an error
		QString msg = "Number of particles must be an unsigned integer.";
		return groFileError(msg);
	}

	//now go through rest of file and store data for each particle
	Residue* res = NULL;
	for(uint i = 0; i < num_particles; i++) {

		//get the next line
		QString line = input.readLine();

		bool resid_ok, pid_ok, x_ok, y_ok, z_ok; //for error check

		//File format for a 'gro' file can be found at: http://manual.gromacs.org/current/online/gro.html
		uint resID = (line.left(5)).toUInt(&resid_ok); //first 5 positions are for residue ID
		QString res_name = (line.mid(5, 5)); //next 5 for residue name
		QString part_name = (line.mid(10, 5)); //next 5 for atom name
		uint partID = (line.mid(15, 5)).toUInt(&pid_ok); //next 5 for atom number
		float x = (line.mid(20, 8)).toFloat(&x_ok); //next 8 for x coordinate
		float y = (line.mid(28, 8)).toFloat(&y_ok); //next 8 for y coordinate
		float z = (line.mid(36, 8)).toFloat(&z_ok); //next 8 for z coordinate

		//need to get rid of white space from strings
		res_name = res_name.trimmed();
		part_name = part_name.trimmed();

		//check if the residue is a lipid
		ResidueType lipid_type = info->residueType(res_name.toStdString());

		//for first particle only
		if(i == 0) {

			res = &(system_.addResidue(res_name, resID, lipid_type));

			//if a lipid, increment the lipid count
			if(info->isLipid(lipid_type))
				num_lipids++;
		}

		//if any problems with the data, show message and exit input
		if(not(resid_ok and pid_ok and x_ok and y_ok and z_ok)) {

			//show message, clean up and return an error
			QString msg = "One or more of the input values on line ";
			msg += (i+2);
			msg += " is invalid.";

			return groFileError(msg);
		}

		//Check for zero value in both residue ID and particle ID - must be > 0
		if(resID == 0 or partID == 0) {

			//show message, clean up and return an error
			QString msg = "Both the residue IDs and the particle IDs must be >= 1.";
			return groFileError(msg);
		}

		//everything is ok, enter all data into residue.
		if(res->resId() != resID) {

			//if a lipid, increment the lipid count
			if(info->isLipid(lipid_type))
				num_lipids++;

			//if not the same as previous residue, need to add previous residue to system and create a new one
			res = &(system_.addResidue(res_name, resID, lipid_type));
		}

		Particle& particle = res->addParticle(part_name, partID);
		particle.setCoord(x, y, z);

		//also add the particle to the system pointer list
		system_.addParticle(&particle);

		//get the AtomType and set it
		string std_part_name = part_name.toStdString();
		string std_res_name = res_name.toStdString();

		//if the system is atomistic, try to find the atom type, if it fails strip the last
		//letter and try again until something works
		AtomType type;
		if(granularity == ATOMISTIC) {

			type = info->atomType(std_part_name);
			while(type == AtomType::UNKNOWN){

				//if all characters removed, not recognized by AtomicInfo
				if(std_part_name.length() == 0) {

					//return an error
					QString msg = "The particle ";
					msg += part_name;
					msg += " in residue ";
					msg += res_name;
					msg += " is unknown to SIMtoEXP, cannot continue.";

					return groFileError(msg);
				}

				std_part_name = std_part_name.substr(0,std_part_name.length()-1);
				type = info->atomType(std_part_name);
			}
		}
		else {

			type = info->atomType(std_res_name, std_part_name);
		}

		//If 'UNKNOWN', this is a serious problem, exit with file code error
		if(type == AtomType::UNKNOWN) {

			//return an error
			QString msg = "The particle ";
			msg += part_name;
			msg += " in residue ";
			msg += res_name;
			msg += " is unknown to SIMtoEXP, cannot continue.";

			return groFileError(msg);
		}

		//otherwise, set the particle's "recognized" name to the string passed to atomType
		particle.setRecognizedName(QString(std_part_name.c_str()));
		particle.setAtomType(type);

		//also want to adjust the particle coordinate file name to start with the recognized portion
		QString simtoexp_name = particle.recognizedName();

        if (simtoexp_name.size() < part_name.size()){
        simtoexp_name += "-";
        simtoexp_name += part_name.right(part_name.size() - particle.recognizedName().size());
        particle.setCoordFileName(simtoexp_name);
        }
		//Finally, want to keep a set of all unique ParticleLabels
		QPair<ParticleLabel,AtomType> elemental_label = res->elementalLabel(particle.id());

		if(not elem_labels.contains(elemental_label.first))
			elem_labels.insert(elemental_label.first, elemental_label.second);

		QPair<ParticleLabel,AtomType> particle_label = res->particleLabel(particle.id());

		if(not particle_labels.contains(particle_label.first))
			particle_labels.insert(particle_label.first, particle_label.second);
	}

	//should be at last line, read in the first three values, ignore the rest, don't deal with any geometry that isn't a box!
	QString value;
	bool x_ok, y_ok, z_ok;
	float x_dim, y_dim, z_dim;

	input >> value;
	x_dim = value.toFloat(&x_ok);
	input >> value;
	y_dim = value.toFloat(&y_ok);
	input >> value;
	z_dim = value.toFloat(&z_ok);

	//check for problems
	if(not(x_ok and y_ok and z_ok)) {

		QString msg = "One or more of the box dimensions is invalid.";
		return groFileError(msg);
	}

	system_.setBoxDimensions(x_dim, y_dim, z_dim);

	//clean up
	file.close();

	//otherwise we're done :)
	gro_file_loaded = true;
	return READ_SUCCESS;
}

//Reads in the given XTC file and calculates the particle probability densities.  Does this for
//a small number of frames at a time in order to use parallel processing.
//TODO: use threads to calculate the ParticleDensity objects while - each thread
//should get its own XtcTrajectory object, then delete the XtcFrames when it is done.
//TODO: do some optimization for number of frames read at a time and number of threads
FileCode TrajectoryData::readXtcFile(const QString& path, QProgressDialog& progress) {

	//first check to see if 'gro' file has been loaded - absolutely must be!
	if(not gro_file_loaded)
		return groFileError("Gromacs coordinate file ('gro') must be loaded before XTC can be read!");

	//TODO: for now, just create a single XtcTrajectory class and read in all of the data at once,
	//once that is up and working, move on to threaded version
	XtcFile xtc_file(path, &system_);

	//try opening file, if doesn't work open a msg box and return error
	if(not xtc_file.openXtc()) {

		QString msg = "Can't open 'xtc' file.";
		return groFileError(msg);
	}

	//Read the header info and pass along to trajectory
	if(not xtc_file.readHeader()) {

		xtc_file.close();
		QString msg = "Problem reading header file.";
		return groFileError(msg);
	}

	//set the progress bar max (that unfortunately only takes qint64)
	qint64 file_size_64 = xtc_file.size();

	//determine the scaling factor by the file size
	int scaling_factor = 1; //default for all files smaller than INT_MAX
	if(file_size_64 > INT_MAX) {

		scaling_factor = (file_size_64 - INT_MAX) / INT_MAX + 10; //add ten for good measure so there are no boundary issues, can easily scale by ten anyway
	}

	int scaled_file_size = file_size_64/scaling_factor;
	progress.setMaximum(scaled_file_size);

	//once file is open and header has been read, do a first loop through the file to calculate the average MDBox
    /*
     * NOTE: no need to do it prior to reading all the data any more. We can do this during main cycle through XTC file.
     * Also, in this implementation of density calculations, we will read XTC file only once instead of twice.
    */
    //calculateAvgBox(xtc_file, progress, scaling_factor);

    Matrix3D cumulative_box;  // container for average MDBox calculation (for average area per lipid)

	//now reading data
	progress.setLabelText("Reading densities from Gromacs XTC file...");

	num_particles = xtc_file.numberOfAtoms();
	first_step = xtc_file.firstStep();
	start_time = xtc_file.startTime();

	//read the first two frames here in order to calculate the time step
	//should always be at least two frames in a trajectory, the start and end!
	//Also, each time a frame is read, need to update the Grids from the recently updated System member
	xtc_file.readNextFrame(); //first frame
	updateGrids();
	updateFrames(xtc_file.currentStep(), xtc_file.currentTime());
    num_frames++;
    //add to the cumulative total
    cumulative_box += system_.boxDimensions();

	xtc_file.readNextFrame(); //second frame
	time_step = xtc_file.timeStep();
	updateGrids();
	updateFrames(xtc_file.currentStep(), xtc_file.currentTime());
    num_frames++;
    //add to the cumulative total
    cumulative_box += system_.boxDimensions();

	//Now go through all the frames and collect the info in the
	while(not xtc_file.atEnd()) {

		int scaled_pos = xtc_file.pos()/scaling_factor;
		progress.setValue(scaled_pos);
		xtc_file.readNextFrame();
		updateGrids();
		updateFrames(xtc_file.currentStep(), xtc_file.currentTime());
        num_frames++;
        //add to the cumulative total
        cumulative_box += system_.boxDimensions();
	}

    //divide by total number of frames to get the average
    double inv_num_frames = 1.0 / num_frames;
    cumulative_box *= inv_num_frames;

    avg_box.setBoxDimensions(cumulative_box);

	progress.setValue(scaled_file_size);

	//clean up the file
	xtc_file.close();

	//made it all the way, return success
	return READ_SUCCESS;
}

//Assumes the System is up to date, and proceeds to do the work of populating both the normal and polar grids
//want to do both at the same time to avoid having to loop through the system twice.
//NOTE: assumes the grids sizes have already been set as well
void TrajectoryData::updateGrids() {

	//first clear the grids of any previous data
	norm_grid.resetIntensities();

    //compute lipid COM (in normal direction)
    double com = 0;
    double lip_atom_total_mass = 0;

    for(Residue residue = system_.firstResidue(); not system_.atEnd(); residue = system_.nextResidue()) {

        if (info->isLipid(info->residueType(residue.name().toStdString()))) {
            for(Particle particle = residue.firstParticle(); not residue.atEnd(); particle = residue.nextParticle()) {

                com += particle.coord().dot(norm)*(info->mass(particle.atomType()));
                lip_atom_total_mass += info->mass(particle.atomType());

            }
        }
    }

    com /= lip_atom_total_mass;

	//if empty of residues, will just skip the loop
	for(Residue residue = system_.firstResidue(); not system_.atEnd(); residue = system_.nextResidue()) {

		for(Particle particle = residue.firstParticle(); not residue.atEnd(); particle = residue.nextParticle()) {

			//TODO: for now only have Normal grid here - need to add polar grid once the class has been developed

            // compute distance from COM to current atom, taking into account pbc
            Vector3D d = system_.systemBox().realSpace().inverse()*(norm*com - particle.coord()); // convert vector (from center to particle) to box coordinates
            for (int i=0; i<3; ++i)
                d[i] -= round(d[i]); // subtract integral number of box vectors
            double position = (system_.systemBox().realSpace()*d).dot(norm); // convert back to cartesian coordinates and then leave only normal component

			ParticleLabel elem_label = residue.elementalLabel(particle.id()).first;
			ParticleLabel particle_label = residue.particleLabel(particle.id()).first;

			norm_grid.incrementIntensity(elem_label, position, 1); //increment by one
			norm_grid.incrementIntensity(particle_label, position, 1); //increment by one
		}
	}
}

//Assumes the grids are up to date and uses them to create and store a Frame object (that of course corresponds to the current frame
//from the simulation file)
void TrajectoryData::updateFrames(uint step, double time) {

	//Do elemental frame, then particle frame

	DensityFrame elem_frame(step, time, norm_axis);
	elem_frame.setSimBox(system_.systemBox());

	//go through all known particle types and obtain the ParticleDensity that corresponds to each
	AtomTypeMap::iterator itr, itr_end = elem_labels.end();
	for(itr = elem_labels.begin(); itr != itr_end; ++itr) {

		NormalDensity density = norm_grid.normalDensity(itr.key(), itr.value());

		//need to set the residue type prior to adding to the frame
		ResidueType res = info->residueType(density.residueName().toStdString());
		density.setResidueType(res);
		elem_frame.addNormalDensity(itr.key(), density);
	}

	//now have all of the data for the current frame

	//before scaling the axis, want to normalize with respect to slab volume
	elem_frame.normalizeNormalValues(norm_grid.binSize());

	//also need to deal with the variable box size
    //elem_frame.scaleNormalAxis(bilayer_thickness);

	//add the current density frame to the running total
	elemental_frame.addNormalDensityFrame(elem_frame);

	DensityFrame p_frame(step, time, norm_axis);
	p_frame.setSimBox(system_.systemBox());

	//go through all known particle types and obtain the ParticleDensity that corresponds to each
	//Can't use handy new C++11 ranged loop, need the iterator itself to get the key
	itr_end = particle_labels.end();
	for(itr = particle_labels.begin(); itr != itr_end; ++itr) {

		NormalDensity density = norm_grid.normalDensity(itr.key(), itr.value());

		//need to set the residue type prior to adding to the frame
		ResidueType res = info->residueType(density.residueName().toStdString());
		density.setResidueType(res);
		p_frame.addNormalDensity(itr.key(), density);
	}

	//now have all of the data for the current frame

	//before scaling the axis, want to normalize with respect to slab volume
	p_frame.normalizeNormalValues(norm_grid.binSize());

	//also need to deal with the variable box size
    //p_frame.scaleNormalAxis(bilayer_thickness);

	//add the current density to the running total
	particle_frame.addNormalDensityFrame(p_frame);
}

/*
 * Scales the total values for the Density Frames after all frames have been read and added.  Returns
 * a container of all different DensityFrame types calculated, details below.
 */
DensityFrames TrajectoryData::ensembleDensities() {



	//now have all DensityFrames added, so just need to divide by total frame number
	elemental_frame.scaleNormalValues(num_frames);
	particle_frame.scaleNormalValues(num_frames);

    //NOTE: not called now
    //
	//some adjustments: center the plots about zero
    //elemental_frame.centerNormalValues();
    //particle_frame.centerNormalValues();

	//one last thing, adjust the hydrogen from water to reflect the fraction of heavy water requested by the user
	elemental_frame.adjustDeuteriumFraction(heavy_water_fraction, info, granularity, true);
	particle_frame.adjustDeuteriumFraction(heavy_water_fraction, info, granularity, false);

	DensityFrames frames;

	//MUST BE ADDED IN ORDER: ELEMENTAL, PARTICLE
	frames.append(elemental_frame);
	frames.append(particle_frame);

	return frames;
}

//Recompiles the QList of DensityFrames by averaging the particles into 'clusters' of frames
//This is intended to remove fluctuations over small time changes so that changes in structure
//that occur over long time periods are easier to see
DensityFrames TrajectoryData::densitiesOverTime(uint cluster_size) {


}

//IMPORTANT: this MUST be called prior to a TrajectoryData object going out of scope, or else
//a horrible memory leak will develop - this will delete all the frames, but more importantly
//all the objects pointed to by the NormalDensity pointers hiding inside them
void TrajectoryData::clearFrames() {

	elemental_frame.clearNormalDensities();
	particle_frame.clearNormalDensities();
}

// NOTE: not called now
//
//private function that goes through the trajectory (assuming the file is open and header read) and calculates
//the average MDBox that will be used to scale the bin widths later on
void TrajectoryData::calculateAvgBox(XtcFile& xtc_file, QProgressDialog& progress, double scaling_factor) {

	//set the progress bar name to tell user avg box is being calculated
	progress.setLabelText("Calculating average box from Gromacs XTC file...");

	Matrix3D cumulative_box;

	//Now go through all the frames and collect the info in the
	while(not xtc_file.atEnd()) {

		int scaled_pos = xtc_file.pos()/scaling_factor;
		progress.setValue(scaled_pos);
		xtc_file.readNextFrame();
		num_frames++;

		//add to the cumulative total
		cumulative_box += system_.boxDimensions();
	}

	//divide by total number of frames to get the average
	double inv_num_frames = 1.0 / num_frames;
	cumulative_box *= inv_num_frames;

	avg_box.setBoxDimensions(cumulative_box);

	//use this to determine which axis to scale by
	switch(norm_axis) {

	case X_AXIS:
		bilayer_thickness = avg_box.x();
		break;

	case Y_AXIS:
		bilayer_thickness = avg_box.y();
		break;

	case Z_AXIS:
		bilayer_thickness = avg_box.z();
		break;

	default:
		cerr << "Error - TrajectoryData::ensembleDensities(): Normal axis not set!" << endl;
	}

	//also set the DensityFrame information now that we have the average box
	elemental_frame.setSimBox(avg_box);
	elemental_frame.setNormalAxis(norm_axis);

	particle_frame.setSimBox(avg_box);
	particle_frame.setNormalAxis(norm_axis);

	//reset xtc file and the progress bar
	xtc_file.resetToStart();
	progress.reset();
}

////private helper function for dealing with file errors in 'gro' files
FileCode TrajectoryData::groFileError(const QString& msg) {

	QMessageBox::critical(0, "Gro file error:", msg);
	system_.clearSystem();
	return FILE_ERROR;
}
