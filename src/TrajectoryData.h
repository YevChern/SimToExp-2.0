/*
 * TrajectoryData.h
 *
 *  Created on: Feb 25, 2013
 *      Author: bholland
 *
 * Reads and stores all the info for a particular system, and one frame from the
 * trajectory at a time in a System object. This is then used to update the appropriate
 * Grid(s) object. Can read through trajectories of the following formats:
 *
 * XTC (with a 'gro' file)
 * DCD (with a 'pdb' file)
 *
 * This is a QObject in order to connect to trajectory reading classes
 * that send out the coordinates in signals; this is done in order to update the System
 * object directly without having to give the System object to the readers, thus
 * removing a stronger dependency.
 *
 */

#ifndef GROMACSDATA_H_
#define GROMACSDATA_H_

#include "System.h"
#include "FileCode.h"
#include "AtomicInfo.h"
#include "NormalGrid.h"
#include "PolarGrid.h"
#include "DensityFrame.h"
#include "Axis.h"
#include "TrajectoryInfo.h"
#include "XtcFile.h"

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QProgressDialog>
#else
    #include <QtGui/QProgressDialog>
#endif

typedef QMap<ParticleLabel,AtomType> AtomTypeMap;
typedef QList<DensityFrame> DensityFrames;

class TrajectoryData : public QObject {

	Q_OBJECT;

	bool gro_file_loaded;

	System system_;
	SimGranularity granularity;
	AtomicInfo* info;
	AtomTypeMap elem_labels; //want them ordered, so need to use QMap
	AtomTypeMap particle_labels;

	//vector that defines the normal to the plane (and therefore also defines the plane) - user input
	Vector3D norm;
	Axis norm_axis;

	//Amount of D2O used in experimental system
	double heavy_water_fraction;

	//Trajectory details
	uint num_particles;
	uint num_lipids;
	uint first_step;
	double start_time;
	double time_step;
	uint num_frames;

	//these can/will grow over an extended trajectory, but will never shrink
	//since the system size shouldn't change too much, this should not be an issue.
	NormalGrid norm_grid; //grid normal (i.e. perpendicular) to the bilayer
	PolarGrid polar_grid; //grid parallel to plane of bilayer

	//all the calculated density data kept by element/residue
	DensityFrame elemental_frame;

	//calculated density data kept by particle/residue - needed for component models
	DensityFrame particle_frame;

	MDBox avg_box; //for ensemble simulations with fluctuating volume
	double bilayer_thickness;

public:

	TrajectoryData(QObject* parent = 0);
	TrajectoryData(AtomicInfo*,QObject* parent = 0);
	virtual ~TrajectoryData();

	//Gromacs formats
	FileCode readGroFile(const QString&);
	FileCode readXtcFile(const QString&,QProgressDialog&);

	//Charmm/NAMD formats (TODO: not implemented yet)
	FileCode readPdbFile(const QString&);
	FileCode readDdcFile(const QString&);

	void setNormalVector(Axis);
	void setNormalGridSize(double);
	void setPolarGridSizes(double,double);

	void setSimGranularity(SimGranularity);
	void setHeavyWaterFraction(double);

	double averageLipidArea() const;
	uint numLipids() const;

	void updateGrids();
	void updateFrames(uint,double);

	DensityFrames ensembleDensities();
	DensityFrames densitiesOverTime(uint);

	void clearFrames();

private:

	void calculateAvgBox(XtcFile&, QProgressDialog&, double);
	FileCode groFileError(const QString&);

};

#endif /* GROCOORDFILE_H_ */
