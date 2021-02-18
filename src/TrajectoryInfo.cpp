/*
 * TrajectoryInfo.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: bholland
 */

#include "TrajectoryInfo.h"

const double DEFAULT_NORMAL_GRID_SPACING = 0.2; //in Angstroms
const double DEFAULT_RADIAL_GRID_SPACING = 0.2; //in Angstroms
const double DEFAULT_ANGULAR_GRID_SPACING = 0.05; // in radians

//Default constructor - sets all values to defaults
TrajectoryInfo::TrajectoryInfo() {

	setDefaultValues();
}

//Destructor
TrajectoryInfo::~TrajectoryInfo() {}

//getters
const QString& TrajectoryInfo::filePath() const {return file_path;}
const QString& TrajectoryInfo::optFilePath() const {return opt_file_path;}
TrajectoryType TrajectoryInfo::trajectoryType() const {return type;}
SimGranularity TrajectoryInfo::simGranularity() const {return granularity;}
Axis TrajectoryInfo::normalAxis() const {return normal_axis;}

double TrajectoryInfo::heavyWaterFraction() const {return heavy_fraction;}

bool TrajectoryInfo::calculateNormalToSurface() const {return calc_normal;}
bool TrajectoryInfo::normalTimeAvg() const {return normal_time_avg;}
bool TrajectoryInfo::normalFunctionOfTime() const {return normal_time_funct;}
uint TrajectoryInfo::normalSamplesPerPoint() const {return normal_num_samples;}
double TrajectoryInfo::normalGridSize() const {return normal_grid_size;}

bool TrajectoryInfo::calculateParallelToSurface() const {return calc_planar;}
bool TrajectoryInfo::planarTimeAvg() const {return planar_time_avg;}
bool TrajectoryInfo::planarFunctionOfTime() const {return planar_time_funct;}
uint TrajectoryInfo::planarSamplesPerPoint() const {return planar_num_samples;}
double TrajectoryInfo::radialGridSize() const {return radial_grid_size;}
double TrajectoryInfo::angularGridSize() const {return angular_grid_size;}

//setters
void TrajectoryInfo::setFilePath(const QString& path) {file_path = path;}
void TrajectoryInfo::setOptFilePath(const QString& path) {opt_file_path = path;}
void TrajectoryInfo::setTrajectoryType(TrajectoryType type) {this->type = type;}
void TrajectoryInfo::setSimGranularity(SimGranularity gran) {granularity = gran;}
void TrajectoryInfo::setNormalAxis(Axis axis) {normal_axis = axis;}

void TrajectoryInfo::setHeavyWaterFraction(double fraction) {heavy_fraction = fraction;}

void TrajectoryInfo::setCalcNormalToSurface(bool value) {calc_normal = value;}
void TrajectoryInfo::setNormalTimeAvg(bool value) {normal_time_avg = value;}
void TrajectoryInfo::setNormalFunctionOfTime(bool value) {normal_time_funct = value;}
void TrajectoryInfo::setNormalSamplesPerPoint(uint num) {normal_num_samples = num;}
void TrajectoryInfo::setNormalGridSize(double size) {normal_grid_size = size;}

void TrajectoryInfo::setCalcParallelToSurface(bool value) {calc_planar = value;}
void TrajectoryInfo::setPlanarTimeAvg(bool value) {planar_time_avg = value;}
void TrajectoryInfo::setPlanarFunctionOfTime(bool value) {planar_time_funct = value;}
void TrajectoryInfo::setPlanarSamplesPerPoint(uint num) {planar_num_samples = num;}
void TrajectoryInfo::setRadialGridSize(double size) {radial_grid_size = size;}
void TrajectoryInfo::setAngularGridSize(double size) {angular_grid_size = size;}

//for setting all of the trajectory values to the default
void TrajectoryInfo::setDefaultValues() {

		type = GROMACS_XTC;
		normal_axis = Z_AXIS;
		granularity = ATOMISTIC;

		heavy_fraction = 0;

		calc_normal = true;
		normal_time_avg = true;
		normal_time_funct = false;
		normal_num_samples = 1;
		normal_grid_size = DEFAULT_NORMAL_GRID_SPACING;

		calc_planar = false;
		planar_time_avg = true;
		planar_time_funct = false;
		planar_num_samples = 1;
		radial_grid_size = DEFAULT_RADIAL_GRID_SPACING;
		angular_grid_size = DEFAULT_ANGULAR_GRID_SPACING;
}
