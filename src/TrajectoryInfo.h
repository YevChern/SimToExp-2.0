/*
 * TrajectoryInfo.h
 *
 *  Created on: Apr 22, 2013
 *      Author: bholland
 *
 *  Class that collects all the information necessary for reading trajectory files, including what
 *  about the trajectory files is of interest, e.g. should the probabilities be calculated parallel
 *  to the surface.
 */

#ifndef TRAJECTORYINFO_H_
#define TRAJECTORYINFO_H_

#include <QtCore/QString>

#include "Axis.h"
#include "SimGranularity.h"

//for choice of file to open
enum TrajectoryType {GROMACS_XTC, NAMD_DCD, AMBER_CRD};

typedef unsigned int uint;

class TrajectoryInfo {

	QString file_path;
	QString opt_file_path;
	TrajectoryType type;

	SimGranularity granularity;
	Axis normal_axis;

	double heavy_fraction; // fraction of heavy water

	bool calc_normal;
	bool normal_time_avg;
	bool normal_time_funct;
	uint normal_num_samples;
	double normal_grid_size;

	bool calc_planar;
	bool planar_time_avg;
	bool planar_time_funct;
	uint planar_num_samples;
	double radial_grid_size;
	double angular_grid_size;

public:

	TrajectoryInfo();
	virtual ~TrajectoryInfo();

	//getters
	const QString& filePath() const;
	const QString& optFilePath() const;
	TrajectoryType trajectoryType() const;
	SimGranularity simGranularity() const;
	Axis normalAxis() const;

	double heavyWaterFraction() const;

	bool calculateNormalToSurface() const;
	bool normalTimeAvg() const;
	bool normalFunctionOfTime() const;
	uint normalSamplesPerPoint() const;
	double normalGridSize() const;

	bool calculateParallelToSurface() const;
	bool planarTimeAvg() const;
	bool planarFunctionOfTime() const;
	uint planarSamplesPerPoint() const;
	double radialGridSize() const;
	double angularGridSize() const;

	//setters
	void setFilePath(const QString&);
	void setOptFilePath(const QString&);
	void setTrajectoryType(TrajectoryType);
	void setSimGranularity(SimGranularity);
	void setNormalAxis(Axis);

	void setHeavyWaterFraction(double);

	void setCalcNormalToSurface(bool);
	void setNormalTimeAvg(bool);
	void setNormalFunctionOfTime(bool);
	void setNormalSamplesPerPoint(uint);
	void setNormalGridSize(double);

	void setCalcParallelToSurface(bool);
	void setPlanarTimeAvg(bool);
	void setPlanarFunctionOfTime(bool);
	void setPlanarSamplesPerPoint(uint);
	void setRadialGridSize(double);
	void setAngularGridSize(double);

	void setDefaultValues();
};

#endif /* TRAJECTORYINFO_H_ */
