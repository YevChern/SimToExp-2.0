/*
 * PolarGrid.h
 *
 *  Created on: Mar 10, 2013
 *      Author: bholland
 */

#ifndef POLARGRID_H_
#define POLARGRID_H_

#include "ParticleIntensities.h"

#include <QtCore/QPointF>

class PolarGrid {

	double r_max;
	double r_min;
	double dr;

	double th_max;
	double th_min;
	double dth;

	QMap<AtomType,uint> total_intensities;

	/*
	 * 2D container for the grid, (r,theta), r = outer key, theta = inner key
	 * The value is the probability density for all of the atom types
	 */
	QMap<QPointF,ParticleIntensities> pol_grid;

public:

	PolarGrid();
	virtual ~PolarGrid();

	int numberOfRadii();
	int numberOfAngles();
	QMap<double,double>& values(double);

	void setRadiusRange(double,double);
	void setAngleRange(double,double);
	void createGrid(double,double);
};

#endif /* POLARGRID_H_ */
