/*
 * NormalGrid.cpp
 *
 *  Created on: Mar 10, 2013
 *      Author: bholland
 */

#include <iostream>
#include <limits>

#include <QtDebug>

#include "NormalGrid.h"
#include "MathUtils.h"
#include <cmath>

using namespace std;

//Default constructor
NormalGrid::NormalGrid() {

	dx = 0;
	info = NULL;

	//also set a single bin at x = 0
	norm_grid.insert(0, ParticleIntensities());

	norm_grid_itr = norm_grid.end();
}

//Constructor that sets the AtomicInfo pointer
NormalGrid::NormalGrid(AtomicInfo* info) {

	dx = 0;
	this->info = info;

	//also set a single bin at x = 0
	norm_grid.insert(0, ParticleIntensities());
	norm_grid_itr = norm_grid.end();
}

//Destructor
NormalGrid::~NormalGrid() {}

//setter for AtomicInfo pointer
void NormalGrid::setAtomicInfo(AtomicInfo* info) {

	this->info = info;
}

//Return the bin size
double NormalGrid::binSize() const {return dx;}

//setter for bin size
void NormalGrid::setBinSize(double dx) {

	this->dx = dx;
}

/*
 * Returns the number of bins currently in the grid
 */
uint NormalGrid::numBins() const {

	return norm_grid.size();
}

//Return the maximum and minimum real space points of the grid - these are the start of the bins though, so grid goes to max() + dx
double NormalGrid::max() const {

	//if the grid is empty, the max must be zero
	if(norm_grid.empty())
		return 0;

	return (--(norm_grid.end())).key();
}

double NormalGrid::min() const {

	//same as max, the min for an empty grid is just zero
	if(norm_grid.empty())
		return 0;

	return norm_grid.begin().key();
}

//Returns the histogram of the probability density (no normalization here, this is done later by volume) for the given particle
//The return value is a copy, but might want to consider a pointer as these object will be reasonably large, although
//memory management later becomes a little tricky
NormalDensity NormalGrid::normalDensity(const ParticleLabel& label, AtomType type) {

	//create and define a ParticleDensity object
	NormalDensity density(label);
	double elecs = info->numberElectrons(type);
	double neuts = info->neutronSL(type);
	density.setAtomInfo(type, elecs, neuts);

	//iterate through grid and create a vector of points of the x-coord and the prob density
	GridIntensityMap::iterator itr, itr_end = norm_grid.end();

	for(itr = norm_grid.begin(); itr != itr_end; ++itr) {

		QPointF point = probDensity(label, itr);
		density.addDensity(point);
	}

	return density;
}

//Increments the intensity of the particle in both containers
void NormalGrid::incrementIntensity(const ParticleLabel& label, double point, uint to_incr) {

	//check to see if the point is represented by the current bins in the map; if not,
	//need to create new bins up to the point
	GridIntensityMap::iterator bin;

	//if point is between min and max+dx, the bin should exist, otherwise create them first
	double lower_min = min();
	double upper_max = max() + dx;
	if(point < lower_min or point > upper_max) {

		//if lower, create bins downward until passed 'point'
		while(point < lower_min) {

			lower_min -= dx;
			norm_grid.insert(lower_min, ParticleIntensities());
		}

		//if higher, create bins upward until passed 'point'
		while(point > upper_max) {

			norm_grid.insert(upper_max, ParticleIntensities());
			upper_max += dx;
		}
	}

	//since the bin should exist, output an error message and return if not within the new bounds
	if(point < lower_min or point > upper_max) {

		cerr << "Bins not created that contain the point value: " << point << endl;
		return;
	}

	//now there should be a bin for the point, find it and increment
	//the bin value will be LARGER, so really in the previous bin
	bin = --(norm_grid.upperBound(point));
	bin.value().incrementIntensity(label, to_incr);
}

/*
 * Returns the intensity for the given ParticleLabel and closest to the provided grid value.  If label does
 * not exist, returns zero so be careful!
 */
uint NormalGrid::intensity(const ParticleLabel& label, double grid_value) const {

	const ParticleIntensities* intensities;

	//get the ParticleIntensities object closest to 'grid_value'
	GridIntensityMap::const_iterator low_bound_itr = norm_grid.lowerBound(grid_value), lower_bound_itr;

	//If the iterator is at the end, this is an error
	if(low_bound_itr == norm_grid.end()) {

		cerr << "Error - NormalGrid::intensity: Grid is empty or have iterated through!" << endl;
		return 0;
	}

	//if key and grid_value not exactly equal, need to check one grid value down (if it exists) to see
	//if it is closer
	if(low_bound_itr.key() != grid_value) {

		if(low_bound_itr != norm_grid.begin()) {

			lower_bound_itr = --low_bound_itr;

			//check which is the closest and set the intensities appropriately
			double low_diff = fabs(low_bound_itr.key() - grid_value);
			double lower_diff = fabs(lower_bound_itr.key() - grid_value);

			if(low_diff < lower_diff) {

				intensities = &low_bound_itr.value();

			} else {

				intensities = &lower_bound_itr.value();
			}


		} else {

			//has to be closest
			intensities = &low_bound_itr.value();
		}

	} else {

		intensities = &low_bound_itr.value();
	}

	//have the intensities, just return the intensity for the given ParticleLabel
	return intensities->intensity(label);
}

//Only resets the ParticleIntensities for the next frame, the grid points are left alone
void NormalGrid::resetIntensities() {

	//Iterate through grid clear all ParticleIntensity objects
	GridIntensityMap::iterator g_itr, g_itr_end = norm_grid.end();
	for(g_itr = norm_grid.begin(); g_itr != g_itr_end; ++g_itr)
		g_itr.value().clear();
}

//Clears the containers for the grid
void NormalGrid::clearGrid() {

	norm_grid.clear();
}

//Private function that returns the intensity for a specific particle and grid point
QPointF NormalGrid::probDensity(const ParticleLabel& label, const GridIntensityMap::iterator& itr) const {

	double intensity = itr.value().intensity(label);

	return QPointF(itr.key(),intensity);
}

/*
 * ----------------------------------------------------------------------------------
 * Iterator methods, primary intention is just for testing, shouldn't need to use them (but
 * whatever if you do I guess)
 * ------------------------------------------------------------------------------------
 */
const ParticleIntensities& NormalGrid::firstParticleInt() {

	norm_grid_itr = norm_grid.begin();

	return norm_grid_itr++.value(); //increment after in order to check for end
}

const ParticleIntensities& NormalGrid::nextParticleInt() {

	return norm_grid_itr++.value(); //increment after in order to check for end
}

bool NormalGrid::atEnd() const {

	if(norm_grid_itr == norm_grid.end())
		return true;
	return false;
}
