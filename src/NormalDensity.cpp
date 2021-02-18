/*
 * NormalDensity.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: bholland
 */  


#include "NormalDensity.h"
#include "AtomicInfo.h"

#include <iostream>
#include <cfloat>

#include <QtCore/QPair>
#include <QtDebug>

typedef QPair<double,double> BinBoundaries;

const double NONZERO_TOLERANCE = 1E-4;

using namespace std;

//default constructor
NormalDensity::NormalDensity() {

   num_electrons = 0;
   neutron_sl = 0;
   atom_type = UNKNOWN;
   res_type = UNKNOWN_RESIDUE;
   all_densities_updated = true;
   densities_itr = densities_.end(); //initialize at end
   bin_size = 0;
   bin_size_updated = false;
}

//Constructor that sets the name
NormalDensity::NormalDensity(const ParticleLabel& p_label) {

   label = p_label;
   num_electrons = 0;
   neutron_sl = 0;
   atom_type = UNKNOWN;
   res_type = UNKNOWN_RESIDUE;
   all_densities_updated = true;
   densities_itr = densities_.end(); //initialize at end
   bin_size = 0;
   bin_size_updated = false;
}

//destructor
NormalDensity::~NormalDensity() {
   
}

//return the particle name (id)
QString NormalDensity::particleName() {return label.particleName();}
const QString& NormalDensity::particleName() const {return label.particleName();}

//return the residue name
QString NormalDensity::residueName() {return label.residueName();}
const QString& NormalDensity::residueName() const {return label.residueName();}

//return the particle label itself
ParticleLabel NormalDensity::particleLabel() {return label;}
const ParticleLabel& NormalDensity::particleLabel() const {return label;}

//A couple of setters for the particle label info (i.e. for changing them independently)
void NormalDensity::setParticleName(const QString& name) {

	label.setParticleName(name);
}

void NormalDensity::setResidueName(const QString& name) {

	label.setResidueName(name);
}

//set the name
void NormalDensity::setLabel(const ParticleLabel& label) {this->label = label;}

//return the element
AtomType NormalDensity::atomType() const {return atom_type;}
ResidueType NormalDensity::residueType() const {return res_type;}

//Returns true if the residue type is a lipid - just uses enum to determine, so all lipids
//must be contiguous in the enum
bool NormalDensity::residueIsLipid() const {

	if(res_type < DHPC or res_type > BOLB)
		return false;
	return true;
}

//set the lipid type
void NormalDensity::setResidueType(ResidueType res) {res_type = res;}

//getters
double NormalDensity::numberOfElectrons() const {return num_electrons;}
double NormalDensity::neutronSL() const {return neutron_sl;}
double NormalDensity::binSize() const {return bin_size;}

//Return a reference to the container of densities (i.e. points); updates the containers if necessary
Points& NormalDensity::densities() {

	return vec_densities;
}

/*
 * Return a copy of the Points in the map, avoids calculating all densities or storing
 * any extra info if not required (i.e. for NormalDensity objects in a DensityFrame).
 */
Points NormalDensity::densitiesFromMap() {

	Points densities;

	//go through vector for densities and create points by multiplying by
	//electron number and neutron SL
	DensitiesMap::iterator itr, itr_end = densities_.end();

	for(itr = densities_.begin(); itr != itr_end; ++itr) {

		double cur_x = itr.key();
		double cur_d = itr.value();

		//standard density container
		densities.append(QPointF(cur_x, cur_d));
	}

	return densities;
}

//return the first and last points (respectively) entered into the vector.  Also
//sets the iterator to the beginning
QPointF NormalDensity::firstPoint() {

	QPointF first;
	densities_itr = densities_.begin();
	first.setX(densities_itr.key());
	first.setY(densities_itr.value());

	//move the iterator ahead before returning - only way to check for 'end'
	if(densities_itr != densities_.end())
		++densities_itr;

	return first;
}

//Return the first point that has a non-zero value - does not affect the iterator
QPointF NormalDensity::firstNonZeroPoint() {

	QPointF first;

	//iterate until reach a non-zero point
	DensitiesMap::iterator itr = densities_.begin(), itr_end = densities_.end();
	while(itr.value() < NONZERO_TOLERANCE and itr != itr_end)
		++itr;

	//if made it to end, probably an error since there should be data to use this method
	if(itr == itr_end)
		cerr << "NormalDensity::getFirstNonZeroPoint: all values non-zero, no data in object?";

	//otherwise we have our point
	first.setX(itr.key());
	first.setY(itr.value());

	return first;
}

//Gets the point currently pointed to by the iterator.  If reaches 'end',
//return origin with an error message to stderr, should check with 'densitiesAtEnd' first.
QPointF NormalDensity::nextPoint() {

	QPointF next;

	if(densities_itr != densities_.end()) {

		next.setX(densities_itr.key());
		next.setY(densities_itr.value());

		//move the iterator ahead before returning, needed for checking for the end
		++densities_itr;

	} else {

		//at end, return origin
		cerr << "Density iterator at end of container, returning origin" << endl;
		next.setX(0);
		next.setY(0);
	}

	return next;
}

QPointF NormalDensity::lastPoint() {

	QPointF last;

	//check to see if there is even a point
	if(densities_.begin() != densities_.end()) {

		DensitiesMap::iterator itr = densities_.end();
		last.setX((--itr).key()); //move back from end
		last.setY(itr.value());

	} else {

		//otherwise just return origin with error message output
		cerr << "No points in density function, returning origin" << endl;
		last.setX(0);
		last.setY(0);
	}
	return last;
}

//Return the last point that has a non-zero value - does not affect the iterator
QPointF NormalDensity::lastNonZeroPoint() {

	QPointF last;

	//iterate until reach a non-zero point
	DensitiesMap::iterator itr = --densities_.end(), itr_beg = densities_.begin();
	while(itr.value() < NONZERO_TOLERANCE and itr != itr_beg)
		--itr;

	//if made it to end, probably an error since there should be data to use this method
	if(itr == itr_beg)
		cerr << "NormalDensity::getLastNonZeroPoint: all values non-zero, no data in object?";

	//otherwise we have our point
	last.setX(itr.key());
	last.setY(itr.value());

	return last;
}


//Checks to see if the iterator is at the end of the container
bool NormalDensity::densitiesAtEnd() const {

	if(densities_itr == densities_.end())
		return true;
	return false;
}

//Boolean to see if density map is empty
bool NormalDensity::isEmpty() const {

	return densities_.empty();
}

/*
 * Equality operators: equal iff ParticleLabel, AtomType, ResidueType are the same. Not concerned with
 * distributions details.
 */
bool NormalDensity::operator== (const NormalDensity& other) const {

	if(other.particleName() == label.particleName() and
	   other.residueName() == label.residueName() and
	   other.atomType() == atom_type and
	   other.residueType() == res_type) {

		return true;
	}

	return false;
}

bool NormalDensity::operator!= (const NormalDensity& other) const {

	return not (*this == other);
}

//Boolean that returns true if all density containers are synched
bool NormalDensity::densitiesCalculated() const {return all_densities_updated;}

//return the symmetrized version of the number densities; ensures container is up to date
Points& NormalDensity::symmetricDensities() {

	return sym_densities;
}

//return the container for electron densities
Points& NormalDensity::electronDensities() {

	return elec_densities;
}

//return the container for neutron densities
Points& NormalDensity::neutronDensities() {

	return neut_densities;
}

/*
 * Adds a new density point to the container
 * CAUTION: Currently does not check to see if 'z' is already in the Map
 */
void NormalDensity::addDensity(double z, double value) {
 
   densities_.insert(z, value);
   all_densities_updated = false;

   //check to see if can set the bin size
   if(not bin_size_updated and densities_.size() > 1)
	   setBinSize();
}

/*
 * Overloaded: same as above except takes in QPointF
 */
void NormalDensity::addDensity(QPointF point) {

   addDensity(point.x(), point.y());
}

//Add the results from a separate NormalDensity object. The values 'z' values will essentially
//never be equal unless a constant volume is used during the simulation, since the densities will
//be scaled to account for the change in volume (e.g. NPT ensemble, most common)
void NormalDensity::addDensities(NormalDensity other) {

	//go through 'other' densities
	Points other_points = other.densitiesFromMap();

	//if the current NormalDensity object is empty, simply add the 'other' densities and return
	if(densities_.isEmpty()) {

		for(auto point : other_points)
			addDensity(point);

		return;
	}


	//need access to a counter, so can't use ranged for loop
	int other_size = other_points.size();
    int this_size = densities_.size();

    //add densities
    for (int i=0; i<other_size; ++i){
        densities_[other_points[i].x()] += other_points[i].y();
    }

    //check if densities_.size is correct after addition
    //(we use map for densities_, so can create extra point(s) unintentionally if keys for x are not exactly the same for this and other densities)
    if ((densities_.size() != other_size) && (densities_.size() != this_size)){
        cerr << "Wrong size of NormalDensity after addition!" << endl;
    }


    /* -------------------------------------------------------------------
     * Old version with variable bin size
     *
    //

	for(int i = 0; i < other_size; ++i) {

		//need containers to keep track of all the bins to add the other bin to
		QVector<double> density_keys;
		QVector<BinBoundaries> bounds;

		const QPointF& point = other_points[i];

		//first thing needed is the bin boundaries of the current point being added
		double other_lowerbin = point.x();
		double other_upperbin = other_lowerbin + other.binSize();

		//use the logarithmic searching of the map container to find the bounds on the current point
		DensitiesMap::iterator f_itr = densities_.upperBound(point.x());

		//get the bounds for that bin - if at end, there is no lower bound in the current
		//container that is GREATER than the lower bound of the argument, so set to the last bin
		if(f_itr == densities_.end())
			--f_itr;

		DensitiesMap::iterator r_itr = f_itr;

		double lowerbin = f_itr.key();
		double upperbin = lowerbin + bin_size;

		//only add the current bin if it overlaps - might have to move first to get a bin that overlaps
		if(upperbin > other_lowerbin and lowerbin < other_upperbin) {

			bounds.append(BinBoundaries(lowerbin,upperbin));
			density_keys.append(f_itr.key());
		}

		//go forward until either at end or move past 'other' upper bound
		DensitiesMap::iterator itr_end = densities_.end();
		++f_itr; //need to move ahead first to check for end
		while(f_itr != itr_end and upperbin < other_upperbin) {

			lowerbin = (f_itr++).key();
			upperbin = lowerbin + bin_size;
			bounds.append(BinBoundaries(lowerbin, upperbin));
			density_keys.append(f_itr.key());
        }
		//now go the other way, stop at begin to avoid "issues"
		DensitiesMap::iterator itr_begin = densities_.begin();
		while(r_itr != itr_begin and lowerbin > other_lowerbin) {

			lowerbin = (--r_itr).key();
			upperbin = lowerbin + bin_size;
			bounds.append(BinBoundaries(lowerbin, upperbin));
			density_keys.append(r_itr.key());
		}

		//deal with first Density if necessary
		if(lowerbin > other_lowerbin) {

			lowerbin = r_itr.key();
			upperbin = lowerbin + bin_size;
			bounds.append(BinBoundaries(lowerbin, upperbin));
			density_keys.append(r_itr.key());
		}

		//now have all bin boundaries and their associated bin iterators, go through
		//and add the appropriate fraction of the current 'other' bin to them
		int bounds_size = bounds.size();

		for(int j = 0; j < bounds_size; ++j) {

			//want to determine how much of the current bounds lie within the 'other' bounds
			//choose the lowest upper bound and the highest lower bound
			double lowest_upper, highest_lower;

			if(bounds[j].first > other_lowerbin)
				highest_lower = bounds[j].first;
			else
				highest_lower = other_lowerbin;

			if(bounds[j].second < other_upperbin)
				lowest_upper = bounds[j].second;
			else
				lowest_upper = other_upperbin;

			//get the fraction of 'other', should obviously be a max of 1
			//also, if negative, means there is zero overlap
			double fraction = (lowest_upper - highest_lower) / other.binSize();
			if(fraction > 1)
				fraction = 1.0;

			else if(fraction < 0)
				fraction = 0;

			//add the fraction of the density to the current bin
			densities_[density_keys[j]] += fraction * point.y();
		}
	}
    *-------------------------------------------------------------------------------------------
    */

	//if the current NormalDensity now has > 1 point, need to check on bin size
   if(not bin_size_updated and densities_.size() > 1)
	   setBinSize();

	//set the flag that the containers are no longer synchronized
	all_densities_updated = false;
}

/*
 * Calculates and returns the geometric centre point along the normal axis for the density distribution.
 */
double NormalDensity::densityCentre() {

	//just need to find the average along the normal axis weighted by the density
	double total = 0, total_weight = 0;
	for(QPointF point = firstPoint(); not densitiesAtEnd(); point = nextPoint()) {

		total += point.x() * point.y();
		total_weight += point.y();
	}

	return total/total_weight;
}

//Calculates the symmetric, electron and neutron densities; should only be called after all number densities have been added
void NormalDensity::calcAllDensities(bool traj_data) {
   
	//clear the vectors first to make sure starting from a clean slate
	vec_densities.clear();
	sym_densities.clear();
	neut_densities.clear();
	elec_densities.clear();

	//go through vector for densities and create points by multiplying by
	//electron number and neutron SL
	DensitiesMap::iterator itr, r_itr = densities_.end(), itr_end = densities_.end();
	--r_itr; //set reverse iterator to the last value

	for(itr = densities_.begin(); itr != itr_end; ++itr) {

		//TODO: include the following only for traj data: + 0.5 * bin_size;
		double cur_x = itr.key();
		if(traj_data)
			cur_x += 0.5 * bin_size;

		double cur_d = itr.value();

		//standard density container
		vec_densities.append(QPointF(cur_x, cur_d));

		//for symmetrized
		double cur_d_sym = (cur_d + r_itr.value())/2;
		sym_densities.append(QPointF(cur_x, cur_d_sym));

		double cur_ed = cur_d * num_electrons;
		double cur_nsld = cur_d * neutron_sl;

		elec_densities.append(QPointF(cur_x, cur_ed));
		neut_densities.append(QPointF(cur_x, cur_nsld));

		//move reverse iterator
		--r_itr;
	}

	//everything is updated and synched
	all_densities_updated = true;
}

/*
 * Deletes all density values from all of the containers
 */
void NormalDensity::clearDensities() {

	//The map of course
	densities_.clear();

	vec_densities.clear();
	sym_densities.clear();
	elec_densities.clear();
	neut_densities.clear();

	//technically, the densities are up to date, but no more bin size
	all_densities_updated = true;
	bin_size_updated = false;
}

//Scales all of the positions along normal axis by 'factor'.
void NormalDensity::scaleIndependentAxis(double factor) {

	//Take all the data out and clear the container
	QList<float> density_keys = densities_.keys();
	QList<float> density_values = densities_.values();
	densities_.clear();

	//go through and put the scaled values back in
	int size = density_keys.size();
	for(int i = 0; i < size; ++i) {

		double scaled_key = density_keys[i] * factor;
		densities_.insert(scaled_key, density_values[i]);
	}

	//need to adjust the bin size
	setBinSize();

	all_densities_updated = false;
}

//Scales all of the values along normal axis by 'factor'.
void NormalDensity::scaleDependentAxis(double factor) {

	//Take all the data out and clear the container
	QList<float> density_keys = densities_.keys();
	QList<float> density_values = densities_.values();
	densities_.clear();

	//go through and put the scaled values back in
	int size = density_keys.size();
	for(int i = 0; i < size; ++i) {

		double scaled_value = density_values[i] * factor;
		double density_key = density_keys[i];
		densities_.insert(density_key, scaled_value);
	}

	//densities no longer up to date
	all_densities_updated = false;
}

//Shift all the density points along the independent axis
void NormalDensity::shiftIndependentAxis(double to_shift) {

	//Take all the data out and clear the container
	QList<float> density_keys = densities_.keys();
	QList<float> density_values = densities_.values();
	densities_.clear();

	//go through and put the scaled values back in
	int size = density_keys.size();
	for(int i = 0; i < size; ++i) {

		double shifted_key = density_keys[i] + to_shift;
		densities_.insert(shifted_key, density_values[i]);
	}

	all_densities_updated = false;
}

//Removes all data below the lower cut-off and above the upper
void NormalDensity::trimIndependentAxis(float lower, float upper) {

	//iterate through from beginning, get rid of all points smaller than 'lower'
	DensitiesMap::iterator itr_low = densities_.begin();
	DensitiesMap::iterator itr_high = --densities_.end();

	//use a float tolerance to avoid floating point silliness
	float tolerance = bin_size * 1e-3;

	while(itr_low.key() < lower and abs(itr_low.key() - lower) > tolerance) {

		//otherwise not there yet, delete and set the iterator to the next value
		itr_low = densities_.erase(itr_low);
	}

	//do same for upper
	while(itr_high.key() > upper and abs(itr_high.key() - upper) > tolerance) {

		itr_high = densities_.erase(itr_high);
		--itr_high; //since it should be 'end'
	}
}

//sets the atomic info
void NormalDensity::setAtomInfo(AtomType type, double ed, double nsld) {
   
   this->atom_type= type;
   num_electrons = ed;
   neutron_sl = nsld;
}

//helper function that sets the bin size; assumes there are at least two entries in the map
void NormalDensity::setBinSize() {

	DensitiesMap::iterator itr = densities_.begin();
	double first_position = itr.key();
	double second_position = (++itr).key();
	bin_size = second_position - first_position;

	bin_size_updated = true;
}
