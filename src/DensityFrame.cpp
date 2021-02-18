/*
 * Frame.cpp
 *
 *  Created on: Mar 13, 2013
 *      Author: bholland
 */

#include "DensityFrame.h"

#include <QtAlgorithms>
#include <QtDebug>

#include <iostream>
#include <fstream>
#include <cfloat>

//how much of the plot past the bilayer to keep
//TODO: find a physically relevant value for this, currently zero is only value that has
//      any physical relevance, i.e. only concerned about scattering of bilayer and water included in it.
const double BILAYER_MARGIN = 2;
const int NUMBER_OF_AXES = 3;

using namespace std;

//Default constructor
DensityFrame::DensityFrame() {

	step_ = 0;
	sim_time = 0;
	normal_axis = AXIS_UNKNOWN;
}

//Constructor that sets the step number for the frame, as well as the time for the simulation
//in picoseconds (i.e. NOT real time)
DensityFrame::DensityFrame(uint step, double time, Axis axis) {

	step_ = step;
	sim_time = time;
	normal_axis = axis;
}

//Destructor
DensityFrame::~DensityFrame() {

	//go through and delete all the NormalDensity objects in the container
	//qDeleteAll(normal_densities);

	//Actually don't, but NEED to delete these upon destruction of DensityFrame object, however
	//it happens... this is really annoying and bad coding practice, but really want to only
	//keep the pointers since making many copies of NormalDensity objects could be very slow
}

//getters
uint DensityFrame::step() const {return step_;}
double DensityFrame::time() const {return sim_time;}
const MDBox& DensityFrame::simBox() const {return sim_box;}
Axis DensityFrame::normalAxis() const {return normal_axis;}

//setters
void DensityFrame::setStep(uint step) {

	step_ = step;
}

void DensityFrame::setSimTime(double time) {

	sim_time = time;
}

void DensityFrame::setSimBox(const MDBox& box) {

	sim_box = box;
}

void DensityFrame::setNormalAxis(Axis axis) {

	normal_axis = axis;
}

//checks if the key is in the container
bool DensityFrame::inNormalContainer(const ParticleLabel& label) const {

	if(normal_densities.contains(label))
		return true;
	return false;
}

//returns the number of NormalDensities in the map
uint DensityFrame::numNormalDensities() const {

	return normal_densities.size();
}

//bool DensityFrame::inPolarContainer(const ParticleLabel& label) const {
//
//	if(polar_pds.contains(label))
//		return true;
//	return false;
//}

//returns the normal ParticleDensity object for the given ParticleLabel.
//if 'pl' isn't contained in the map, it will be added and an empty ParticleDensity object will be returned
//This is bad, do NOT want it added this way, so make sure it is in the map before calling this
NormalDensity& DensityFrame::normalDensity(const ParticleLabel& pl) {

	//check if it contains it first, return error message for debugging if not - should still be faster than making a copy
	if(not normal_densities.contains(pl))
		cerr << "DensityFrame::normalDensity: label not in map, default NormalDensity created and added - not good!" << endl;

	return normal_densities[pl];
}

//QMap won't return a constant reference, kind of annoying, but don't need to worry about adding a default object to the container at least
NormalDensity DensityFrame::normalDensity(const ParticleLabel& pl) const {

	return normal_densities.value(pl);
}

//Iterator method for getting NormalDensity objects from map
NormalDensity* DensityFrame::firstNormalDensity() {

	//set iterator to beginning
	normal_map_itr = normal_densities.begin();

	//check to see if frame is empty - if so it will crash, so return a debug message
	if(normal_map_itr == normal_densities.end())
		cerr << "DensityFrame::firstNormalDensity: Map is empty, nothing to return" << endl;

	return &normal_map_itr.value();
}

//Returns next NormalDensity, but does not check if at end, use 'normalAtEnd()' as check
NormalDensity* DensityFrame::nextNormalDensity() {

	//first move ahead
	++normal_map_itr;

	//check if at end, if so return a dummy value - this will never get return in a properly constructed for loop
	if(normal_map_itr == normal_densities.end())
		return &dummy_density;

	return &normal_map_itr.value();
}

//checks to see if iterator is at end of container
bool DensityFrame::normalAtEnd() const {

	if(normal_map_itr == normal_densities.end())
		return true;
	return false;
}

//return the polar ParticleDensity object for the given ParticleLabel
//if 'pl' isn't contained in the map, it will be added and an empty ParticleDensity object will be returned
//PolarDensity& DensityFrame::polarParticleDensity(const ParticleLabel& pl) {
//
//	return radial_pds.value(pl);
//}

//Adds the given ParticleDensity for the associated ParticleLabel.  If the ParticleLabel
//already exists in the container, does nothing.
void DensityFrame::addNormalDensity(const ParticleLabel& label, const NormalDensity& pd) {

	if(not normal_densities.contains(label)) {
		normal_densities.insert(label, pd);
	}
}

//Deletes all of the objects pointed to by the NormalDensity pointers, and clears the map
//completely.
void DensityFrame::clearNormalDensities() {

	normal_densities.clear();
}

//same as above for radial densities
//void DensityFrame::addPolarDensity(const ParticleLabel& label, PolarDensity* pd) {
//
//	if(not polar_pds.contains(label)) {
//		polar_pds.insert(label, pd);
//	}
//}

//Returns a list of the NormalDensity objects, i.e. no pointers
NormalDensities DensityFrame::normalDensities() const {

	return normal_densities.values();
}

//Return the entire map of densities... intended to be used by other DensityFrame objects
NormalDensityLabelMap& DensityFrame::normalDensityMap() {return normal_densities;}
const NormalDensityLabelMap& DensityFrame::normalDensityMap() const {return normal_densities;}

//Handy function for adding the data from another DensityFrame object to the current one. Will look
//to see if the Particle exists in this frame, and if so adds to it. If the current DensityFrame
//is completely empty, the list is simply added to the current object
void DensityFrame::addNormalDensityFrame(const DensityFrame& other) {

	bool normal_empty = normal_densities.isEmpty();
//	bool polar_empty = polar_pds.isEmtpy();

	//first check to see if the current lists are empty
	if(normal_empty) {

		//set and return, TODO: add polar as well, will need to adjust this region a bit
		normal_densities = other.normalDensityMap();
		return;
	}

	//loop through the 'other' DensityFrame and look for the Particles
	NormalDensityLabelMap::const_iterator itr, itr_end = other.normalDensityMap().end();
	NormalDensityLabelMap::iterator found;

	for(itr = other.normalDensityMap().begin(); itr != itr_end; ++itr) {

		ParticleLabel label = itr.key();

		found = normal_densities.find(label);
		if(found != normal_densities.end()) {

			//if ParticleLabel found, add the normal density from other to the current
			found.value().addDensities(itr.value());
		}
	}
}

/*
 * Scales the NormalDensity objects in the current frame using the provided axis and average thickness
 *  along that axis so that the densities are properly comparable.
 *
 *  The reason for doing this is that in an experiment, the volume fluctuations of the system will actually be
 *  negligible.  But due to the small size of simulated systems, the fluctuations become large enough relative
 *  to the box size that they need be considered.
 */
void DensityFrame::scaleNormalAxis(double avg_thickness) {

	//the normal axis must be set, just make it a 'cerr' as program will crash anyway
	if(normal_axis == AXIS_UNKNOWN) {

		cerr << "Error - DensityFrame::scaleNormalAxis: the normal axis is not set!" << endl;
	}

	double scaling_factor = avg_thickness / sim_box.realDiagonal()[normal_axis];

	for(auto& density : normal_densities) {

		density.scaleIndependentAxis(scaling_factor);
	}
}

//Scales all NormalDensity objects by the inverse of the given value.  Intended to be used on a DensityFrame
//that has been calculated as an average over numerous DensityFrame objects
void DensityFrame::scaleNormalValues(uint factor) {

	double scaling_factor = 1.0 / (double) factor;

	for(auto& density : normal_densities) {

		density.scaleDependentAxis(scaling_factor);
	}
}

//Normalizes the density values by the slab volume.
void DensityFrame::normalizeNormalValues(double bin_size) {

	//the normal axis must be set, just make it a 'cerr' as program will crash anyway
	if(normal_axis == AXIS_UNKNOWN) {

		cerr << "Error - DensityFrame::scaleNormalAxis: the normal axis is not set!" << endl;
	}

	//for now the assumption is that the 'box' is in fact just a normal box (only diagonal in the tensor)
	Vector3D box_dimensions = sim_box.realDiagonal();

	//want the two values that are NOT the axis
	double dim1 = 0, dim2 = 0;
	for(int i = 0; i < NUMBER_OF_AXES; ++i) {

		if(i == normal_axis)
			continue;

		if(dim1 == 0)
			dim1 = box_dimensions[i];
		else
			dim2 = box_dimensions[i];
	}

	//get the factor and go through all densities
	double norm_factor = 1.0 / (bin_size * dim1 * dim2);

	for(auto& density : normal_densities) {

		density.scaleDependentAxis(norm_factor);
	}
}

/*
 * Assumes that the bilayer was centered prior to analysis with SIMtoEXP and centers the data about zero using
 * the box size
 * TODO: this is for GROMACS only, NAMD does not reset everything > 0
 */
void DensityFrame::centerNormalValues() {

	//the normal axis must be set, just make it a 'cerr' as program will crash anyway
	if(normal_axis == AXIS_UNKNOWN) {

		cerr << "Error - DensityFrame::scaleNormalAxis: the normal axis is not set!" << endl;
	}

	//Find the center of all lipid particles and use this as the shift
	double total = 0;
	int num_particles = 0;

	for(auto& density : normal_densities){

		ResidueType type = density.residueType();
		if(type != WATER and type != ION and type != HEAVY_WATER) {

			total += density.densityCentre();
			++num_particles;
		}
	}

	double to_shift = -total/num_particles;

	//go through all the densities again and shift them

	for(auto& density : normal_densities) {

		density.shiftIndependentAxis(to_shift);
	}
}

//TODO: this is unnecessary and should really be changed to only trim the erroneous low density
//      water at the box edges - the amount of bulk water has not effect on the structure factors
//      NOT CURRENTLY CALLED
void DensityFrame::trimNormalDensities() {

	//find the min and max independent variables that have non-zero values
	//out of all residues that are NON-water
	double max = -DBL_MAX;
	double min = DBL_MAX;

	for(auto& density : normal_densities) {

		//if not part of solution, look for ends
		if(density.residueType() != WATER and density.residueType() != ION and density.residueType() != HEAVY_WATER) {

			QPointF first = density.firstNonZeroPoint();
			QPointF last = density.lastNonZeroPoint();
			double cur_min = first.x();
			double cur_max = last.x();

			if(cur_min < min and first.y() > 0)
				min = cur_min;

			if(cur_max > max and last.y() > 0)
				max = cur_max;
		}
	}

	//assume centered about zero (for now anyway) and only keep the value farthest from zero to keep symmetric
	if(abs(min) > abs(max))
		max = abs(min);
	else
		min = -max;

	//now add margins
	min -= BILAYER_MARGIN;
	max += BILAYER_MARGIN;

	//go through and trim the NormalDensity objects
	for(auto& density : normal_densities) {

		density.trimIndependentAxis(min, max);
	}
}

/*
 * Takes the hydrogen from water in the system and reduces it by (1 - heavy_fraction) while adding the
 * requested amount of deuterium. Assumes of course that there is in fact water in the system!
 *
 * This will only change hydrogen to deuterium, so if the system was actually just simulated with D2O,
 * simply leave the fraction to zero!
 */
void DensityFrame::adjustDeuteriumFraction(double heavy_fraction, AtomicInfo* info, SimGranularity granularity, bool isElemental) {

	//TODO: should really change particle labels to include AtomType and ResidueType enumerators,
	//      and use these below instead of the residue and particle names, much more robust and safer
	//      currently do not use AtomicInfo, but it is included for this purpose

	//if the fraction is zero, move on with our lives
	if(heavy_fraction < 1E-12)
		return;

	//boolean for error reporting
	bool water_hydrogen_found = false;

	QString deut_name = "D-W";
	int particle_number = 1;

	//first get the hydrogen from water in the system, if not there output an error message and return
	for(NormalDensity* density = firstNormalDensity(); not normalAtEnd(); density = nextNormalDensity()) {

		ParticleLabel label = density->particleLabel();

		QString res_name = density->residueName().toLower();
		QString particle_name = density->particleName().toLower();

		//if the simulation is Martini based, changes are based on water beads
		if(granularity == MARTINI) {

			if(res_name == "w" or res_name == "pw") {

				//create a deuterium density with correct fraction and add it to the frame
				NormalDensity deuterium = *density;
				deuterium.setParticleName("D2O");
				deuterium.setResidueName("D2O");
				deuterium.setAtomInfo(MART_D2O, info->numberElectrons(MART_D2O), info->neutronSL(MART_D2O));
				deuterium.setResidueType(HEAVY_WATER);

				deuterium.scaleDependentAxis(heavy_fraction);
				normal_densities.insert(deuterium.particleLabel(), deuterium);

				//adjust the fraction for the water - don't delete the hydrogen though to allow the fraction to be changed later
				density->scaleDependentAxis(1 - heavy_fraction);

				water_hydrogen_found = true;
			}

		} else if(granularity == ATOMISTIC) {

			if(res_name == "sol" or res_name == "h2o" or res_name == "water") {

				//this should work for both elemental and coordinate file particle names
				if(particle_name[0] == 'h') {

					//create a deuterium density with correct fraction and add it to the frame
					NormalDensity deuterium = *density;

					if(isElemental)
						deuterium.setParticleName("D");

					else {

						QString cur_name = deut_name + QString::number(particle_number);
						deuterium.setParticleName(cur_name);
						particle_number++;
					}

					deuterium.setResidueName(density->residueName());
					deuterium.setAtomInfo(DEUTERIUM, info->numberElectrons(DEUTERIUM), info->neutronSL(DEUTERIUM));
					deuterium.setResidueType(HEAVY_WATER);

					deuterium.scaleDependentAxis(heavy_fraction);
					normal_densities.insert(deuterium.particleLabel(), deuterium);

					//adjust the fraction for the water if necessary, otherwise delete the H2O if fraction = 1
					density->scaleDependentAxis(1 - heavy_fraction);

					water_hydrogen_found = true;

					//if elemental, found the water hydrogen, can move on
					if(isElemental)
						return;
				}
			}
		}
	}

	//if no water hydrogen found, report an error
	if(not water_hydrogen_found) {

		cerr << "DensityFrame::adjustDeuteriumFraction - no hydrogen founds as part of a water molecule, no deuterium added" << endl;
	}
}

/*
 * For debugging purposes only... puts all the NormalDensity info into a single string for dumping to a file
 */
QString DensityFrame::dumpNormalDensityInfo() {

	QString to_dump;

	//go through all the NormalDensity objects
	for(NormalDensity* density = firstNormalDensity(); not normalAtEnd(); density = nextNormalDensity()) {

		to_dump += density->particleName();
		to_dump += "\n";

		//output all density point values
		for(QPointF p = density->firstPoint(); not density->densitiesAtEnd(); p = density->nextPoint()) {

			to_dump += QString::number(p.x());
			to_dump += "  ";
			to_dump += QString::number(p.y());
			to_dump += "\n";
		}

		to_dump += "\n";
	}

	return to_dump;
}
