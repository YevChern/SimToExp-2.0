/*
 * System.cpp
 *
 *  Created on: Feb 25, 2013
 *      Author: bholland
 */

//project includes
#include "System.h"

//Qt includes
#include <QtDebug>

//c++ includes
#include <iostream>

using namespace std;

//Default constructor
System::System() {

	cur_particle = 0;

	//start traversing itr at end
	res_itr = residues_.end();
	dummy_res.setName("Dummy residue");
}

//Constructor that sets the system name upon creation
System::System(const QString& name) {

	cur_particle = 0;
	name_ = name;
	res_itr = residues_.end();

	dummy_res.setName("Dummy residue");
}

//Destructor - don't worry about particle pointers, they are in a separate container anyway
System::~System() {}

//getters
const QString& System::name() const {return name_;}

//To return references in a QMap, need to get iterator first
Residue& System::operator [](uint index) {

	ResidueMap::iterator itr = residues_.find(index);
	if(itr == residues_.end())
		cerr << "System::operator []: residue index " << index << " not found" << endl;

	return itr.value();
}

const Residue& System::operator [](uint index) const {

	ResidueMap::const_iterator itr = residues_.find(index);
	if(itr == residues_.end())
		cerr << "System::operator []: residue index " << index << " not found" << endl;

	return itr.value();
}

//Iterator methods for traversing the map, useful if either don't know the ids, also
//O(N) for traversal instead of O(NlogN) from the [] operator
Residue& System::firstResidue() {

	//set the iterator, check to see if at end already and output debugging error if so
	res_itr = residues_.begin();
	if(res_itr == residues_.end())
		cerr << "System::firstResidue(): system is empty, no residues." << endl;

	//move iterator to "next" after returning
	return res_itr.value();
}

Residue& System::nextResidue() {

	++res_itr;

	//If already at end, return a dummy residue
	if(res_itr == residues_.end())
		return dummy_res;

	//return value
	return res_itr.value();
}

//bool to check if at end of map - this is for the traversing iterator, not the addition iterator
//Guaranteed to be at 'end' upon construction and after addition of residue
bool System::atEnd() const {

	if(res_itr == residues_.end())
		return true;
	return false;
}

//just sets the name
void System::setName(const QString& name) {

	name_ = name;
}

//creates and adds a Residue, returning the object just created by reference
Residue& System::addResidue(const QString& name, uint id, ResidueType lipid) {

	//each time a residue is added, need to reset the traversing iterator to the 'new' end
	//not sure if this is necessary, 'end' could be a special case that doesn't
	//get invalidated, but do it anyway
	Residue& res = residues_.insert(id, Residue(name,id)).value();
	res_itr = residues_.end();

	return res;
}

/*
 * Add a particle pointer to end of the container - should be added in same order as in coord file
 */
void System::addParticle(Particle* p) {

	particles_.append(p);
}

//Checks to see if there are residues in the container
bool System::isEmpty() const {return residues_.empty();}

//goes through all residues and adds up the number of particles from each
uint System::numberOfParticles() const {

	uint num_particles = 0;

	for(const Residue& res : residues_)
		num_particles += res.numberOfParticles();

	return num_particles;
}

//return the number of residues in the system
uint System::numberOfResidues() const {

	return residues_.size();
}

//Return the MDBox object
const MDBox& System::systemBox() const {return box;}

//Only sets the diagonal of the box, the rest of the tensor is set to zero (i.e. it is a true 'box')
void System::setBoxDimensions(double x, double y, double z) {

	box.setBoxDimensions(x, y, z);
}

//sets the box dimensions to the given values
void System::setBoxDimensions(const Vector3D& x, const Vector3D& y, const Vector3D& z) {

	box.setBoxDimensions(x, y, z);
}

//return the box dimensions
const Matrix3D& System::boxDimensions() const {return box.realSpace();}

//clears entire system of all information
void System::clearSystem() {

	name_ = "";
	residues_.clear(); //no dynamically allocated values

	//reset dimensions
	box = MDBox();
}

//slots for receiving data updates
//updates the MDBox with the current value - only changes for NPT ensemble of course
void System::updateMDBox(const Matrix3D& box) {

	this->box.setBoxDimensions(box);
}

//updates the current coordinate
void System::updateNextCoord(const Vector3D& coord) {

	//check to see if current particle is greater than total number, if so reset to zero
	if(cur_particle == particles_.size())
		cur_particle = 0;

	//Update the current coordinate
	particles_[cur_particle]->setCoord(coord);

	//move index ahead for next time
	cur_particle++;
}
