/*
 * Particle.cpp
 *
 *  Created on: Feb 25, 2013
 *      Author: bholland
 */

#include "Particle.h"

using namespace std;

//Default constructor
Particle::Particle() {

	type = UNKNOWN;
	id_ = 0;
}

//Constructor that sets the coordinate file name and id
Particle::Particle(const QString& name, uint id) {

	coord_file_name = name;
	type = UNKNOWN;
	id_ = id;
}

//Destructor
Particle::~Particle() {}

//getters
const QString& Particle::coordFileName() const {return coord_file_name;}
const QString& Particle::recognizedName() const {return recognized_name;}
AtomType Particle::atomType() const {return type;}
uint Particle::id() const {return id_;}
const Vector3D& Particle::coord() const {return coord_;}

//setters
void Particle::setAtomType(AtomType type) {this->type = type;}
void Particle::setID(uint id) {id_ = id;}

void Particle::setCoordFileName(const QString& name) {

	coord_file_name = name;
}

void Particle::setRecognizedName(const QString& name) {

	recognized_name = name;
}

//set the coordinate using individual values
void Particle::setCoord(double x, double y, double z) {

	coord_.set(x,y,z);
}

//set the coordinates using a 3D vector
void Particle::setCoord(const Vector3D& coord) {

	coord_.set(coord);
}


//Less-than operator, needed for QMap in Residue. Order determined only by ID.
bool Particle::operator <(const Particle& other) const {

	if(id_ < other.id())
		return true;
	return false;
}

//equality operators - equal iff name, id and type are the same... does NOT check coord
bool Particle::operator ==(const Particle& other) const {

	if(recognized_name == other.recognizedName() and type == other.atomType() and id_ == other.id())
		return true;
	return false;
}

bool Particle::operator !=(const Particle& other) const {

	return not (*this == other);
}
