/*
 * Residue.cpp
 *
 *  Created on: Feb 27, 2013
 *      Author: bholland
 */

#include "Residue.h"

#include <iostream>

#include <QtCore/QtDebug>

using namespace std;

//Default constructor
Residue::Residue() {

	name_ = "";
	res_id = 0;
	cur_particle = particles_.end();

	dummy_p.setRecognizedName("Dummy particle");
}

//Constructor that sets the name and id
Residue::Residue(const QString& name, uint id) {

	name_ = name;
	res_id = id;
	cur_particle = particles_.end();

	dummy_p.setRecognizedName("Dummy particle");
}

//Destructor
Residue::~Residue() {}

//getters
const QString& Residue::name() const {return name_;}
uint Residue::resId() const {return res_id;}

const Particle& Residue::particle (uint id) const {

	//Do an check and output an error message so if crashes due to null pointer, should know why
	if(not particles_.contains(id))
		cerr << "Residue::particle: particle ID not in residue, NULL pointer likely returned" << endl;

	return particles_.constFind(id).value();
}

//Creates and returns a ParticleLabel for the given particle index
QPair<ParticleLabel,AtomType> Residue::particleLabel(uint index) const {

	Particle particle = particles_[index];
	ParticleLabel label(name_, particle.coordFileName());
	return QPair<ParticleLabel,AtomType>(label, particle.atomType());
}

/*
 * Creates and returns a ParticleLabel for the given particle index that is common for
 * all particles with the same name.  For a ParticleLabel that distinguishes between different
 * particle of the same type, use 'particleLabel'
 */
QPair<ParticleLabel,AtomType> Residue::elementalLabel(uint index) const {

	Particle particle = particles_[index];
	ParticleLabel label(name_, particle.recognizedName());
	return QPair<ParticleLabel,AtomType>(label, particle.atomType());
}

//just return the number of particles in the residue
uint Residue::numberOfParticles() const {

	return particles_.size();
}

//setters
void Residue::setName(const QString& name) {name_ = name;}
void Residue::setResID(uint id) {res_id = id;}

//Adds a Particle to the Residue, returning the particle just added by reference
Particle& Residue::addParticle(const QString& name, uint id) {

	Particle& particle = particles_.insert(id, Particle(name,id)).value();

	//iterator is no longer valid, put at end
	cur_particle = particles_.end();

	return particle;
}

//Adds a copy of the given particle
void Residue::addParticle(const Particle& particle) {

	particles_.insert(particle.id(), particle);

	//iterator is no longer valid, put at end
	cur_particle = particles_.end();
}

//bool that check if there are any Particle objects in the residue
bool Residue::isEmpty() const {

	return particles_.isEmpty();
}

//Iterative functions that go through the residue in order
//-----------------------------------------------------
//Return the first particle in the map
Particle& Residue::firstParticle() {

	//set to beginning, if also end, report error for debugging but let crash
	cur_particle = particles_.begin();
	if(cur_particle == particles_.end())
		cerr << "Residue::firstParticle: residue is empty, nothing to return" << endl;

	//return and move the iterator ahead for the "next" particle
	return cur_particle.value();
}

//Return the next Particle iteration in the map
Particle& Residue::nextParticle() {

	//move particle ahead first
	++cur_particle;

	//if already at end, return a dummy particle
	if(cur_particle == particles_.end())
		return dummy_p;

	//return and move the iterator ahead for the "next" particle
	return cur_particle.value();
}

//Just checks to see if the iterator is at the end of the container
bool Residue::atEnd() const {

	if(cur_particle == particles_.end())
		return true;
	return false;
}

//Operators
//Array style getters - they use 'find' so they won't silently insert anything
//LogN time in the map, so use iterators if going through all Particles objects
Particle& Residue::operator [](uint index) {

	return particles_.find(index).value();
}

const Particle& Residue::operator [](uint index) const {

	return particles_.constFind(index).value();
}

//equality operator - equal iff name and res_id are same... does NOT check through particles
bool Residue::operator ==(const Residue& other) const {

	if(name_ == other.name() and res_id == other.resId())
		return true;
	return false;
}

bool Residue::operator !=(const Residue& other) const {

	return not (*this == other);
}
