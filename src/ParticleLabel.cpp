/*
 * ParticleLabel.cpp
 *
 *  Created on: Mar 13, 2013
 *      Author: bholland
 */

#include "ParticleLabel.h"

//Default constructor
ParticleLabel::ParticleLabel() {}

//constructor that sets both names in the label
ParticleLabel::ParticleLabel(const QString& res_name, const QString& particle_name) {

	this->res_name = res_name;
	this->particle_name = particle_name;
}

//Destructor
ParticleLabel::~ParticleLabel() {}

//getters
const QString& ParticleLabel::residueName() const {return res_name;}
const QString& ParticleLabel::particleName() const {return particle_name;}

//setters
void ParticleLabel::setResidueName(const QString& name) {

	res_name = name;
}

void ParticleLabel::setParticleName(const QString& name) {

	particle_name = name;
}

//less than operator - for sorting; first by residue, then by particle
bool ParticleLabel::operator<(ParticleLabel other) const {

	if(res_name < other.residueName())
		return true;
	else if(res_name > other.residueName())
		return false;

	//otherwise residue name are equal
	if(particle_name < other.particleName())
		return true;
	return false;
}

//equals operator, only here for completeness. ParticleLabels are equal
//iff both the residue and particle names are equal - case sensitive!
bool ParticleLabel::operator ==(ParticleLabel other) const {

	if(res_name == other.residueName() and particle_name == other.particleName())
		return true;
	return false;
}
