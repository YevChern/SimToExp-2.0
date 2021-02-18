/*
 * ParticleDensities.cpp
 *
 *  Created on: Mar 10, 2013
 *      Author: bholland
 */

#include "ParticleIntensities.h"

//Default constructor
ParticleIntensities::ParticleIntensities() {}

//Destructor
ParticleIntensities::~ParticleIntensities() {}

//Returns the intensity for the given Particle; CAUTION: quietly returns zero if the Particle does not exist
uint ParticleIntensities::intensity(const ParticleLabel& label) const {

	return intensities.value(label);
}

//Increments the intensity of the given AtomType by 'to_incr'; default is 1
void ParticleIntensities::incrementIntensity(const ParticleLabel& label, uint to_incr) {

	//check for 'type' first, if not there yet, just set to 'to_incr'
	QMap<ParticleLabel,uint>::iterator itr = intensities.find(label);

	if(itr == intensities.end())
		intensities[label] = to_incr;
	else
		itr.value() += to_incr;
}

/*
 * Returns the number of particles (i.e. ParticleLabels) in the container
 */
uint ParticleIntensities::numParticles() {

	return intensities.size();
}

//Clears the entire intensity map
void ParticleIntensities::clear() {

	intensities.clear();
}
