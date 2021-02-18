/*
 * ParticleIntensities.h
 *
 *  Created on: Mar 10, 2013
 *      Author: bholland
 *
 * For keeping track of intensities read from trajectories; intended to be used
 * for a single 'bin', i.e. spatial container along z or in (x,y)
 */

#ifndef PARTICLEINTENSITIES_H_
#define PARTICLEINTENSITIES_H_

#include <QtCore/QMap>

#include "ParticleLabel.h"

typedef QMap<ParticleLabel,uint> IntensityMap;

class ParticleIntensities {

	IntensityMap intensities;

public:

	ParticleIntensities();
	virtual ~ParticleIntensities();

	uint intensity(const ParticleLabel&) const;
	void incrementIntensity(const ParticleLabel&, uint to_incr = 1);

	uint numParticles();

	void clear();
};

#endif /* PARTICLEINTENSITIES_H_ */
