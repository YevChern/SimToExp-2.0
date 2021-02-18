/*
 * NormalGrid.h
 *
 *  Created on: Mar 10, 2013
 *      Author: bholland
 *
 *  A 1D "grid" for intensity values in the direction normal to the plane of the bilayer (as defined by the
 *  user).  User provides the bin size, and the grid will grow automatically according to the points that are
 *  in the coordinate file.
 *
 */

#ifndef NORMALGRID_H_
#define NORMALGRID_H_

#include <QtCore/QMap>
#include <QtCore/QPointF>

#include "ParticleIntensities.h"
#include "AtomicInfo.h"
#include "NormalDensity.h"

typedef QMap<double,ParticleIntensities> GridIntensityMap;

class NormalGrid {

	//grid is defined as the center of the bins! Use 'dx' to determine if within the bin
	double dx;

	GridIntensityMap norm_grid;
	GridIntensityMap::iterator norm_grid_itr; //intended for testing use only

	AtomicInfo* info;

public:

	NormalGrid();
	NormalGrid(AtomicInfo*);
	virtual ~NormalGrid();

	void setAtomicInfo(AtomicInfo*);

	double max() const;
	double min() const;
	NormalDensity normalDensity(const ParticleLabel&, AtomType);

	double binSize() const;
	void setBinSize(double);
	uint numBins() const;

	void incrementIntensity(const ParticleLabel&, double, uint to_incr = 1);
	uint intensity(const ParticleLabel&, double) const;

	void resetIntensities();
	void clearGrid();

	//mainly intended for testing
	const ParticleIntensities& firstParticleInt();
	const ParticleIntensities& nextParticleInt();
	bool atEnd() const;

private:

	QPointF probDensity(const ParticleLabel&,const GridIntensityMap::iterator&) const;
};

#endif /* NORMALGRID_H_ */
