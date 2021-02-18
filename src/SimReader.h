/*
 * SimReader.h
 *
 *  Created on: Nov 14, 2012
 *      Author: bholland
 */

#ifndef SIMREADER_H_
#define SIMREADER_H_

#include "NormalDensity.h"
#include "Component.h"
#include "FileCode.h"

#include <QtCore/QList>
#include <QtCore/QString>
#include <QtCore/QStringList>

typedef QMap<uint,QString> ParticleOrderMap;
typedef QMap<QString,NormalDensity> NormalDensityMap;
typedef QList<Component> Components;
typedef QList<NormalDensity> NormalDensities;

class SimReader {

	double area; //surface area of lipids
	uint nlip; //number of lipids
	bool normalize_volume; //true if normalization set

	double thickness;

	//Note: originally used a simple QList to store Densities, but wanted to avoid slow
	//lookup (i.e. O(N)) for a Density through the name, so use a map to order the names,
	//the side effect is the order the particles were entered is lost. Can create a second
	//map to track this, but does the order really matter, i.e. is alphabetical okay?

	QStringList first_line; //stores the names read in from the first line of the file
	NormalDensityMap sim_data; //data read in from the sim file
	NormalDensityMap::iterator sim_data_itr;
	NormalDensity dummy_density; //for the iterators, not to be used

	Components components; //data read in from cmp file

public:

	SimReader();
	virtual ~SimReader();

	double lipidArea() const;
	uint numLipids() const;

	NormalDensity& getParticleDensity(const QString&);
	const NormalDensity& getParticleDensity(const QString&) const;

	NormalDensities getParticleDensities() const;
	void clearSimData();

	NormalDensity* getFirstParticleDensity();
	NormalDensity* getNextParticleDensity();
	bool isEmpty() const;
	bool atEnd() const;

	bool particleLoaded(const QString&) const;

	Components& getComponents();
	const Components& getComponents() const;

	FileCode readSimFile(const QString&);
	FileCode readComponentFile(const QString&);

private:

	void normalizeData();
};

#endif /* SIMREADER_H_ */
