/*
 * Frame.h
 *
 *  Created on: Mar 13, 2013
 *      Author: bholland
 *
 *  Intended to represent the calculated densities for a 'frame' of an extended trajectory file.  This way,
 *  Frame objects are calculated one at a time and then stored, allowing the remainder of the coordinate data to be overwritten.  This
 *  removes the need to load the entire trajectory into memory, and should allow for the analysis of very large file sizes.
 */

#ifndef DENSITYFRAME_H_
#define DENSITYFRAME_H_

#include <QtCore/QMap>

#include "NormalDensity.h"
#include "ParticleLabel.h"
#include "AtomicInfo.h"
#include "MDBox.h"
#include "Axis.h"
#include "SimGranularity.h"

typedef QList<NormalDensity> NormalDensities;
typedef QMap<ParticleLabel,NormalDensity> NormalDensityLabelMap;
typedef QMap<ResidueType,uint> LipidNumberMap;

class DensityFrame {

	//Values that define the frame
	uint step_;
	double sim_time;
	MDBox sim_box;
	Axis normal_axis;

	NormalDensityLabelMap normal_densities;
	NormalDensityLabelMap::iterator normal_map_itr;
//	QMap<ParticleLabel,PolarDensity> polar_pds;

	//dummy value for iterators
	NormalDensity dummy_density;

public:

	DensityFrame();
	DensityFrame(uint,double,Axis);
	virtual ~DensityFrame();

	uint step() const;
	double time() const;
	const MDBox& simBox() const;
	Axis normalAxis() const;

	void setStep(uint);
	void setSimTime(double);
	void setSimBox(const MDBox&);
	void setNormalAxis(Axis);

	bool inNormalContainer(const ParticleLabel&) const;
//	bool inPolarContainer(const ParticleLabel&) const;

	NormalDensity& normalDensity(const ParticleLabel&);
	NormalDensity normalDensity(const ParticleLabel&) const;

	NormalDensity* firstNormalDensity();
	NormalDensity* nextNormalDensity();
	bool normalAtEnd() const;

	uint numNormalDensities() const;

//	PolarDensity& polarParticleDensity(const ParticleLabel&);

	void addNormalDensity(const ParticleLabel&, const NormalDensity&);
//	void addPolarDensity(const ParticleLabel&, const PolarDensity&);

	NormalDensities normalDensities() const;
	void clearNormalDensities();

	NormalDensityLabelMap& normalDensityMap();
	const NormalDensityLabelMap& normalDensityMap() const;
//	const QMap<ParticleLabel,PolarDensity>& polarDensities();

	void addNormalDensityFrame(const DensityFrame&);

	void scaleNormalAxis(double);
	void scaleNormalValues(uint);

	void normalizeNormalValues(double);
	void centerNormalValues();
	void trimNormalDensities();

	void adjustDeuteriumFraction(double,AtomicInfo*,SimGranularity,bool);

	//debugging
	QString dumpNormalDensityInfo();
};

#endif /* DENSITYFRAME_H_ */
