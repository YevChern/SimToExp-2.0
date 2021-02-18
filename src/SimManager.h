/*
 * SimManager.h
 *
 *  Created on: Oct 2, 2012
 *      Author: bholland
 */

#ifndef SIMMANAGER_H_
#define SIMMANAGER_H_

#include <QtCore/QMap>
#include <QtCore/QRectF>

#include "AtomicInfo.h"
#include "FourierFormFactor.h"
#include "SimReader.h"
#include "Axis.h"
#include "DensityFrame.h"

typedef QVector<double> Doubles;
typedef QVector<QPointF> Points;
typedef QList<AtomType> AtomTypes;
typedef QList<DensityFrame> DensityFrames;
typedef QVector<Doubles> ComponentMatrix;
typedef QMap<AtomType,Points> TotalDensities;
typedef QList<NormalDensity> NormalDensities;
typedef QList<Component> Components;

const uint elemental_index = 0;
const uint particle_index = 1;

class SimManager {

	SimReader reader;
	AtomTypes atom_types;

	DensityFrames ensemble_frames;
	DensityFrames time_function_frames;

	ComponentMatrix cmp_matrix; //for calculating volume
	Doubles cmp_vector; //for calculating volume

	TotalDensities total_elec_densities;
	TotalDensities total_neut_densities;

	FourierFormFactor* xray_fourier_FF;
	FourierFormFactor* neut_fourier_FF;

	//handy to have pointer
	AtomicInfo* info;

	bool sim_data_loaded;
	bool sim_data_from_traj;
	bool form_factors_calculated;
	bool atomTypesSet;
	bool compFileLoaded;
	bool volCalculated;

	QRectF num_densities_rect;
	QRectF ed_rect;
	QRectF nsld_rect;
	QRectF xrayFF_rect;
	QRectF neutFF_rect;
	
	double volume_rms;

	//values to be set if read from trajectory
	double traj_lipid_area;
	uint traj_num_lipids;

public:

	SimManager();
	virtual ~SimManager();

	void setAtomInfo(AtomicInfo&); //also create FF objects

	FileCode readSimFile(const QString&);
	FileCode readComponentFile(const QString&);
	void setSimDataFromTrajectory(const DensityFrames&);
//	void setSimDataFromTrajectory(const DensityFrames&);
	bool dataLoaded() const;

	bool simDataIsFromTrajectory() const;
	double lipidArea() const;
	uint numLipids() const;
	void setLipidArea(double);
	void setLipidNumber(uint);

	NormalDensities getElementalData() const;
	NormalDensities getParticleData() const;

	Components& getComponents();
	const Components& getComponents() const;
	uint numberComponents() const;
	Points getTotalComponentProbs();

	const QRectF& numDensitiesRect();
	const QRectF& electronDensityRect();
	const QRectF& neutronSLDensityRect();
	QRectF xrayFormFactorRect() const;
	QRectF neutronFormFactorRect() const;

	const AtomTypes& getAtomTypes();
	uint numberAtomTypes() const;
	bool particleLoaded(const QString&) const;

	Points& getTotalElectronDensity(AtomType);
	const Points& getTotalElectronDensity(AtomType) const;
	Points& getTotalNeutronDensity(AtomType);
	const Points& getTotalNeutronDensity(AtomType) const;

	//for setting up fourier transforms
	void setXrayQvalues(double,double,double);
	void setNeutronQvalues(double,double,double);

	//sets both transforms
	void setSymmetrized(bool);
	void setSolventScatteringDensities(double,double);

	void calcFourierTransforms(QWidget*);

	Points& getXrayFormFactors();
	const Points& getXrayFormFactors() const;
	Points& getNeutronFormFactors();
	const Points& getNeutronFormFactors() const;

	bool formFactorsCalculated() const;
	void clearFourier();
	void removeData();

	bool componentFileLoaded() const;
	bool volumeCalculated() const;
	
	void solveVolumeEquations();
	double volumeRMS() const;
	void removeComponents();

private:

	void calcTotalDensities();
	void calcTrajDataTotalDensities();
	void insertDensities(NormalDensity*);

	void alignDensityFrame(DensityFrame& frame);

	void setNormalDensityAtomInfo(NormalDensity&);

	void populateVolumeEquations();
	void calcVolumeRMS(const Doubles&);
};

#endif /* SIMMANAGER_H_ */
