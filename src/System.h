/*
 * System.h
 *
 *  Created on: Feb 25, 2013
 *      Author: bholland
 */

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <QtCore/QObject>

//project includes
#include "Residue.h"
#include "MDBox.h"
#include "AtomicInfo.h"

typedef QMap<uint,Residue> ResidueMap;
typedef QVector<Particle*> ParticlePointers;
typedef QMap<ResidueType,uint> LipidNumberMap;

class System {

	QString name_;
	ResidueMap residues_;
	ParticlePointers particles_;

	//For tracking the "current particle"
	ResidueMap::iterator res_itr;
	int cur_particle; //just an index for the QVector

	//For iterator returns to avoid null
	Residue dummy_res;

	//A system is defined by a volume of space and the stuff in that volume
	MDBox box;

public:

	System();
	System(const QString&);
	virtual ~System();

	const QString& name() const;
	void setName(const QString&);

	uint numberOfParticles() const;

	Residue& operator[] (uint);
	const Residue& operator[] (uint) const;
	uint numberOfResidues() const;

	//also use iterators, faster and easier in many situations
	Residue& firstResidue();
	Residue& nextResidue();
	bool atEnd() const;

	Residue& addResidue(const QString&,uint,ResidueType lipid = UNKNOWN_RESIDUE);
	void addParticle(Particle*);
	bool isEmpty() const;

	const MDBox& systemBox() const;
	void setBoxDimensions(double,double,double);
	void setBoxDimensions(const Vector3D&,const Vector3D&,const Vector3D&);
	const Matrix3D& boxDimensions() const;

	void clearSystem();

	void updateNextCoord(const Vector3D&);
	void updateMDBox(const Matrix3D&);
};

#endif /* SYSTEM_H_ */
