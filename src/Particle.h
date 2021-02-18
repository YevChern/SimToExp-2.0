/*
 * Particle.h
 *
 *  Created on: Feb 25, 2013
 *      Author: bholland
 *
 * Defines a atomistic or coarse grain particle, also storing its extended coordinates over
 * a trajectory.
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "Vector3D.h"
#include "AtomType.h"

#include <QtCore/QString>
#include <QtCore/QVector>

typedef unsigned int uint;

class Particle {

	QString coord_file_name;
	QString recognized_name;
	AtomType type;
	uint id_;

	Vector3D coord_;

public:

	Particle();
	Particle(const QString&, uint);
	virtual ~Particle();

	//getters
	const QString& coordFileName() const; //the name "as-is" in the coordinate file
	const QString& recognizedName() const; //the name for the atom type recognized by SIMtoEXP
	AtomType atomType() const;
	uint id() const;
	const Vector3D& coord() const;

	//setters
	void setCoordFileName(const QString&);
	void setRecognizedName(const QString&);
	void setAtomType(AtomType);
	void setID(uint);

	void setCoord(double,double,double);
	void setCoord(const Vector3D&);

	//eoperators
	bool operator< (const Particle&) const;
	bool operator== (const Particle&) const;
	bool operator!= (const Particle&) const;
};

#endif /* PARTICLE_H_ */
