/*
 * Residue.h
 *
 *  Created on: Feb 27, 2013
 *      Author: bholland
 */

#ifndef RESIDUE_H_
#define RESIDUE_H_

#include "Particle.h"
#include "ParticleLabel.h"

#include <QtCore/QPair>
#include <QtCore/QMap>

typedef QMap<uint,Particle> ParticleMap;

class Residue {

	QString name_;
	uint res_id;
	ParticleMap particles_;
	ParticleMap::iterator cur_particle;

	//dummy particle for returning at end of iterator
	Particle dummy_p;

public:

	Residue();
	Residue(const QString&,uint);
	virtual ~Residue();

	//getters
	const QString& name() const;
	uint resId() const;
	const Particle& particle(uint) const;
	QPair<ParticleLabel,AtomType> particleLabel(uint) const;
	QPair<ParticleLabel,AtomType> elementalLabel(uint) const;

	uint numberOfParticles() const;
	bool isEmpty() const;

	//setters
	void setName(const QString&);
	void setResID(uint);
	Particle& addParticle(const QString&, uint);
	void addParticle(const Particle&);

	//iteration
	Particle& firstParticle();
	Particle& nextParticle();
	bool atEnd() const;

	//operators
	Particle& operator[] (uint);
	const Particle& operator[] (uint) const;
	bool operator== (const Residue&) const;
	bool operator!= (const Residue&) const;
};

#endif /* RESIDUE_H_ */
