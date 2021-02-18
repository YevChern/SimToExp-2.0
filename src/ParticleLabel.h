/*
 * ParticleLabel.h
 *
 *  Created on: Mar 13, 2013
 *      Author: bholland
 *
 *  Simple class used to define (or label) a specific atom within a molecule.
 *  Each of the same molecule type will have the same labels, and thus the intensities
 *  can be stored properly using this class as a map key.
 */

#ifndef PARTICLELABEL_H_
#define PARTICLELABEL_H_

#include <QtCore/QString>

#include "AtomType.h"

class ParticleLabel {

	QString res_name;
	QString particle_name;

public:

	ParticleLabel();
	ParticleLabel(const QString&,const QString&);
	virtual ~ParticleLabel();

	const QString& residueName() const;
	const QString& particleName() const;

	void setResidueName(const QString&);
	void setParticleName(const QString&);

	bool operator<(ParticleLabel) const;
	bool operator==(ParticleLabel) const;
};

#endif /* PARTICLELABEL_H_ */
