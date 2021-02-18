/*
 * ExpReader.h
 *
 *  Created on: Nov 15, 2012
 *      Author: bholland
 *
 * A lot of the Manager functions are simply passed through to the reader which is somewhat superfluous, but this provides
 * a nice logical break from the Qt environment since no connections are made directly to the reader,
 * and are all handled by the ExpManager.
 */

#ifndef EXPREADER_H_
#define EXPREADER_H_

#include <QtCore/QString>
#include <QtCore/QList>
#include <QtCore/QVector>
#include <QtCore/QPointF>

#include "ExpData.h"
#include "FileCode.h"
#include "ScatteringType.h"

typedef unsigned int uint;
typedef QList<ExpData*> ExpDataPointers;
typedef QPair<QString,QVector<QPointF> > QStringPointsPair;
typedef QList<QPair<QString,QVector<QPointF> > > QStringPointsPairs;

class ExpManager; //for friendship

class ExpReader {

	//for loading in experimental form factors
	ExpDataPointers xray_data;
	ExpDataPointers neutron_data;

	//for loading in experimentally derived densities
	QStringPointsPairs elec_densities;
	QStringPointsPairs neut_densities;

public:

	friend class ExpManager;

	ExpReader();
	virtual ~ExpReader();

	//getters for objects in the containers
	ExpData* getExpDataObject(uint,ScatteringType);
	const ExpData* getExpDataObject(uint,ScatteringType) const;
	QStringPointsPair& getElectronDensity(uint);
	const QStringPointsPair& getElectronDensity(uint) const;
	QStringPointsPair& getNeutronSLDensity(uint);
	const QStringPointsPair& getNeutronSLDensity(uint) const;

	FileCode readExpFile(const QString&,ScatteringType);
	FileCode readExpDensityFile(const QString&,ScatteringType);

	bool formFactorEmpty(ScatteringType) const;

	uint numberDataSets(ScatteringType) const;
	uint numberElectronDensitySets() const;
	uint numberNeutronSLDensitySets() const;

	void removeData(int,ScatteringType);
	void removeElectronDensityData();
	void removeNeutronSLDensityData();
};

#endif /* EXPREADER_H_ */
