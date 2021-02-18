/*
 * OutputWriter.h
 *
 *  Created on: Jan 13, 2014
 *      Author: bholland
 */

#ifndef OUTPUTWRITER_H_
#define OUTPUTWRITER_H_

//C++ stdlib includes
#include <fstream>

//Qt includes
#include <QtCore/QString>
#include <QtCore/QList>

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QWidget>
#else
    #include <QtGui/QWidget>
#endif

//SIMtoEXP includes
#include "NormalDensity.h"
#include "Component.h"

typedef QList<NormalDensity> NormalDensities;
typedef QList<Component> Components;

enum OutputType {NUMBER_DENSITY, ELECTRON_DENSITY, NEUTRON_DENSITY, XRAY_STRUCTURE_FACTOR, NEUTRON_STRUCTURE_FACTOR, VOLUMETRIC_FIT};

/*
 * Used for all writing/saving to file.  Might want to make separate classes later, but for now
 * just do it all here.
 */
class OutputWriter {

	//only use a single output stream since all files will be opened and closed upon writing
	std::ofstream out_stream;

	//the parent window the writer is working for
	QWidget* parent;

	//header for all files
	QString header;

public:

	OutputWriter(QWidget* parent = 0);
	virtual ~OutputWriter();

	void setParent(QWidget*);

	void writeNumberDensities(const NormalDensities&, uint, double);
	void writeComponentNumberDensities(const Components&);

	void writeElectronDensities(const NormalDensities&);
	void writeComponentElectronDensities(const Components&);

	void writeNeutronDensities(const NormalDensities&);
	void writeComponentNeutronDensities(const Components&);

	void writeXrayStructureFactors(const QVector<QPointF>&);
	void writeNeutronStructureFactors(const QVector<QPointF>&);
	void writeVolumetricFit(const QVector<Points>&);

private:

	QString getFilePath(OutputType) const;
	void setHeader();
	QString addSpaces(const QString&,uint) const;
	void writeDensityData(const QVector<Points>&, const NormalDensities&);
	void writeComponentData(const QVector<Points>&, const Components&);
	void writeData(const QVector<Points>&, double);
	void writeStructureFactorData(const Points&);

	void writeComponents(OutputType, const Components&);
};

#endif /* OUTPUTWRITER_H_ */
