/*
 * ExpParser.h
 *
 *  Created on: Oct 19, 2012
 *      Author: bholland
 * 
 *  Read the experimental data files and stores the data in 
 */

#ifndef EXPPARSER_H_
#define EXPPARSER_H_

#include "ExpReader.h"

typedef QVector<QPointF> Points;
typedef QPair<QString,QVector<QPointF> > QStringPointsPair;

class ExpManager : public QWidget {

	Q_OBJECT

	ExpReader reader;

	bool ed_data_loaded;
	bool nsld_data_loaded;

	//bounding rectangles
	QRectF ed_rect;
	QRectF nsld_rect;
	QRectF xrayFF_rect;
	QRectF neutFF_rect;

public:

	ExpManager(QWidget* parent = 0);
	virtual ~ExpManager();

	FileCode readExpFile(const QString&,ScatteringType);
	FileCode readExpDensityFile(const QString&,ScatteringType);

	ExpData& getXrayDataObject(uint);
	const ExpData& getXrayDataObject(uint) const;
	ExpData& getNeutronDataObject(uint);
	const ExpData& getNeutronDataObject(uint) const;

	Points& getXrayData(uint);
	const Points& getXrayData(uint) const;
	Points& getNeutronData(uint);
	const Points& getNeutronData(uint) const;

	QStringPointsPair& getElectronDensityData(uint);
	const QStringPointsPair& getElectronDensityData(uint) const;
	QStringPointsPair& getNeutronSLDensityData(uint);
	const QStringPointsPair& getNeutronSLDensityData(uint) const;

	bool electronDensityLoaded() const;
	bool neutronSLDensityLoaded() const;
	bool xrayFormFactorEmpty() const;
	bool neutronFormFactorEmpty() const;

	const QRectF& electronDensityRect();
	const QRectF& neutronSLDensityRect();
	const QRectF& xrayFormFactorRect();
	const QRectF& neutronFormFactorRect();

	double getScalingFactor(uint,ScatteringType) const;
	double getChiSquared(uint,ScatteringType) const;

	uint getNumberXray() const;
	uint getNumberNeutron() const;
	uint getNumberElectronDensity() const;
	uint getNumberNeutronSLDensity() const;

	void removeXrayData(uint);
	void removeNeutronData(uint);

	void removeElectronDensityData();
	void removeNeutronSLDensityData();

	void scaleData(const Points&,const Points&); //for automatic scaling
	void scaleData(const Points&,ScatteringType); //automatic for single type
	void scaleXrayData(const Doubles&); //for providing scaling factors
	void scaleNeutronData(const Doubles&); //for providing scaling factors

public slots:

	void resetAllData();

signals:

	void xrayDataScaled();
	void neutronDataScaled();
};

#endif /* EXPPARSER_H_ */
