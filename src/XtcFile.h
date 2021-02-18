/*
 * XtcFile.h
 *
 *  Created on: May 3, 2013
 *      Author: bholland
 *
 *  Represents the xtc file and all of the information that defines it. Inherits from QFile
 *  and adds the extra functionality required to deal with an XTC file specifically.
 */

#ifndef XTCFILE_H_
#define XTCFILE_H_

#include <QtCore/QString>
#include <QtCore/QFile>
#include <QtCore/QDataStream>
#include <QtCore/QVector>

#include "XtcFile.h"
#include "System.h"
#include "Matrix3D.h"


typedef unsigned char uchar;

class XtcFile : public QFile {

	Q_OBJECT;

	QDataStream stream;
	System* system;

	//data from xtc header - need to be very specific with the 'int' size
	int num_atoms;
	float _precision;

	int first_step;			// Number of the first step
	float start_time;       // Time of first step
	double time_step;
	bool time_step_calculated;

	int cur_step;
	float cur_time;

	//container for magic integers
	QVector<int> magic_ints;

	int coord_block_size; //size of coordinate block in bytes

	//Buffer containers for reading in from the file
	QVector<char> file_data;
	QVector<int> int_buffer;

	//magic int references - in general should be constants, so left as uppercase
	int FIRSTIDX;
	int LASTIDX; //will get reset in the constructor to size of magic int container


public:

	XtcFile();
	XtcFile(const QString&, System*);

	//getters and setters
	uint numberOfAtoms() const;
	float precision() const;
	uint firstStep() const;
	float startTime() const;
	double timeStep() const;
	uint currentStep() const;
	float currentTime() const;

	void setSystem(System*);

	const QVector<int>& magicInts() const;

	bool openXtc();
	bool readNextFrame();
	bool readHeader();
	void resetToStart();

	void currentMDBox(Matrix3D);

private:

	void errorStream(const QString&) const;
	uint readCompressedCoords();

	uint bitsizeQInt32(int) const;
	uint bitsizeQInt32(const QVector<int>&) const;
	int receiveBits(quint32, quint32&, quint32&, quint32&);
	void receiveInts(quint32, quint32,QVector<int>&,quint32&,quint32&,quint32&);
	QVector<char>& readBytes(int);

	void initializeMembers();
	void swapValues(int&,int&);
	void setCoords(double,double,double);
};

#endif /* XTCFILE_H_ */
