/*
 * MDBox.h
 *
 *  Created on: Apr 27, 2013
 *      Author: bholland
 *
 * Represents the actual container that the system go into.  The container can be any shape one can create using
 * three 3D vectors , i.e. if the matrix is diagonal then the container is a box, otherwise it will be a different shape
 * that depends on the vectors.  Usually with periodic boundary conditions, only a few different shapes are allowable, e.g.
 * box, rhombus, dodecahedron, etc.
 */

#ifndef MDBOX_H_
#define MDBOX_H_

#include "Matrix3D.h"
#include "Axis.h"

class MDBox {

	Matrix3D real;	// real space
	Matrix3D reciprocal;	// reciprocal space vectors

public:

	MDBox();
	MDBox(const Matrix3D&);
	MDBox(const Vector3D&, const Vector3D&, const Vector3D&);
	MDBox(double, double, double);
	MDBox(double);
	virtual ~MDBox();

	const Matrix3D& realSpace() const;
	const Matrix3D& reciprocalSpace() const;

	//real space row vectors that together define the spatial tensor
	Vector3D a() const;
	Vector3D b() const;
	Vector3D c() const;

	//inverted space row vectors
	Vector3D u() const;
	Vector3D v() const;
	Vector3D w() const;

	//real space diagonals - for a box these define the space
	double x() const;
	double y() const;
	double z() const;

	double diagThickness(Axis) const;

	bool isBox() const;
	bool isCube() const;

	Vector3D realDiagonal() const;
	Vector3D reciprocalDiagonal() const;

	void setBoxDimensions(double,double,double);
	void setBoxDimensions(const Vector3D&, const Vector3D&, const Vector3D&);
	void setBoxDimensions(const Matrix3D&);

	void scale(double);

	double volume() const;

	bool operator== (const MDBox&) const;
	bool operator!= (const MDBox&) const;

private:

	// Initializer
	void init(const Vector3D&, const Vector3D&, const Vector3D&);

};

#endif /* MDBOX_H_ */
