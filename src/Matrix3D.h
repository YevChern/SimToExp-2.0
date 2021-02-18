/*
 * Matrix3D.h
 *
 *  Created on: Apr 28, 2013
 *      Author: bholland
 */

#ifndef MATRIX3D_H_
#define MATRIX3D_H_

#include "Vector3D.h"

#include <QtCore/QPair>

typedef unsigned int uint;
typedef QPair<uint,uint> Element;

const uint NUM_MATRIX_ELEMENTS = 9;
const uint NUM_ROWS = 3;
const uint NUM_COLUMNS = 3;

class Matrix3D {

	double data[NUM_MATRIX_ELEMENTS];
	double err_out;

public:

	Matrix3D();
	Matrix3D(double);
	Matrix3D(double, double, double);
	Matrix3D(const Vector3D&, const Vector3D&, const Vector3D&);
	virtual ~Matrix3D();

	void zero();

	double& element(uint,uint);
	double& element(const Element&);
	const double& element(uint,uint) const;
	const double& element(const Element&) const;

	Vector3D diagonal() const;
	Vector3D column(uint) const;
	Vector3D row(uint) const;
	Matrix3D transpose() const;
	Matrix3D inverse() const;

	double determinant() const;
	double trace() const;

	Matrix3D abs() const;
	double sum() const;

	bool operator== (const Matrix3D&) const;
	bool operator!= (const Matrix3D&) const;

	Matrix3D operator+ (const Matrix3D&) const;
	Matrix3D operator- (const Matrix3D&) const;
	Matrix3D operator* (const Matrix3D&) const;
	Matrix3D operator* (double) const;
	Vector3D operator* (const Vector3D&) const;

	Matrix3D& operator+= (const Matrix3D&);
	Matrix3D& operator-= (const Matrix3D&);
	Matrix3D& operator*= (const Matrix3D&);
	Matrix3D& operator*= (double);

	//useful static members that return specific matrices
	static Matrix3D identityMatrix();
	static Matrix3D diagonalMatrix(double);
	static Matrix3D outerProduct(const Vector3D&,const Vector3D&);
};

#endif /* MATRIX3D_H_ */
