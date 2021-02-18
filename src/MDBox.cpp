/*
 * MDBox.cpp
 *
 *  Created on: Apr 27, 2013
 *      Author: bholland
 */

#include "MDBox.h"
#include "MathUtils.h"
#include <cmath>

typedef QPair<uint,uint> Element;

const uint REAL_A_ROW = 0;
const uint REAL_B_ROW = 1;
const uint REAL_C_ROW = 2;
const uint RECIP_U_ROW = 0;
const uint RECIP_V_ROW = 1;
const uint RECIP_W_ROW = 2;

const Element X_INDEX = Element(0,0);
const Element Y_INDEX = Element(1,1);
const Element Z_INDEX = Element(2,2);

//Default constructor - just initialize everything to zero
MDBox::MDBox() {

	Vector3D a, b, c; //all initialize to zero
	init(a, b, c);
}

//Constructor that uses a Matrix object to set the 'real' space components
MDBox::MDBox(const Matrix3D& box) {

	real = box;
	reciprocal = box.inverse();
}

//Constructor that sets the real space vectors to the given values
MDBox::MDBox(const Vector3D& a, const Vector3D& b, const Vector3D& c) {

	init(a, b, c);
}

//Constructor that sets the given values as the diagonal of the real space vectors, i.e. creates a box
MDBox::MDBox(double x, double y, double z) {

	Vector3D a(x, 0, 0), b(0, y, 0), c(0, 0, z);
	init(a, b, c);
}

//Constructor that sets the given value as the diagonal of the real space vectors, i.e. creates a cubic box
MDBox::MDBox(double x) {

	Vector3D a(x, 0, 0), b(0, x, 0), c(0, 0, x);
	init(a, b, c);
}

//Destructor
MDBox::~MDBox() {}

//Getters
const Matrix3D& MDBox::realSpace() const {return real;}
const Matrix3D& MDBox::reciprocalSpace() const {return reciprocal;}

Vector3D MDBox::a() const {return real.row(REAL_A_ROW);}
Vector3D MDBox::b() const {return real.row(REAL_B_ROW);}
Vector3D MDBox::c() const {return real.row(REAL_C_ROW);}

Vector3D MDBox::u() const {return reciprocal.row(RECIP_U_ROW);}
Vector3D MDBox::v() const {return reciprocal.row(RECIP_V_ROW);}
Vector3D MDBox::w() const {return reciprocal.row(RECIP_W_ROW);}

double MDBox::x() const {return real.element(X_INDEX);}
double MDBox::y() const {return real.element(Y_INDEX);}
double MDBox::z() const {return real.element(Z_INDEX);}

//just a getter like x(), y() and z(), but determines which to return from the axis
double MDBox::diagThickness (Axis axis) const {

	switch(axis) {

	case X_AXIS:
		return real.element(X_INDEX);

	case Y_AXIS:
		return real.element(Y_INDEX);

	case Z_AXIS:
		return real.element(Z_INDEX);

	default:
		return 0; //should never happen
	}
}

//set the diagonal of the box dimensions to the given values
void MDBox::setBoxDimensions(double x, double y, double z) {

	Vector3D a(x, 0, 0), b(0, y, 0), c(0, 0, z);
	init(a, b, c);
}

//Set the box to the given vectors
void MDBox::setBoxDimensions(const Vector3D& a, const Vector3D& b, const Vector3D& c) {

	init(a, b, c);
}

//Set the box to the given spatial tensor
void MDBox::setBoxDimensions(const Matrix3D& box) {

	real = box;
}

//Scales every value in the box
void MDBox::scale(double factor) {

	real *= factor;
	reciprocal = real.inverse();
}

//calculates and returns the volume of the box
double MDBox::volume() const {

	return fabs(a().cross(b()).dot(c()));
}

//Checks to see if the container has all right angles, i.e. is the matrix diagonal
bool MDBox::isBox() const {

	//go through non-diagonals and check for zero
	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			if(i != j and real.element(i,j) != 0)
				return false;
		}
	}

	return true;
}

//Checks to see if the container has all right angles and if all sides are the same length
bool MDBox::isCube() const {

	//first check to see if a box
	if(not isBox())
		return false;

	//then check diagonal elements
	if(real.element(X_INDEX) == real.element(Y_INDEX) and
	   real.element(X_INDEX) == real.element(Z_INDEX) and
	   real.element(Y_INDEX) == real.element(Z_INDEX))
		return true;

	return false;
}

//Return the diagonal of the box matrix
Vector3D MDBox::realDiagonal() const {

	return real.diagonal();
}

//Return the inverse diagonal (i.e. inverted element by element) of the matrix
Vector3D MDBox::reciprocalDiagonal() const {

	return reciprocal.diagonal();
}

// Initializer
void MDBox::init(const Vector3D& a, const Vector3D& b, const Vector3D& c) {

	real = Matrix3D(a, b, c);
	reciprocal = real.inverse();
}

/*
 * Equality operator - equal iff real space matrices are equal -> NO tolerance accepted
 */
bool MDBox::operator ==(const MDBox& other) const {

	if(other.realSpace() == real)
		return true;
	return false;
}

/*
 * Inequality operator - see above
 */
bool MDBox::operator !=(const MDBox& other) const {

	return not (*this == other);
}
