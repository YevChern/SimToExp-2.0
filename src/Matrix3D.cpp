/*
 * Matrix3D.cpp
 *
 *  Created on: Apr 28, 2013
 *      Author: bholland
 */

#include <iostream>
#include <cfloat>
#include <cmath>

#include "Matrix3D.h"

using namespace std;

//index constants that map position in the matrix to the data array
const uint elem_11 = 0;
const uint elem_12 = 1;
const uint elem_13 = 2;
const uint elem_21 = 3;
const uint elem_22 = 4;
const uint elem_23 = 5;
const uint elem_31 = 6;
const uint elem_32 = 7;
const uint elem_33 = 8;

//Default constructor
Matrix3D::Matrix3D() {

	for(uint i = 0; i < NUM_MATRIX_ELEMENTS; i++)
		data[i] = 0.0;

	err_out = -DBL_MAX;
}

//Constructor that sets all elements to the same value
Matrix3D::Matrix3D(double value) {

	for(uint i = 0; i < NUM_MATRIX_ELEMENTS; i++)
		data[i] = value;

	err_out = -DBL_MAX;
}

//Constructor that sets the diagonal, rest to zero
Matrix3D::Matrix3D(double a11, double a22, double a33) {

	*this = Matrix3D(0.0);
	data[elem_11] = a11;
	data[elem_22] = a22;
	data[elem_33] = a33;

	err_out = -DBL_MAX;
}

//Constructor that sets the rows of the matrix to the given vectors
Matrix3D::Matrix3D(const Vector3D& r1, const Vector3D& r2, const Vector3D& r3) {

	data[elem_11] = r1.x();
	data[elem_12] = r1.y();
	data[elem_13] = r1.z();
	data[elem_21] = r2.x();
	data[elem_22] = r2.y();
	data[elem_23] = r2.z();
	data[elem_31] = r3.x();
	data[elem_32] = r3.y();
	data[elem_33] = r3.z();

	err_out = -DBL_MAX;
}


//Destructor
Matrix3D::~Matrix3D() {}

//Sets all elements to zero
void Matrix3D::zero() {

	for(uint i = 0; i < NUM_MATRIX_ELEMENTS; i++)
		data[i] = 0.0;
}

//Returns a reference to the element at (i,j), so can be manipulated
double& Matrix3D::element(uint i, uint j) {

	//check to make sure i and j are within bound
	//if not -> output error message and return -ve DBL_MAX
	if(i < NUM_ROWS and j < NUM_COLUMNS)
		return data[i*3 + j];

	cerr << "One or more element indexes is out of bounds, must be < " << NUM_ROWS << endl;
	return err_out;
}

//Returns a reference to the element at (i,j), so can be manipulated
double& Matrix3D::element(const Element& index) {

	//check to make sure i and j are within bound
	//if not -> output error message and return -ve DBL_MAX
	if(index.first < NUM_ROWS and index.second < NUM_COLUMNS)
		return data[index.first*3 + index.second];

	cerr << "One or more element indexes is out of bounds, must be < " << NUM_ROWS << endl;
	return err_out;
}


//Returns the value at (i,j)
const double& Matrix3D::element(uint i, uint j) const {

	//check to make sure i and j are within bound
	//if not -> output error message and return -ve DBL_MAX
	if(i < NUM_ROWS and j < NUM_COLUMNS)
		return data[i*3 + j];

	cerr << "One or more element indexes is out of bounds, must be < " << NUM_ROWS << endl;
	return err_out;
}

//Returns the value at (i,j)
const double& Matrix3D::element(const Element& index) const {

	//check to make sure i and j are within bound
	//if not -> output error message and return -ve DBL_MAX
	if(index.first < NUM_ROWS and index.second < NUM_COLUMNS)
		return data[index.first*3 + index.second];

	cerr << "One or more element indexes is out of bounds, must be < " << NUM_ROWS << endl;
	return err_out;
}

//Return a vector representing the diagonal of the matrix, i.e. (a11, a22, a33)
Vector3D Matrix3D::diagonal() const {

	return Vector3D(data[elem_11],data[elem_22],data[elem_33]);
}

//Return the row vector at 'i'
Vector3D Matrix3D::row(uint i) const {

	double a1 = element(i, 0);
	double a2 = element(i, 1);
	double a3 = element(i, 2);

	return Vector3D(a1, a2, a3);
}

//Return the column vector at 'j'
Vector3D Matrix3D::column(uint j) const {

	double a1 = element(0, j);
	double a2 = element(1, j);
	double a3 = element(2, j);

	return Vector3D(a1, a2, a3);
}

//Returns the transpose of the current matrix
Matrix3D Matrix3D::transpose() const {

	Matrix3D trans;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			trans.element(i, j) = element(j, i);
		}
	}

	return trans;
}

/*
 * Returns the mathematical inverse matrix of the current object, where A*A^-1 = I, i.e. the
 * matrix times it's inverse gives the identity matrix.  General solutions will use some sort of
 * decomposition, but for a 3x3, still easier to use the analytical solution
 */
Matrix3D Matrix3D::inverse() const {

	Matrix3D inverse;

	inverse.element(0,0) = data[elem_22]*data[elem_33] - data[elem_23]*data[elem_32];
	inverse.element(0,1) = data[elem_13]*data[elem_32] - data[elem_12]*data[elem_33];
	inverse.element(0,2) = data[elem_12]*data[elem_23] - data[elem_13]*data[elem_22];
	inverse.element(1,0) = data[elem_23]*data[elem_31] - data[elem_21]*data[elem_33];
	inverse.element(1,1) = data[elem_11]*data[elem_33] - data[elem_13]*data[elem_31];
	inverse.element(1,2) = data[elem_13]*data[elem_21] - data[elem_11]*data[elem_23];
	inverse.element(2,0) = data[elem_21]*data[elem_32] - data[elem_22]*data[elem_31];
	inverse.element(2,1) = data[elem_12]*data[elem_31] - data[elem_11]*data[elem_32];
	inverse.element(2,2) = data[elem_11]*data[elem_22] - data[elem_12]*data[elem_21];

	double inv_det = 1/determinant();
	return inverse * inv_det;
}

//Returns the determinant of the current matrix
double Matrix3D::determinant() const {

	double det_11 = data[elem_11] * (data[elem_22]*data[elem_33] - data[elem_23]*data[elem_32]);
	double det_12 = data[elem_12] * (data[elem_21]*data[elem_33] - data[elem_23]*data[elem_31]);
	double det_13 = data[elem_13] * (data[elem_21]*data[elem_32] - data[elem_22]*data[elem_31]);

	return det_11 - det_12 + det_13;
}

//Returns the trace of the current matrix, i.e. the sum of the diagonal
//still returns a value for non-square matrix (that doesn't make sense), but this should never happen for a 3D matrix (tensor).
double Matrix3D::trace() const {

	double sum = 0;
	for(uint i = 0; i < NUM_COLUMNS; i++)
		sum += element(i,i);

	return sum;
}

//Returns a matrix with the absolute value of each element of the current matrix
Matrix3D Matrix3D::abs() const {

	Matrix3D abs;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			abs.element(i, j) = fabs(element(i, j));
		}
	}

	return abs;
}

//Return the sum of all elements
double Matrix3D::sum() const {

	double sum = 0;
	for(uint i = 0; i < NUM_MATRIX_ELEMENTS; i++)
		sum += data[i];

	return sum;
}

//Equality operators- only equal if all elements are
bool Matrix3D::operator== (const Matrix3D& other) const {

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			if(other.element(i, j) != element(i, j))
				return false;
		}
	}

	//if made it all way through, all equal
	return true;
}

bool Matrix3D::operator!= (const Matrix3D& other) const {

	return not ((*this) == other);
}

//Operators for adding and subtracting the elements of the given matrix with the current one
//---------------------------------------------------------
Matrix3D Matrix3D::operator+ (const Matrix3D& other) const {

	Matrix3D sum;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			sum.element(i, j) = element(i, j) + other.element(i, j);
		}
	}

	return sum;
}

Matrix3D Matrix3D::operator- (const Matrix3D& other) const {

	Matrix3D diff;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			diff.element(i, j) = element(i, j) - other.element(i, j);
		}
	}

	return diff;
}

//Operators that multiplies the two matrices using standard matrix multiplication
Matrix3D Matrix3D::operator* (const Matrix3D& other) const {

	Matrix3D product;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			product.element(i, j) = element(i, 0) * other.element(0, j) +
									element(i, 1) * other.element(1, j) +
									element(i, 2) * other.element(2, j);
		}
	}

	return product;
}

//Operators that returns a matrix scaled by the given value
Matrix3D Matrix3D::operator* (double scalar) const {

	Matrix3D product;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			product.element(i,j) = element(i,j) * scalar;
		}
	}

	return product;
}

//Operator that returns the vector result of multiplying the matrix by a column vector
Vector3D Matrix3D::operator* (const Vector3D& vector) const {

	Vector3D product;

	for(uint i = 0; i < NUM_ROWS; i++) {

			product[i] = element(i, 0) * vector.x() +
					     element(i, 1) * vector.y() +
					     element(i, 2) * vector.z();
	}

	return product;
}

//----------------------------------------------------------

//Adds the elements of the given matrix to the elements of the current one
Matrix3D& Matrix3D::operator+= (const Matrix3D& other) {

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			element(i,j) += other.element(i,j);
		}
	}

	return (*this);
}

//Subtracts the elements of the given matrix from the elements of the current one
Matrix3D& Matrix3D::operator-= (const Matrix3D& other) {

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			element(i,j) -= other.element(i,j);
		}
	}

	return (*this);
}

//Performs a matrix multiplication and sets the current matrix to the result
Matrix3D& Matrix3D::operator*= (const Matrix3D& other) {

	Matrix3D product = (*this) * other;
	*this = product;

	return (*this);
}

//Scales the current matrix by the given value
Matrix3D& Matrix3D::operator*= (double scalar) {

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			element(i,j) *= scalar;
		}
	}

	return (*this);
}

//Some useful static functions
Matrix3D Matrix3D::identityMatrix() {

	return Matrix3D(1,1,1);
}

Matrix3D Matrix3D::diagonalMatrix(double value) {

	return Matrix3D(value, value, value);
}

Matrix3D Matrix3D::outerProduct(const Vector3D& v1, const Vector3D& v2) {

	Matrix3D result;

	for(uint i = 0; i < NUM_ROWS; i++) {
		for(uint j = 0; j < NUM_COLUMNS; j++) {

			result.element(i,j) = v1[i] * v2[j];
		}
	}

	return result;
}
