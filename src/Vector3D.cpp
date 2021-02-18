#include "Vector3D.h"

#include <cfloat>
#include <string>
#include <cmath>
#include <iostream>

using namespace std;

const uint output_precision = 17;

//Default constructor - sets vector to origin
Vector3D::Vector3D() {

	_x = 0.0;
	_y = 0.0;
	_z = 0.0;
}

//Constructor that sets the vector to the given values
Vector3D::Vector3D(double x, double y, double z) {

	_x = x;
	_y = y;
	_z = z;
}

//Constructor that sets all components to the same value
Vector3D::Vector3D(double value) {

	_x = _y = _z = value;
}

//Destructor
Vector3D::~Vector3D() {};

//getters
double Vector3D::x() const {return _x;}
double Vector3D::y() const {return _y;}
double Vector3D::z() const {return _z;}

//setters
void Vector3D::setX(double x) {_x = x; }
void Vector3D::setY(double y) {_y = y; }
void Vector3D::setZ(double z) {_z = z; }

void Vector3D::set(const Vector3D& other) {

	_x = other.x();
	_y = other.y();
	_z = other.z();
}

void Vector3D::set(double x, double y, double z) {

	_x = x;
	_y = y;
	_z = z;
}

//set the vector to the origin
void Vector3D::zero() {

	_x = 0.0;
	_y = 0.0;
	_z = 0.0;
}

//calculates and return the dot/scalar product of 'this' vector with the given vector
double Vector3D::dot(const Vector3D& other) const {

	return _x*other.x() + _y*other.y() + _z*other.z();
}

//calculates and returns the cross/vector product of 'this' vector with the given vector
Vector3D Vector3D::cross(const Vector3D& other) const {

	return Vector3D(_y * other.z() - _z * other.y(),
					_z * other.x() - _x * other.z(),
					_x * other.y() - _y * other.x());
}

//returns the dot product of the vector with itself
double Vector3D::norm2() const {

	return dot(*this);
}

//returns the norm of the current vector
double Vector3D::norm() const {

	return std::sqrt(norm2());
}

//return a vector that is the absolute value of the current vector
Vector3D Vector3D::abs() const {

	return Vector3D(fabs(_x), fabs(_y), fabs(_z));
}

//return a vector that has the inverse of the current components
Vector3D Vector3D::inv() const {

	return Vector3D(1.0/_x, 1.0/_y, 1.0/_z);
}

//Adds the given vector to the current one
void Vector3D::add(const Vector3D& other) {

	(*this) += other;
}

//Subtracts the given vector from the current one
void Vector3D::sub(const Vector3D& other) {

	(*this) -= other;
}

//Multiplies the components of the given vector and the current one
void Vector3D::mult(const Vector3D& other) {

	(*this) *= other;
}

//Scales the current vector by the given value
void Vector3D::scale(double c) {

	_x *= c;
	_y *= c;
	_z *= c;
}

//Converts the components of the current vector to their respective ceilings
void Vector3D::ceiling() {

	_x = (int)((_x < 0) ? (_x) : (_x + 1.0));
	_y = (int)((_y < 0) ? (_y) : (_y + 1.0));
	_z = (int)((_z < 0) ? (_z) : (_z + 1.0));
}

//Converts the components of the current vector to their respective floors
void Vector3D::floor() {

	_x = (int)((_x < 0) ? (_x - 1.0) : (_x));
	_y = (int)((_y < 0) ? (_y - 1.0) : (_y));
	_z = (int)((_z < 0) ? (_z - 1.0) : (_z));
}

//Rounds the current vector to the nearest integer
void Vector3D::round() {

	_x = (int)((_x < 0) ? (_x - 0.5) : (_x + 0.5));
	_y = (int)((_y < 0) ? (_y - 0.5) : (_y + 0.5));
	_z = (int)((_z < 0) ? (_z - 0.5) : (_z + 0.5));
}

//Returns a the unit vector that corresponds to the current vector
Vector3D Vector3D::normalize() const {

	return (*this) * (1/norm());
}

//Returns the length of the projected vector of the current vector along the given axis
double Vector3D::scalarProjection(const Vector3D& axis) const {

	return dot(axis.normalize());
}

//Returns the projection of the current vector along the given axis
Vector3D Vector3D::vectorProjection(const Vector3D& axis) const {

	Vector3D axis_uvect = axis.normalize();
	return axis_uvect * dot(axis_uvect);
}

//Returns the sum of the component values
double Vector3D::sum() const {

	return _x + _y + _z;
}

//Operators for addition, subtraction and multiplication
Vector3D Vector3D::operator+ (const Vector3D& other) const {

	double x = _x + other.x();
	double y = _y + other.y();
	double z = _z + other.z();
	return Vector3D(x, y, z);
}

Vector3D Vector3D::operator- (const Vector3D& other) const {

	double x = _x - other.x();
	double y = _y - other.y();
	double z = _z - other.z();
	return Vector3D(x, y, z);
}

Vector3D Vector3D::operator* (const Vector3D& other) const {

	double x = _x * other.x();
	double y = _y * other.y();
	double z = _z * other.z();
	return Vector3D(x, y, z);
}

Vector3D Vector3D::operator* (double c) const {

	double x = _x * c;
	double y = _y * c;
	double z = _z * c;
	return Vector3D(x, y, z);
}

//Operators for adding, subtracting and multiplying with the current vector and setting the result in the current vector
Vector3D& Vector3D::operator+= (const Vector3D& other) {

	_x += other.x();
	_y += other.y();
	_z += other.z();
	return *this;
}

Vector3D& Vector3D::operator-= (const Vector3D& other) {

	_x -= other.x();
	_y -= other.y();
	_z -= other.z();
	return *this;
}

Vector3D& Vector3D::operator*= (double c) {

	_x *= c;
	_y *= c;
	_z *= c;
	return *this;
}

Vector3D& Vector3D::operator*= (const Vector3D& other) {

	_x *= other.x();
	_y *= other.y();
	_z *= other.z();
	return *this;
}

Vector3D& Vector3D::operator/= (double c) {

	scale(1/c);
	return *this;
}

Vector3D& Vector3D::operator/= (const Vector3D& other) {

	_x /= other.x();
	_y /= other.y();
	_z /= other.z();
	return *this;
}

//Equality operators - all components must be equal
bool Vector3D::operator== (const Vector3D& other) const {

	return (_x == other.x() and _y == other.y() and _z == other.z());
}

bool Vector3D::operator!= (const Vector3D& other) const {

	return not (*this == other);
}

//Array style operator; any index out of bounds will return an error message
double& Vector3D::operator[] (uint index) {

	switch(index) {

	case 0:
		return _x;

    case 1:
    	return _y;

    case 2:
		return _z;

    default:
    	cerr << "Vector3D index out of bounds! Returned _x anyway" << endl;
    	return _x;
	}
}

//Constant version of array style operator; any index out of bounds will return an error message
double Vector3D::operator[] (uint index) const {

	switch(index) {

	case 0:
		return _x;

    case 1:
    	return _y;

    case 2:
		return _z;

    default:
    	cerr << "Vector3D index out of bounds!" << endl;
    	return -DBL_MAX;
	}
}

//Boolean for origin, or all zeros
bool Vector3D::isOrigin() const {

	return (_x == 0.0 and _y == 0.0 and _z == 0.0);
}

//Returns the angle between the two vectors
double Vector3D::angle(const Vector3D& other) const {

	if(isOrigin() or other.isOrigin())
		return 0.0;

	double length = norm();
	double other_length = other.norm();

	double cosine = dot(other) / (length * other_length);

	return acos(cosine);
}
