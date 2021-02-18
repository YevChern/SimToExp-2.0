/*
 * Vector3D.h
 *
 *  Created on: Apr 24, 2013
 *      Author: Dmitri Rozmanov and B. Holland
 *      
 * Cleaned up version of DR's 3D vector class
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

typedef unsigned int uint;

class Vector3D {

	double _x;
	double _y;
	double _z;

public:

	Vector3D();
	Vector3D(double, double, double);
	Vector3D(double);
	virtual ~Vector3D();

	double x() const;
	double y() const;
	double z() const;

	void setX(double);
	void setY(double);
	void setZ(double);

	void set(const Vector3D&);
	void set(double, double, double);
	void zero();

	double dot(const Vector3D&) const;
	Vector3D cross(const Vector3D&) const;

	double norm() const; //i.e. the magnitude or length
	double norm2() const;
	Vector3D abs() const;
	Vector3D inv() const;

	void add(const Vector3D&);
	void sub(const Vector3D&);
	void mult(const Vector3D&);
	void scale(double);

	void round();
	void ceiling();
	void floor();

	Vector3D normalize() const;
	double scalarProjection(const Vector3D&) const;
	Vector3D vectorProjection(const Vector3D&) const;

	double sum() const;

	Vector3D operator+ (const Vector3D&) const;
	Vector3D operator- (const Vector3D&) const;
	Vector3D operator* (const Vector3D&) const;
	Vector3D operator* (double) const;

	Vector3D& operator+= (const Vector3D&);
	Vector3D& operator-= (const Vector3D&);
	Vector3D& operator*= (double);
	Vector3D& operator*= (const Vector3D&);
	Vector3D& operator/= (double);
	Vector3D& operator/= (const Vector3D&);

	bool operator== (const Vector3D&) const;
	bool operator!= (const Vector3D&) const;
	bool isOrigin() const;

	double& operator[] (uint);
	double operator[] (uint) const;

	double angle(const Vector3D&) const;
};

#endif /* VECTOR3D_H_ */
