
#include "LinearSolver.h"

#include <QtCore/QDebug>

#include <iostream>
#include <cmath>

using namespace std;

const double MIN_MATRIX_VALUE = 1E-30; //anything too small is likely a zero

namespace lsolver {

/* For solving a system of linear equations of type Ax = b (for x) using Gaussian elimination,  where 'A' is a matrix,
 * 'x' and 'b' are column vectors. This function uses partial pivoting to minimize errors introduced by finite precision.
 * Returns a vector of the solutions to the linear equations, which is empty if the matrix 'A' is singular. Also,
 * instead of crashing, if the matrix and vector sizes do not match up, a vector of {-1} is returned as an error
 * message.
 *
 * This follows Algorithms 2.4 (pg. 73) and 2.2 (pg. 66) of
 * "Scientific Computing: An Introductory Surver", 2nd Ed, Michael T. Heath
 */
Doubles gauss_elim(Equations& A, Doubles& b) {

	Doubles x;

	//first need to make sure the matrix column/row size = QVector length of b; this also guarantees a square matrix
	if(b.empty() or A.empty() or not (A.size() == b.size()) or not (A[0].size() == b.size())) {

		x.push_back(-1);
		qDebug() << "A size: " << A.size();
		qDebug() << "A[0] size: " << A[0].size();
		qDebug() << "b size: " << b.size();

		return x;
	}

	//if things are okay, start by performing elimination to create triangular matrix
	//outer loop is over columns, skip last column though
	for(int k = 0; k < A.size()-1; k++) {

		//search for largest pivot
		double pivot = A[k][k];
		int pivot_row = k;

		for(int i = k+1; i < A.size(); i++) {

			if(fabs(A[i][k]) > fabs(pivot)) {
				pivot = A[i][k];
				pivot_row = i;
			}
		}

		//if pivot has changed, swap rows
		if(pivot_row != k) {

			//change row in matrix
			QVector<double> temp = A[k];
			A[k] = A[pivot_row];
			A[pivot_row] = temp;

			//change row in vector
			double tmp = b[k];
			b[k] = b[pivot_row];
			b[pivot_row] = tmp;
		}

		//if the pivot is still zero, we have ourselves a singular matrix, nothing more to do!
		if(fabs(A[k][k]) <= MIN_MATRIX_VALUE)
			return x;

		//with matrix set up properly, compute multipliers
		QVector<double> mult;
		for(int i = k+1; i < A.size(); i++)
			mult.push_back(A[i][k] / A[k][k]);

		//now apply the transformation, but only to the portion of the matrix affected by it
		//outer loop is columns, inner is rows
		uint m;
		for(int j = k+1; j < A.size(); j++) {

			m = 0; //for multipliers
			for(int i = k+1; i < A.size(); i++) {

				A[i][j] -= mult[m] * A[k][j];
				m++;
			}
		}

		//apply transformation to 'b' as well
		m = 0;
		for(int i = k+1; i < b.size(); i++) {

			b[i] -= mult[m] * b[k];
			m++;
		}
	}

	//Matrix should now be transformed into an upper triangular - all values below the diagonal should
	//be zero, but are simply ignored instead.  Now use back-substitution to solve system
	x = Doubles (A.size()); //initialize with zeros
	for(int j = A.size()-1; j >= 0; j--) {

		//calculate solution at x[j]
		x[j] = b[j] / A[j][j];
		//qDebug() << "x: " << x[j];

		//update the vector for next row
		for(int i = 0; i < j; i++)
			b[i] -= A[i][j] * x[j];
	}

	//and we're done
	return x;
}

//A test system with solution: x1 = 3, x2 = -2, x3 = 4;
TestSystem testSystem() {

	Equations A(3);
	Doubles b;

	A[0].push_back(2);
	A[0].push_back(3);
	A[0].push_back(1);
	A[1].push_back(1);
	A[1].push_back(-1);
	A[1].push_back(1);
	A[2].push_back(5);
	A[2].push_back(2);
	A[2].push_back(2);

	b.push_back(4);
	b.push_back(9);
	b.push_back(19);

	return TestSystem(A,b);
}
}
