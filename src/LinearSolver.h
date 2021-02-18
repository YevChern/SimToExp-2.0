/*
 * LinearSolver.h
 *
 * Created on: Oct, 2012
 * Author: Bryan Holland
 *
 *
 *
 * Description:
 *
 * Contains functions for solving various linear algebra equations / problems.
 *
 */

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

#include <QtCore/QVector>
#include <QtCore/QPair>

typedef unsigned int uint;
typedef QVector<double> Doubles;
typedef QVector<QVector<double> > Equations;
typedef QPair<QVector<QVector<double> >,QVector<double> > TestSystem;

namespace lsolver {

   Doubles gauss_elim(Equations&, Doubles&);
   TestSystem testSystem();
}

#endif /* LINEARSOLVER_H_ */
 
