/*
* FourierFormFactor.h
*
*  Created on: Oct 5, 2012
*      Author: bholland
*
* Function for determining the final form factor over a range of 'Q' (given by the margins)
*
* Replaces the 'FourierTransform' struct from previous version (file = fourier.h/cpp)
*
*/

#ifndef FOURIERFORMFACTOR_H_
#define FOURIERFORMFACTOR_H_ 

#include <QtCore/QVector>
#include <QtCore/QMap>
#include <QtCore/QPointF>
#include <QtCore/QString>
#include <QtCore/QMutex>

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QProgressDialog>
#else
    #include <QtGui/QProgressDialog>

#endif

#include "ScatteringType.h"
#include "XrayAtomicFormFactor.h"
#include "NormalDensity.h"

typedef QVector<double> Doubles;
typedef QVector<QPointF> Points;
typedef QMap<AtomType,XrayAtomicFormFactor> AtomTypeXrayMap;
typedef QList<NormalDensity> NormalDensities;
typedef QList<AtomType> AtomTypes;

class FourierFormFactor {
   
   //members
   Doubles x;  //positions along function in real space
   Doubles q; //positions in reciprical space (some redundancy with QPointF below, but easier this way?)
   double solvent_sd; //solvent scattering density (get from GUI)
   double q_min,q_max,q_step;
   bool symmetrized;
   
   ScatteringType ff_type;
   AtomicInfo* info; //useful to keep a pointer to the info

   AtomTypeXrayMap xray_atomic_FFs; //replaces array of AFF structs.
   Points form_factors; //final form factors
   
   QMutex ff_mutex; //for multithreading safety

   public:
      
      FourierFormFactor(ScatteringType,AtomicInfo&);
      virtual ~FourierFormFactor();
      
      //these five are set up for the calculation
      void setSymmetrized(bool);
      void setQvalues(double,double,double);
      void setSolventScatteringDensity(double);
      void setRealSpaceValues(const Points&);
      void createAtomicFormFactors(const AtomTypes&);
      
      void calcFormFactors(NormalDensities,QProgressDialog&,uint);
      
      double qMin() const;
      double qMax() const;
      Points& getFormFactors();
      const Points& getFormFactors() const;
      uint getNumberQValues() const;

      void clear();
      
   private:
   
      void calcPoints();
      void sortPoints();
};


#endif /* FOURIERFORMFACTOR_H_ */
