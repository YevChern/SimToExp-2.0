/*
 * NormalDensity.h
 *
 *  Created on: Oct 2, 2012
 *      Author: bholland
 */

#ifndef NORMALDENSITY_H_
#define NORMALDENSITY_H_

#include <QtCore/QString>
#include <QtCore/QVector>
#include <QtCore/QPointF>
#include <QtCore/QMap>

#include "AtomType.h"
#include "ParticleLabel.h"

typedef QMap<float,float> DensitiesMap;
typedef QVector<QPointF> Points;

class NormalDensity {
   
   ParticleLabel label;
   AtomType atom_type;
   ResidueType res_type;

   double num_electrons; //number of electrons for 'type'
   double neutron_sl; //neutron scattering length for 'type'
   
   DensitiesMap densities_; //want them ordered by 'z' position
   DensitiesMap::iterator densities_itr;

   //The following are store for quick look up; use a boolean to keep track if up to date
   //Important: since these are for final output, all points along independent axis
   //are shifted half a bin size to the right - this is because the keys in the density map
   //refer to the lower end (left side) of the bin, but should plot center of the bin
   Points vec_densities; //for convenience
   Points sym_densities; //for convenience
   Points elec_densities; //stored for convenience
   Points neut_densities; //stored for convenience
   
   double bin_size;

   bool all_densities_updated;
   bool bin_size_updated;

   public:
      
   	  NormalDensity();
   	  NormalDensity(const ParticleLabel&);
      virtual ~NormalDensity();
      
      QString particleName();
      const QString& particleName() const;
      QString residueName();
      const QString& residueName() const;

      ParticleLabel particleLabel();
      const ParticleLabel& particleLabel() const;

      void setParticleName(const QString&);
      void setResidueName(const QString&);

      AtomType atomType() const;
      ResidueType residueType() const;
      bool residueIsLipid() const;

      double numberOfElectrons() const;
      double neutronSL() const;
      bool densitiesCalculated() const;
      double binSize() const;

      void setLabel(const ParticleLabel&);
      void setAtomInfo(AtomType,double,double);
      void setResidueType(ResidueType);
      
      //iterators for map
      QPointF firstPoint(); //only for number densities
      QPointF firstNonZeroPoint();
      QPointF nextPoint(); //only for number densities
      QPointF lastPoint();  //only for number densities
      QPointF lastNonZeroPoint();
      bool densitiesAtEnd() const;
      bool isEmpty() const;

      //direct access to Point containers
      //should NOT BE USED with trajectory NormalDensity objects, only the final result for convenience,
      //as this can add a LOT of extra memory usage
      Points& densities();
      Points& symmetricDensities();
      Points& electronDensities();
      Points& neutronDensities();
      
      //gets the number densities from the map - only use this for trajectory objects
      Points densitiesFromMap();

      void addDensity(double,double); //only need to add them for now
      void addDensity(QPointF);
      void addDensities(NormalDensity);

      double densityCentre();

      void scaleIndependentAxis(double);
      void scaleDependentAxis(double);
      void shiftIndependentAxis(double);
      void trimIndependentAxis(float,float);

      void calcAllDensities(bool);
      void clearDensities();

      bool operator== (const NormalDensity&) const;
      bool operator!= (const NormalDensity&) const;

   private:

      void setBinSize();
};

#endif /* PARTICLEDENSITY_H_ */ 
