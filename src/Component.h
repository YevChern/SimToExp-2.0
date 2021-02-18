/*
 * Component.h
 *
 *  Created on: Oct 29, 2012
 *      Author: bholland
 * 
 *  from NK's struct
 * 
 * 	int id;
	char *name;
	char atom[MaxAtom][10];
	int Natom;
	double nat;		// total number of atoms in component
	double nel;		// total number of electrons in component
	double nsl;		// total neutron scattering length in component
	double ngr;		// number of groups per component
	double natpgr;	// number of atoms per component group = nat/ngr
	double nelpgr;	// number of electrons per component group = nel/ngr
	double nslpgr;	// neutron scattering length per component group = nsl/ngr
	double volume;	// of the component group
	short int fix;
//	double xsim[MaxPoint];
    double ysim[MaxPoint];
	double ysimbak[MaxPoint];
    double yprob[MaxPoint];
	double yed[MaxPoint];
	double ynsld[MaxPoint];
	
 */

#ifndef COMPONENT_H_
#define COMPONENT_H_

#include "NormalDensity.h"

typedef QVector<QPointF> Points;
typedef QList<NormalDensity> NormalDensities;

class Component { 

   QString cmp_name;
   NormalDensities particles;
   
   Points num_densities;
   Points elec_densities;
   Points neutron_SL_densities;
   Points probabilities;
   
   bool probs_calculated;

   double total_electrons;
   double total_neutron_SL;
   double number_groups;
   double cmp_volume;
   
   public:
   
      Component();
      Component(const QString&,double);
      virtual ~Component();
      
      void addParticle(const NormalDensity&);
      void calcDensities();
      void setVolume(double);
      
      //getters
      QString name();
      const QString& name() const;
      double electrons() const; //number of electrons
      double neutronSL() const; //neutron scattering length
      double numberOfGroups() const;
      double volume() const;
      uint numberOfAtoms() const;
      
      bool isFixed() const;
      
      Points& getNumberDensities();
      const Points& getNumberDensities() const;
      Points& getElectronDensities();
      const Points& getElectronDensities() const;
      Points& getNeutronSLDensities();
      const Points& getNeutronSLDensities() const;
      Points& getProbabilities();
      
      //calculated values
      double atomsPerGroup() const;
      double electronsPerGroup() const;
      double neutronSLPerGroup() const;

   private:

      void calcProbabilities();
};

#endif /* COMPONENT_H_ */
