/*
 * AtomicInfo.h
 *
 *  Created on: Oct 11, 2012
 *      Author: bholland
 * 
 * Very simple class that stores atom info; very useful for looking up values!
 * 
 */

#ifndef ATOMICINFO_H_
#define ATOMICINFO_H_ 

#include "AtomType.h"
#include "HSVColour.h"

#include <map>
#include <vector>
#include <string>

typedef std::pair<std::string,AtomType> AtomTypePair;
typedef std::map<std::string,AtomType> InfoAtomTypeMap;

typedef std::pair<std::string,ResidueType> ResidueStringPair;
typedef std::map<std::string,ResidueType> ResidueStringMap;

typedef std::map<std::string,std::string> ResidueTypeMap;
typedef std::map<AtomType,double> ElectronDensityMap, NeutronSLMap, MassMap;
typedef std::pair<std::string,std::string> StringPair;
typedef std::map<AtomType,StringPair> NameAndSymbolMap;
typedef std::map<AtomType,HSVColour> HSVColourMap;
typedef std::map<AtomType,std::vector<double> > AtomicFFMap;
typedef std::pair<AtomType,std::vector<double> > AtomicFFPair;

class AtomicInfo {
   
   InfoAtomTypeMap types;
   ResidueStringMap res_types;
   ResidueTypeMap lipid_residues;
   ElectronDensityMap electrons;
   ElectronDensityMap electrons_default;
   NeutronSLMap neutron_SLs;
   NeutronSLMap neutron_SLs_default;
   MassMap atomic_mass;
   NameAndSymbolMap names_and_symbols;
   HSVColourMap colours;
   AtomicFFMap aff_constants;
   
   public:
      
      AtomicInfo();
      virtual ~AtomicInfo();
      
      AtomType atomType(const std::string&);
      AtomType atomType(std::string, std::string&);
      int numberOfAtomTypes() const;
      std::vector<AtomType> allAtomTypes() const;

      ResidueType residueType(const std::string&);

      double mass(AtomType) const;
      double numberElectrons(AtomType) const;
      double neutronSL(AtomType) const;
      void setNumberElectrons(AtomType,double);
      void setNeutronSL(AtomType,double);
      void setUserDefinedName(AtomType,const std::string&);
      bool setUserDefinedID(AtomType,const std::string&);
      void resetScatteringDefaults();

      StringPair nameAndSymbol(AtomType) const;
      HSVColour colour(AtomType) const;
      const std::vector<double>& atomicFormFactorConstants(AtomType) const;
      double UAFormFactorConstant(AtomType) const;
      
      bool isUnitedAtom(AtomType) const;
      bool isLipid(ResidueType) const;
      
   private:
      
      void buildTypeMap();
      void buildLipidResidueMap();
      void buildMassMap();
      void buildElectronMap();
      void buildNeutronMap();
      void buildNameMap();
      void buildColourMap();
      void buildAFFconstantMap();

      void addColours(uint&,std::vector<HSVColour>&,double,std::vector<AtomType>&);
};

#endif /* ATOMICINFO_H_  */
