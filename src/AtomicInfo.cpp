
#include "AtomicInfo.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

//because I am lazy
typedef unsigned int uint;

//Default constructor - builds the database upon construction
AtomicInfo::AtomicInfo() {

	//build all of the databases
	buildTypeMap();
	buildLipidResidueMap();
    buildMassMap();
	buildElectronMap();
	buildNeutronMap();
	buildNameMap();
	buildColourMap();
	buildAFFconstantMap();
}

//Destructor
AtomicInfo::~AtomicInfo() {}

//Returns the atom type for the given string, case sensitive!  If not found, returns 'UNKNOWN'.
AtomType AtomicInfo::atomType(const string& ref) {

	InfoAtomTypeMap::iterator itr = types.find(ref);

	//if found, return it
	if(itr != types.end())
		return itr->second;

	return UNKNOWN;
}

//Overloaded AtomType search method; this searches for the residue amongst the lipid residues and
//if found adds the appropriate prefix to the particle name.  If nothing found, return 'UNKNOWN'
AtomType AtomicInfo::atomType(string residue, string& particle) {

	//first try to find particle, if found return right away, if not try prefix
	AtomType type = atomType(particle);
	if(type != UNKNOWN)
		return type;

	ResidueTypeMap::iterator res_itr = lipid_residues.find(residue);

	//if found, add the prefix to particle
	if(res_itr != lipid_residues.end())
		particle = res_itr->second + "_" + particle;

	//now look for the particle
	return atomType(particle);
}

//Returns the number of AtomTypes in the map
int AtomicInfo::numberOfAtomTypes() const {

	return types.size();
}

//Returns a vector container with all of the AtomType enums available for use
vector<AtomType> AtomicInfo::allAtomTypes() const {

	//use a set first to ensure the types are only added once
	set<AtomType> type_set;

	//iterate through the map, adding to the vector for each value
	InfoAtomTypeMap::const_iterator itr, itr_end = types.end();
	for(itr = types.begin(); itr != itr_end; ++itr) {

		type_set.insert(itr->second);
	}

	//now that all unique types have been added, dump them in a vector

	set<AtomType>::iterator set_itr, set_itr_end = type_set.end();
	vector<AtomType> type_list;
	for(set_itr = type_set.begin(); set_itr != set_itr_end; ++set_itr) {

		type_list.push_back(*set_itr);
	}

	return type_list;
}

//Returns the lipid type for the given string, case sensitive!  If not found, returns 'UNKNOWN_LIPID'.
ResidueType AtomicInfo::residueType(const string& ref) {

	ResidueStringMap::iterator itr = res_types.find(ref);

	//if found, return it
	if(itr != res_types.end())
		return itr->second;

	return UNKNOWN_RESIDUE;
}

//returns mass associated with the given AtomType. If 'AtomType'
//is not found, returns 1
double AtomicInfo::mass(AtomType type) const {

    ElectronDensityMap::const_iterator itr = atomic_mass.find(type);

    //if found, return it
    if(itr != atomic_mass.end())
        return itr->second;

    return 1;
}

//returns the number of electrons associated with the given AtomType. If 'AtomType'
//is not found, returns 0
double AtomicInfo::numberElectrons(AtomType type) const {

	ElectronDensityMap::const_iterator itr = electrons.find(type);

	//if found, return it
	if(itr != electrons.end())
		return itr->second;

	return 0;
}

//returns the neutron scattering length associated with the given AtomType. If 'AtomType'
//is not found, returns 0
double AtomicInfo::neutronSL(AtomType type) const {

	NeutronSLMap::const_iterator itr = neutron_SLs.find(type);

	//if found, return it
	if(itr != neutron_SLs.end())
		return itr->second;

	return 0;
}

//Changes the number of electrons for the 'AtomType' to the value provided, partial charges allowed
void AtomicInfo::setNumberElectrons(AtomType type, double num_elecs) {

	//if 'type' not found, do nothing
	ElectronDensityMap::iterator itr = electrons.find(type);

	//if found, set to num_elecs
	if(itr != electrons.end())
		itr->second = num_elecs;
}

//Changes the neutron scattering length for the 'AtomType' to the value provided
void AtomicInfo::setNeutronSL(AtomType type, double neutron_SL) {

	//if 'type' not found, do nothing
	NeutronSLMap::iterator itr = neutron_SLs.find(type);

	//if found, set to num_elecs
	if(itr != neutron_SLs.end())
		itr->second = neutron_SL;
}

//Allows the user to change the names of the 'user defined' particles.  If the 'type' is not
//found, quietly returns without doing anything.
void AtomicInfo::setUserDefinedName(AtomType type, const string& name) {

	//get the current id
	NameAndSymbolMap::iterator itr = names_and_symbols.find(type);
	if(itr == names_and_symbols.end())
		return;

	//otherwise make the change to the new values
	names_and_symbols[type] = StringPair(name,itr->second.second);
}

//Allows the user to change the ID of the 'user defined' particles, but not with one that already exists
bool AtomicInfo::setUserDefinedID(AtomType type, const string& id) {

	//check to make sure the ID doesn't already exist, if it does return false
	if(types.find(id) != types.end())
		return false;

	//get the current name
	const string& name = names_and_symbols[type].first;

	//find the key 'type' and set it to the new user defined values
	names_and_symbols[type] = StringPair(name,id);
	return true;
}


//Changes the current scattering lengths back to the default values provided when the class is constructed
void AtomicInfo::resetScatteringDefaults() {

	//remove the NON-user defined particle values
	electrons.erase(electrons.begin(), electrons.find(USER_1));
	neutron_SLs.erase(neutron_SLs.begin(), neutron_SLs.find(USER_1));

	//only want to add back the default values for NON-user defined particles
	electrons.insert(electrons_default.begin(), electrons_default.end());
	neutron_SLs.insert(neutron_SLs_default.begin(), neutron_SLs_default.end());
}

//return the name and symbol of the given AtomType. If AtomType not found,
//returns an empty string and '?'
StringPair AtomicInfo::nameAndSymbol(AtomType type) const {

	NameAndSymbolMap::const_iterator itr = names_and_symbols.find(type);

	//if found, return it
	if(itr != names_and_symbols.end())
		return itr->second;

	//else return blank string and question mark
	return StringPair("","?");
}

//returns an HSV representation of the assigned colour in a standard C++ vector object.  Values
//range from 0 to 1, unlike the common integer version that uses 0 to 255.
HSVColour AtomicInfo::colour(AtomType type) const {

	map<AtomType,HSVColour>::const_iterator itr = colours.find(type);

	//if found, return it
	if(itr != colours.end())
		return itr->second;

	//else return white
	return HSVColour(0.0,0.0,1.0);
}

//returns pointer to the set of constants for calculating the form factor for
//the given AtomType.  Returns an empty vector if not found
const vector<double>& AtomicInfo::atomicFormFactorConstants(AtomType type) const {

	map<AtomType,vector<double> >::const_iterator itr = aff_constants.find(type);

    //if not found, put out a cerr message, will probably crash here if not SYSTM_TOTAL
    if((itr == aff_constants.end()) && (type != SYSTEM_TOTAL))
		cerr << "XrayAtomicFormFactor::getAtomicFormFactorConstants: AtomType not found in map, iterator pointing to end" << endl;

	//else return
	return itr->second;
}

//returns the 'c' constant for the UA form factor calculation. 'c' is same
//as the electron number!
double AtomicInfo::UAFormFactorConstant(AtomType type) const {

	return numberElectrons(type);
}

//Returns true if the AtomType is a united atom particle,
bool AtomicInfo::isUnitedAtom(AtomType type) const {

	//All UAs (including MARTINI now) come together in the enumerator, so use this to determine
	if(type < UA_METHYL or type > MART_D)
		return false;
	return true;
}

//Return true if the give residue is a lipid
bool AtomicInfo::isLipid(ResidueType type) const  {

	if(type < DHPC or type > BOLB)
		return false;
	return true;
}

//helper function - creates map for retrieving the AtomType from the string that represents it in the .sim file and/or trajectory files
void AtomicInfo::buildTypeMap() {

	//atoms
	types.insert(AtomTypePair("H", HYDROGEN));
	types.insert(AtomTypePair("D", DEUTERIUM));
	types.insert(AtomTypePair("C", CARBON));
	types.insert(AtomTypePair("N", NITROGEN));
	types.insert(AtomTypePair("O", OXYGEN));
	types.insert(AtomTypePair("P", PHOSPHORUS));
	//ions - should never have atomic sodium/chlorine -> would need to revisit definitions below if so
	types.insert(AtomTypePair("NA", SODIUM_PLUS));
	types.insert(AtomTypePair("Na", SODIUM_PLUS));
	types.insert(AtomTypePair("NA+", SODIUM_PLUS));
	types.insert(AtomTypePair("Na+", SODIUM_PLUS));
	types.insert(AtomTypePair("CL", CHLORIDE));
	types.insert(AtomTypePair("Cl", CHLORIDE));
	types.insert(AtomTypePair("CL-", CHLORIDE));
	types.insert(AtomTypePair("Cl-", CHLORIDE));

	//united atoms
	types.insert(AtomTypePair("CH", UA_METHINE));
	types.insert(AtomTypePair("CH2", UA_METHYLENE));
	types.insert(AtomTypePair("CH3", UA_METHYL));
	types.insert(AtomTypePair("H2O", UA_WATER));
	types.insert(AtomTypePair("D2O", UA_HEAVY_WATER));

	//martini
	types.insert(AtomTypePair("NH3", MART_NH3));
	types.insert(AtomTypePair("NC1", MART_NC3));
	types.insert(AtomTypePair("NC2", MART_NC3));
	types.insert(AtomTypePair("NC3", MART_NC3));
	types.insert(AtomTypePair("PO1", MART_PO4));
	types.insert(AtomTypePair("PO2", MART_PO4));
	types.insert(AtomTypePair("PO4", MART_PO4));
	types.insert(AtomTypePair("CNO", MART_CNO));
	types.insert(AtomTypePair("GL0", MART_GL0));
	types.insert(AtomTypePair("GL1", MART_GL1));
	types.insert(AtomTypePair("GL2", MART_GL2));
	types.insert(AtomTypePair("GL3", MART_GL3));
	types.insert(AtomTypePair("GL4", MART_GL4));

	//Lipids - tails, unfortunately each molecule can have its own representation of the "atom" types
	//DP
	types.insert(AtomTypePair("DP_C1A", MART_PC1));
	types.insert(AtomTypePair("DP_C2A", MART_PC2));
	types.insert(AtomTypePair("DP_C3A", MART_PC3));
	types.insert(AtomTypePair("DP_C4A", MART_PC4));
	types.insert(AtomTypePair("DP_C1B", MART_PC1));
	types.insert(AtomTypePair("DP_C2B", MART_PC2));
	types.insert(AtomTypePair("DP_C3B", MART_PC3));
	types.insert(AtomTypePair("DP_C4B", MART_PC4));

	//DH
	types.insert(AtomTypePair("DH_C1A", MART_HC1));
	types.insert(AtomTypePair("DH_C2A", MART_HC2));
	types.insert(AtomTypePair("DH_C1B", MART_HC1));
	types.insert(AtomTypePair("DH_C2B", MART_HC2));

	//DL
	types.insert(AtomTypePair("DL_C1A", MART_LC1));
	types.insert(AtomTypePair("DL_C2A", MART_LC2));
	types.insert(AtomTypePair("DL_C3A", MART_LC3));
	types.insert(AtomTypePair("DL_C1B", MART_LC1));
	types.insert(AtomTypePair("DL_C2B", MART_LC2));
	types.insert(AtomTypePair("DL_C3B", MART_LC3));

	//DM
	types.insert(AtomTypePair("DM_C1A", MART_MC1));
	types.insert(AtomTypePair("DM_C2A", MART_MC2));
	types.insert(AtomTypePair("DM_C3A", MART_MC3));
	types.insert(AtomTypePair("DM_C1B", MART_MC1));
	types.insert(AtomTypePair("DM_C2B", MART_MC2));
	types.insert(AtomTypePair("DM_C3B", MART_MC3));

	//DS
	types.insert(AtomTypePair("DS_C1A", MART_SC1));
	types.insert(AtomTypePair("DS_C2A", MART_SC2));
	types.insert(AtomTypePair("DS_C3A", MART_SC3));
	types.insert(AtomTypePair("DS_C4A", MART_SC4));
	types.insert(AtomTypePair("DS_C5A", MART_SC5));
	types.insert(AtomTypePair("DS_C1B", MART_SC1));
	types.insert(AtomTypePair("DS_C2B", MART_SC2));
	types.insert(AtomTypePair("DS_C3B", MART_SC3));
	types.insert(AtomTypePair("DS_C4B", MART_SC4));
	types.insert(AtomTypePair("DS_C5B", MART_SC5));

	//PO
	types.insert(AtomTypePair("PO_C1A", MART_PC1));
	types.insert(AtomTypePair("PO_C2A", MART_PC2));
	types.insert(AtomTypePair("PO_C3A", MART_PC3));
	types.insert(AtomTypePair("PO_C4A", MART_PC4));
	types.insert(AtomTypePair("PO_C1B", MART_OC1));
	types.insert(AtomTypePair("PO_C2B", MART_OC2));
	types.insert(AtomTypePair("PO_D3B", MART_OD3));
	types.insert(AtomTypePair("PO_C4B", MART_OC4));
	types.insert(AtomTypePair("PO_C5B", MART_OC5));

	//DO
	types.insert(AtomTypePair("DO_C1A", MART_OC1));
	types.insert(AtomTypePair("DO_C2A", MART_OC2));
	types.insert(AtomTypePair("DO_D3A", MART_OD3));
	types.insert(AtomTypePair("DO_C4A", MART_OC4));
	types.insert(AtomTypePair("DO_C5A", MART_OC5));
	types.insert(AtomTypePair("DO_C1B", MART_OC1));
	types.insert(AtomTypePair("DO_C2B", MART_OC2));
	types.insert(AtomTypePair("DO_D3B", MART_OD3));
	types.insert(AtomTypePair("DO_C4B", MART_OC4));
	types.insert(AtomTypePair("DO_C5B", MART_OC5));

	//DA
	types.insert(AtomTypePair("DA_D1A", MART_AD1));
	types.insert(AtomTypePair("DA_D2A", MART_AD2));
	types.insert(AtomTypePair("DA_D3A", MART_AD3));
	types.insert(AtomTypePair("DA_D4A", MART_AD4));
	types.insert(AtomTypePair("DA_C5A", MART_AC5));
	types.insert(AtomTypePair("DA_D1B", MART_AD1));
	types.insert(AtomTypePair("DA_D2B", MART_AD2));
	types.insert(AtomTypePair("DA_D3B", MART_AD3));
	types.insert(AtomTypePair("DA_D4B", MART_AD4));
	types.insert(AtomTypePair("DA_C5B", MART_AC5));

	//DU
	types.insert(AtomTypePair("DU_C1A", MART_UC1));
	types.insert(AtomTypePair("DU_D2A", MART_UD2));
	types.insert(AtomTypePair("DU_D3A", MART_UD3));
	types.insert(AtomTypePair("DU_C4A", MART_UC4));
	types.insert(AtomTypePair("DU_C1B", MART_UC1));
	types.insert(AtomTypePair("DU_D2B", MART_UD2));
	types.insert(AtomTypePair("DU_D3B", MART_UD3));
	types.insert(AtomTypePair("DU_C4B", MART_UC4));

	//BOLA
	types.insert(AtomTypePair("BOLA_C1A", MART_BC1));
	types.insert(AtomTypePair("BOLA_C2A", MART_BC2));
	types.insert(AtomTypePair("BOLA_C3A", MART_BC3));
	types.insert(AtomTypePair("BOLA_C4A", MART_BC4));
	types.insert(AtomTypePair("BOLA_C1B", MART_BC1));
	types.insert(AtomTypePair("BOLA_C2B", MART_BC2));
	types.insert(AtomTypePair("BOLA_C3B", MART_BC3));
	types.insert(AtomTypePair("BOLA_C4B", MART_BC4));
	types.insert(AtomTypePair("BOLA_C1C", MART_BC1));
	types.insert(AtomTypePair("BOLA_C2C", MART_BC2));
	types.insert(AtomTypePair("BOLA_C3C", MART_BC3));
	types.insert(AtomTypePair("BOLA_C4C", MART_BC4));
	types.insert(AtomTypePair("BOLA_C1D", MART_BC1));
	types.insert(AtomTypePair("BOLA_C2D", MART_BC2));
	types.insert(AtomTypePair("BOLA_C3D", MART_BC3));
	types.insert(AtomTypePair("BOLA_C4D", MART_BC4));

	//BOLB
	types.insert(AtomTypePair("BOLB_C1A", MART_BC1));
	types.insert(AtomTypePair("BOLB_C2A", MART_BC2));
	types.insert(AtomTypePair("BOLB_C3A", MART_BC3));
	types.insert(AtomTypePair("BOLB_C4A", MART_BC4));
	types.insert(AtomTypePair("BOLB_C1B", MART_PC1));
	types.insert(AtomTypePair("BOLB_C2B", MART_PC2));
	types.insert(AtomTypePair("BOLB_C3B", MART_PC3));
	types.insert(AtomTypePair("BOLB_C4B", MART_PC4));
	types.insert(AtomTypePair("BOLB_C1C", MART_BC1));
	types.insert(AtomTypePair("BOLB_C2C", MART_BC2));
	types.insert(AtomTypePair("BOLB_C3C", MART_BC3));
	types.insert(AtomTypePair("BOLB_C4C", MART_BC4));
	types.insert(AtomTypePair("BOLB_C1D", MART_PC1));
	types.insert(AtomTypePair("BOLB_C2D", MART_PC2));
	types.insert(AtomTypePair("BOLB_C3D", MART_PC3));
	types.insert(AtomTypePair("BOLB_C4D", MART_PC4));

	//water and polarizable water
	//note: MART_D are just charges, only MART_POL should be considered later
	types.insert(AtomTypePair("W", MART_WATER));
	types.insert(AtomTypePair("PW", MART_POL));
	types.insert(AtomTypePair("WP", MART_D));
	types.insert(AtomTypePair("WM", MART_D));
	types.insert(AtomTypePair("MART_D2O", MART_D2O));

	//set the user defined types to start with a number, can be changed later
	types.insert(AtomTypePair("1", USER_1));
	types.insert(AtomTypePair("2", USER_2));
	types.insert(AtomTypePair("3", USER_3));
	types.insert(AtomTypePair("4", USER_4));
	types.insert(AtomTypePair("5", USER_5));
	types.insert(AtomTypePair("6", USER_6));
	types.insert(AtomTypePair("7", USER_7));
	types.insert(AtomTypePair("8", USER_8));
	types.insert(AtomTypePair("9", USER_9));

	types.insert(AtomTypePair("NA", SODIUM_PLUS));
	types.insert(AtomTypePair("NA+", SODIUM_PLUS));
	types.insert(AtomTypePair("CL", CHLORIDE));
	types.insert(AtomTypePair("CL-", CHLORIDE));

	types.insert(AtomTypePair("System Total", SYSTEM_TOTAL));

	//done with Atoms, also build lipid type map
	res_types.insert(ResidueStringPair("DHPC", DHPC));
	res_types.insert(ResidueStringPair("DLPC", DLPC));
	res_types.insert(ResidueStringPair("DPPC", DPPC));
	res_types.insert(ResidueStringPair("DMPC", DMPC));
	res_types.insert(ResidueStringPair("DSPC", DSPC));
	res_types.insert(ResidueStringPair("DAPC", DAPC));
	res_types.insert(ResidueStringPair("DUPC", DUPC));
	res_types.insert(ResidueStringPair("DOPC", DOPC));
	res_types.insert(ResidueStringPair("POPC", POPC));

	res_types.insert(ResidueStringPair("DPPE", DPPE));
	res_types.insert(ResidueStringPair("DHPE", DHPE));
	res_types.insert(ResidueStringPair("DLPE", DLPE));
	res_types.insert(ResidueStringPair("DMPE", DMPE));
	res_types.insert(ResidueStringPair("DSPE", DSPE));
	res_types.insert(ResidueStringPair("DOPE", DOPE));
	res_types.insert(ResidueStringPair("POPE", POPE));

	res_types.insert(ResidueStringPair("PPCS", PPCS));

	res_types.insert(ResidueStringPair("DOPG", DOPG));
	res_types.insert(ResidueStringPair("POPG", POPG));
	res_types.insert(ResidueStringPair("DOPS", DOPS));
	res_types.insert(ResidueStringPair("POPS", POPS));

	res_types.insert(ResidueStringPair("BOLA", BOLA));
	res_types.insert(ResidueStringPair("BOLB", BOLB));

	res_types.insert(ResidueStringPair("W", WATER));
	res_types.insert(ResidueStringPair("PW", WATER));
	res_types.insert(ResidueStringPair("SOL", WATER));
	res_types.insert(ResidueStringPair("D2O", HEAVY_WATER));

	res_types.insert(ResidueStringPair("ION", ION));
	res_types.insert(ResidueStringPair("NA", ION));
	res_types.insert(ResidueStringPair("NA+", ION));
	res_types.insert(ResidueStringPair("CL", ION));
	res_types.insert(ResidueStringPair("CL-", ION));
}

//helper function - creates container that maps a given residue name (string) to the prefix that
//is used to uniquely identify the appropriate AtomType. Example: for POPC, the prefix is PO, since
//in general MARTINI only has different values for the same particle names (e.g. C1A) in the lipid tail region.
//NOTE: can't simply strip the first two chars, some prefixes are specific to the residue
void AtomicInfo::buildLipidResidueMap() {

	//inserted in same order as in martini_v2.0_lipids.itp
	lipid_residues.insert(StringPair("DPPC","DP"));
	lipid_residues.insert(StringPair("DHPC","DH"));
	lipid_residues.insert(StringPair("DLPC","DL"));
	lipid_residues.insert(StringPair("DMPC","DM"));
	lipid_residues.insert(StringPair("DSPC","DS"));
	lipid_residues.insert(StringPair("POPC","PO"));
	lipid_residues.insert(StringPair("DOPC","DO"));
	lipid_residues.insert(StringPair("DAPC","DA"));
	lipid_residues.insert(StringPair("DUPC","DA"));

	lipid_residues.insert(StringPair("DPPE","DP"));
	lipid_residues.insert(StringPair("DHPE","DH"));
	lipid_residues.insert(StringPair("DLPE","DL"));
	lipid_residues.insert(StringPair("DMPE","DM"));
	lipid_residues.insert(StringPair("DSPE","DS"));
	lipid_residues.insert(StringPair("POPE","PO"));
	lipid_residues.insert(StringPair("DOPE","DO"));

	lipid_residues.insert(StringPair("DOPG","DO"));
	lipid_residues.insert(StringPair("POPG","PO"));
	lipid_residues.insert(StringPair("DOPS","DO"));
	lipid_residues.insert(StringPair("POPS","PO"));

	lipid_residues.insert(StringPair("BOLA","BOLA"));
	lipid_residues.insert(StringPair("BOLB","BOLB"));
}

//creates map for retrieving the mass for the given AtomType
//only for atomic representation for now..
void AtomicInfo::buildMassMap() {

    //both maps initially contain the default values, so set default and use these
    //to loop through and set the standard container

    //atoms and ions
    atomic_mass.insert(pair<AtomType,double>(HYDROGEN,1.008));
    atomic_mass.insert(pair<AtomType,double>(DEUTERIUM,2.014));
    atomic_mass.insert(pair<AtomType,double>(CARBON,12.011));
    atomic_mass.insert(pair<AtomType,double>(NITROGEN,14.007));
    atomic_mass.insert(pair<AtomType,double>(OXYGEN,15.999));
    atomic_mass.insert(pair<AtomType,double>(FLUORINE,18.998));
    atomic_mass.insert(pair<AtomType,double>(SODIUM,22.99));
    atomic_mass.insert(pair<AtomType,double>(MAGNESIUM,24.305));
    atomic_mass.insert(pair<AtomType,double>(PHOSPHORUS,30.973762));
    atomic_mass.insert(pair<AtomType,double>(CHLORINE,35.45));
    atomic_mass.insert(pair<AtomType,double>(POTASSIUM,39.0983));
    atomic_mass.insert(pair<AtomType,double>(CALCIUM,40.078));
}

//helper function - creates map for retrieving the number of electrons for the given AtomType
//values are doubles in case someone wants to use partial charges
void AtomicInfo::buildElectronMap() {

	//both maps initially contain the default values, so set default and use these
	//to loop through and set the standard container

	//atoms and ions
	electrons_default.insert(pair<AtomType,double>(HYDROGEN,1));
	electrons_default.insert(pair<AtomType,double>(DEUTERIUM,1));
	electrons_default.insert(pair<AtomType,double>(CARBON,6));
	electrons_default.insert(pair<AtomType,double>(NITROGEN,7));
	electrons_default.insert(pair<AtomType,double>(OXYGEN,8));
	electrons_default.insert(pair<AtomType,double>(FLUORINE,9));
	electrons_default.insert(pair<AtomType,double>(FLUORIDE,10));
	electrons_default.insert(pair<AtomType,double>(SODIUM_PLUS,10));
	electrons_default.insert(pair<AtomType,double>(MAGNESIUM_TWOPLUS,10));
	electrons_default.insert(pair<AtomType,double>(SODIUM,11));
	electrons_default.insert(pair<AtomType,double>(MAGNESIUM,12));
	electrons_default.insert(pair<AtomType,double>(PHOSPHORUS,15));
	electrons_default.insert(pair<AtomType,double>(CHLORINE,17));
	electrons_default.insert(pair<AtomType,double>(CHLORIDE,18));
	electrons_default.insert(pair<AtomType,double>(POTASSIUM_PLUS,18));
	electrons_default.insert(pair<AtomType,double>(CALCIUM_TWOPLUS,18));
	electrons_default.insert(pair<AtomType,double>(POTASSIUM,19));
	electrons_default.insert(pair<AtomType,double>(CALCIUM,20));

	//united atom particles
	electrons_default.insert(pair<AtomType,double>(UA_METHINE,7)); //CH
	electrons_default.insert(pair<AtomType,double>(UA_METHYLENE,8)); //CH2
	electrons_default.insert(pair<AtomType,double>(UA_METHYL,9)); //CH3
	electrons_default.insert(pair<AtomType,double>(UA_WATER,10)); //H2O
	electrons_default.insert(pair<AtomType,double>(UA_HEAVY_WATER,10)); //D2O

	//martini lipid particles
	electrons_default.insert(pair<AtomType,double>(MART_NH3,26)); //NH3+
	electrons_default.insert(pair<AtomType,double>(MART_NC3,50)); //NC3+
	electrons_default.insert(pair<AtomType,double>(MART_PO4,47)); //PO4-
	electrons_default.insert(pair<AtomType,double>(MART_CNO,40)); //Serine (e.g. from POPS)
	electrons_default.insert(pair<AtomType,double>(MART_GL0,41)); //Glycerol for headgroups (e.g. POPG)
	electrons_default.insert(pair<AtomType,double>(MART_GL1,23)); //Part 1 of diglyceride
	electrons_default.insert(pair<AtomType,double>(MART_GL2,16)); //Part 2 of diglyceride
	electrons_default.insert(pair<AtomType,double>(MART_GL3,23)); //Part 1 of diglyceride, other side of "bola" (i.e. cyclic) lipid
	electrons_default.insert(pair<AtomType,double>(MART_GL4,16)); //Part 2 of diglyceride, other side of "bola" (i.e. cyclic) lipid

	//Lipids - tails -> Number of atoms in each group is different for each molecule, thus why it is molecule specific
	//PALMITOYL
	electrons_default.insert(pair<AtomType,double>(MART_PC1,38)); //O=C-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_PC2,32)); //CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_PC3,32)); //CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_PC4,33)); //CH2-CH2-CH2-CH3

	//HEXANOYL
	electrons_default.insert(pair<AtomType,double>(MART_HC1,30)); //O=C-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_HC2,25)); //CH2-CH2-CH3

	//LAUROYL
	electrons_default.insert(pair<AtomType,double>(MART_LC1,38)); //O=C-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_LC2,32)); //CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_LC3,33)); //CH2-CH2-CH2-CH3

	//MYRISTOYL - number of electrons is the average over the three possible bead definitions
	electrons_default.insert(pair<AtomType,double>(MART_MC1,43 + 1/3.0)); //O=C-CH2-CH2-CH2           O=C-CH2-CH2-CH2-CH2      O=C-CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_MC2,37 + 1/3.0)); //CH2-CH2-CH2-CH2-CH2   or  CH2-CH2-CH2-CH2      or  CH2-CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_MC3,38 + 1/3.0)); //CH2-CH2-CH2-CH2-CH3       CH2-CH2-CH2-CH2-CH3      CH2-CH2-CH2-CH3

	//STEAROYL - number of electrons is the average over the two possible symmetric bead definitions
	electrons_default.insert(pair<AtomType,double>(MART_SC1,34)); //O=C-CH2-CH2          O=C-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_SC2,28)); //CH2-CH2-CH2-CH2      CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_SC3,32)); //CH2-CH2-CH2-CH2  or  CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_SC4,28)); //CH2-CH2-CH2-CH2      CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_SC5,29)); //CH2-CH2-CH3          CH2-CH2-CH2-CH3

	//OLEOYL
	electrons_default.insert(pair<AtomType,double>(MART_OC1,34)); //O=C-CH2-CH2          O=C-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_OC2,28)); //CH2-CH2-CH2-CH2      CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_OD3,30)); //CH2-CH=CH-CH2    or  CH2-CH=CH-CH2
	electrons_default.insert(pair<AtomType,double>(MART_OC4,28)); //CH2-CH2-CH2-CH2      CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_OC5,29)); //CH2-CH2-CH3          CH2-CH2-CH2-CH3

	//ARACHIDONOYL
	electrons_default.insert(pair<AtomType,double>(MART_AD1,52)); //O=C-CH2-CH2-CH2-CH=CH
	electrons_default.insert(pair<AtomType,double>(MART_AD2,22)); //CH2-CH=CH
	electrons_default.insert(pair<AtomType,double>(MART_AD3,22)); //CH2-CH=CH
	electrons_default.insert(pair<AtomType,double>(MART_AD4,30)); //CH2-CH=CH-CH2
	electrons_default.insert(pair<AtomType,double>(MART_AC5,33)); //CH2-CH2-CH2-CH3

	//LINOLEYL
	electrons_default.insert(pair<AtomType,double>(MART_UC1,46)); //O=C-CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_UD2,38)); //CH2-CH2-CH2-CH=CH
	electrons_default.insert(pair<AtomType,double>(MART_UD3,30)); //CH2-CH=CH-CH2
	electrons_default.insert(pair<AtomType,double>(MART_UC4,33)); //CH2-CH2-CH2-CH3

	//Cyclic PALMITOYL - (e.g. BOLA in martini itp file)
	electrons_default.insert(pair<AtomType,double>(MART_BC1,38)); //O=C-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_BC2,32)); //CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_BC3,32)); //CH2-CH2-CH2-CH2
	electrons_default.insert(pair<AtomType,double>(MART_BC4,32)); //CH2-CH2-CH2-CH2

	//water and polarizable water
	//note: MART_D are just charges, only MART_POL should be considered later
	electrons_default.insert(pair<AtomType,double>(MART_WATER,40)); //supposed to replace 4 water molecule
	electrons_default.insert(pair<AtomType,double>(MART_POL,40)); //TODO: make sure this is 4 H2O as well
	electrons_default.insert(pair<AtomType,double>(MART_D,0));
	electrons_default.insert(pair<AtomType,double>(MART_D,0));
	electrons_default.insert(pair<AtomType,double>(MART_D2O,40));


	ElectronDensityMap::iterator itr, itr_end = electrons_default.end();
	for(itr = electrons_default.begin(); itr != itr_end; ++itr) {

		electrons.insert(*itr);
	}

	//special user defined types need to be added last, no default values
	electrons.insert(pair<AtomType,double>(USER_1,0));
	electrons.insert(pair<AtomType,double>(USER_2,0));
	electrons.insert(pair<AtomType,double>(USER_3,0));
	electrons.insert(pair<AtomType,double>(USER_4,0));
	electrons.insert(pair<AtomType,double>(USER_5,0));
	electrons.insert(pair<AtomType,double>(USER_6,0));
	electrons.insert(pair<AtomType,double>(USER_7,0));
	electrons.insert(pair<AtomType,double>(USER_8,0));
	electrons.insert(pair<AtomType,double>(USER_9,0));
}

//helper function - creates map for retrieving the neutron scattering length for the given AtomType
void AtomicInfo::buildNeutronMap() {

	//both maps initially contain the default values
	neutron_SLs_default.insert(pair<AtomType,double>(HYDROGEN,-3.7409E-5)); //taken from: "Koester, L., Nistier, W.: Z. Phys. A 272 (1975) 189."
	neutron_SLs_default.insert(pair<AtomType,double>(DEUTERIUM,6.67E-5));
	neutron_SLs_default.insert(pair<AtomType,double>(CARBON,6.6484E-5)); //taken from: "Koester, L., Nistier, W.: Z. Phys. A 272 (1975) 189."
	neutron_SLs_default.insert(pair<AtomType,double>(NITROGEN,9.36E-5)); //taken from: "Koester, L., Knopf, K., Waschkowski, W.: Z. Phys. A 277 (1976) 77."
	neutron_SLs_default.insert(pair<AtomType,double>(OXYGEN,5.805E-5)); //taken from: "Koester, L., Knopf, K., Waschkowski, W.: Z. Phys. A 292 (1979) 95."
	neutron_SLs_default.insert(pair<AtomType,double>(SODIUM,3.63E-5)); //taken from: "Koester, L., 1972" - unknown paper, possible "Koester, L., Knopf, K.: Z. Naturforschung 27A (1972) 901."
	neutron_SLs_default.insert(pair<AtomType,double>(SODIUM_PLUS,3.63E-5));
	neutron_SLs_default.insert(pair<AtomType,double>(PHOSPHORUS,5.17E-5));
	neutron_SLs_default.insert(pair<AtomType,double>(CHLORINE,9.5792E-5)); //includes both isotopes (35,37), taken from: "Koester, L., Nistier, W.: Z. Phys. A 272 (1975) 189."
	neutron_SLs_default.insert(pair<AtomType,double>(CHLORIDE,9.57926E-5));
	neutron_SLs_default.insert(pair<AtomType,double>(UA_METHYL,-4.5743E-05));
	neutron_SLs_default.insert(pair<AtomType,double>(UA_METHYLENE,-0.8334E-5)); //United atom values simply the sum of the contents
	neutron_SLs_default.insert(pair<AtomType,double>(UA_METHINE,2.9075E-05));
	neutron_SLs_default.insert(pair<AtomType,double>(UA_WATER,-1.6768E-05));
	neutron_SLs_default.insert(pair<AtomType,double>(UA_HEAVY_WATER,19.145E-5));

	//use defaults to insert into regular container
	map<AtomType,double>::iterator itr, itr_end = neutron_SLs_default.end();
	for(itr = neutron_SLs_default.begin(); itr != itr_end; ++itr) {

		neutron_SLs.insert(*itr);
	}

	//user defined not in default values
	neutron_SLs.insert(pair<AtomType,double>(USER_1,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_2,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_3,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_4,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_5,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_6,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_7,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_8,0));
	neutron_SLs.insert(pair<AtomType,double>(USER_9,0));
}

//helper function - creates map for the string name and the corresponding short string representation for the AtomType
void AtomicInfo::buildNameMap() {

	//atomic
	names_and_symbols.insert(pair<AtomType,StringPair >(HYDROGEN,StringPair("Hydrogen","H")));
	names_and_symbols.insert(pair<AtomType,StringPair >(DEUTERIUM,StringPair("Deuterium","D")));
	names_and_symbols.insert(pair<AtomType,StringPair >(CARBON,StringPair("Carbon","C")));
	names_and_symbols.insert(pair<AtomType,StringPair >(NITROGEN,StringPair("Nitrogen","N")));
	names_and_symbols.insert(pair<AtomType,StringPair >(OXYGEN,StringPair("Oxygen","O")));
	names_and_symbols.insert(pair<AtomType,StringPair >(SODIUM,StringPair("Sodium","Na")));
	names_and_symbols.insert(pair<AtomType,StringPair >(SODIUM_PLUS,StringPair("Sodium-ion","Na+")));
	names_and_symbols.insert(pair<AtomType,StringPair >(PHOSPHORUS,StringPair("Phosphorus","P")));
	names_and_symbols.insert(pair<AtomType,StringPair >(CHLORINE,StringPair("Chlorine","Cl")));
	names_and_symbols.insert(pair<AtomType,StringPair >(CHLORIDE,StringPair("Chloride-ion","Cl-")));

	//united atoms
	names_and_symbols.insert(pair<AtomType,StringPair >(UA_METHINE,StringPair("Methine","CH")));
	names_and_symbols.insert(pair<AtomType,StringPair >(UA_METHYLENE,StringPair("Methylene","CH2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(UA_METHYL,StringPair("Methyl","CH3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(UA_WATER,StringPair("Water","H2O")));
	names_and_symbols.insert(pair<AtomType,StringPair >(UA_HEAVY_WATER,StringPair("Heavy-water","D2O")));

	//MARTINI particles
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_NC3,StringPair("Choline","NC3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_NH3,StringPair("Ammonium","NH3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_PO4,StringPair("Phosphate","PO4")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_CNO,StringPair("Serine","CNO")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_GL0,StringPair("Glycerol","GL0")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_GL1,StringPair("Diglyceride-1","GL1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_GL2,StringPair("Diglyceride-2","GL2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_GL3,StringPair("Diglyceride-3","GL3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_GL4,StringPair("Diglyceride-4","GL4")));

	//Lipids - tails, unfortunately are specific to each molecule!
	//Palmitoyl (e.g. DPPC) - 16:0 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_PC1,StringPair("Palmitoyl-1, 4:0","PC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_PC2,StringPair("Palmitoyl-2, 4:0","PC2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_PC3,StringPair("Palmitoyl-3, 4:0","PC3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_PC4,StringPair("Palmitoyl-4, 4:0","PC4")));

	//Hexanoyl (e.g. DHPC) - 6:0 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_HC1,StringPair("Hexanoyl-1, 3:0","HC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_HC2,StringPair("Hexanoyl-2, 3:0","HC2")));

	//Lauroyl (e.g. DLPC) - 12:0 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_LC1,StringPair("Lauroyl-1, 4:0","LC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_LC2,StringPair("Lauroyl-2, 4:0","LC2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_LC3,StringPair("Lauroyl-3, 4:0","LC3")));

	//Myristoyl (e.g. DMPC) - 14:0 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_MC1,StringPair("Myristoyl-1, 4.33:0","MC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_MC2,StringPair("Myristoyl-2, 4.33:0","MC2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_MC3,StringPair("Myristoyl-3, 4.33:0","MC3")));

	//Stearoyl (e.g. DSPC) - 18:0 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_SC1,StringPair("Stearoyl-1, 3.5:0","SC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_SC2,StringPair("Stearoyl-2, 3.5:0","SC2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_SC3,StringPair("Stearoyl-3, 4:0","SC3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_SC4,StringPair("Stearoyl-4, 3.5:0","SC4")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_SC5,StringPair("Stearoyl-5, 3.5:0","SC5")));

	//Oleoyl (e.g. POPC) - 18:1 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_OC1,StringPair("Oleoyl-1, 3.5:0","OC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_OC2,StringPair("Oleoyl-2, 3.5:0","OC2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_OD3,StringPair("Oleoyl-3, 4:1","OD3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_OC4,StringPair("Oleoyl-4, 3.5:0","OC4")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_OC5,StringPair("Oleoyl-5, 3.5:0","OC5")));

	//Arachidonoyl - 20:4 carbon
	//TODO: At some point in the distant future, test to see if for beads AD4 and AC5 the form factor
	//is better represented by 4:1/4:0 or 3:1/5:0
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_AD1,StringPair("Arachidonoyl-1, 6:1","AD1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_AD2,StringPair("Arachidonoyl-2, 3:1","AD2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_AD3,StringPair("Arachidonoyl-3, 3:1","AD3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_AD4,StringPair("Arachidonoyl-4, 4:1","AD4")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_AC5,StringPair("Arachidonoyl-5, 4:0","AC5")));

	//Linoleyl - 18:2 carbon
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_UC1,StringPair("Linoleyl-1, 5:0","UC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_UD2,StringPair("Linoleyl-2, 5:1","UD2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_UD3,StringPair("Linoleyl-3, 4:1","UD3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_UC4,StringPair("Linoleyl-4, 4:0","UC4")));

	//BOLA - cyclic di-DP -> last carbon is attached to other DP chain
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_BC1,StringPair("Cyc-palmitoyl-1, 4:0","BC1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_BC2,StringPair("Cyc-palmitoyl-2, 4:0","BC2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_BC3,StringPair("Cyc-palmitoyl-3, 4:0","BC3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_BC4,StringPair("Cyc-palmitoyl-4, 4:0","BC4")));

	//martini water and polarizable types
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_WATER,StringPair("Water","W")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_POL,StringPair("Polar-water","PW")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_D,StringPair("Water-charge","D")));
	names_and_symbols.insert(pair<AtomType,StringPair >(MART_D2O,StringPair("Heavy-water","D2O")));

	//System total
	names_and_symbols.insert(pair<AtomType,StringPair >(SYSTEM_TOTAL,StringPair("System total","")));

	//user defined
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_1,StringPair("User-1","1")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_2,StringPair("User-2","2")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_3,StringPair("User-3","3")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_4,StringPair("User-4","4")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_5,StringPair("User-5","5")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_6,StringPair("User-6","6")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_7,StringPair("User-7","7")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_8,StringPair("User-8","8")));
	names_and_symbols.insert(pair<AtomType,StringPair >(USER_9,StringPair("User-9","9")));
}

//helper function - creates map for storing the associated HSV colour for each AtomType.
void AtomicInfo::buildColourMap() {

	int frequency = 3; //number of colours per revolution of wheel
	double delta_hue = 1.0 / (double) frequency;

	//keep track of channels
	double hue = 0, value = 0.5;

	vector<AtomType> types = allAtomTypes();
	uint num_colours, num_types;
	num_colours = num_types = types.size();

	//Create the colours and put in a container, then apply them to the types
	vector<HSVColour> colours;

	//need to make first revolution of colours
	for(int i = 0; i < frequency; ++i) {

		HSVColour cur_colour(hue, 1.0, value);
		colours.push_back(cur_colour);
		this->colours.insert(pair<AtomType,HSVColour>(types[num_types - num_colours], cur_colour));
		hue += delta_hue;
		--num_colours;
	}

	//get colours from recursive helper function
	addColours(num_colours, colours, value, types);
}

//Recursive function called by buildColourMap that does the work of building the list of
void AtomicInfo::addColours(uint& num_colours, vector<HSVColour>& colours, double value, vector<AtomType>& types) {

	//go through 'colours' and add new colours in between the previous colours' hues and using 'value'
	uint index = 1;
	uint init_size = colours.size();
	uint num_types = types.size();

	while(num_colours > 0) {

		//want the hue to be in between, i.e. the average
		double hue = (colours[index].hue() + colours[index-1].hue()) * 0.5;
		HSVColour cur_colour(hue, 1.0, value);
		colours.push_back(cur_colour);
		this->colours.insert(pair<AtomType,HSVColour>(types[num_types - num_colours], cur_colour));

		--num_colours;
		++index;

		//if we go around the wheel, decrease the 'value' and recursively call the current function
		if(index == init_size) {

			//create a hue that loops around and uses the first and last values
			hue = (colours[0].hue() + 1 + colours[index-1].hue()) * 0.5;
			cur_colour = HSVColour(hue, 1.0, value);
			colours.push_back(cur_colour);
			this->colours.insert(pair<AtomType,HSVColour>(types[num_types - num_colours], cur_colour));

			//decrement value based on number of colours to create - never want it above 0.75
			uint total_colours = num_colours + colours.size();
			value += 0.25 / (double) total_colours;

			//go into function again for next revolution of colour wheel, but need to sort first
			sort(colours.begin(), colours.end());
			addColours(num_colours, colours, value, types);
		}
	}

	//returns through all the recursive calls
	return;
}

/*
 * Helper function - creates map for storing the Xray atomic form factor constants for elemental particles
 *
 * Not sure how most of the coefficients were derived, presumably from some fitting of the Hartree-Fock model
 * or possibly something slightly more sophisticated. As a check, the carbon values were used to plot F(q) against the
 * values from the Cromer reference below out to q = 3 Ang^-1 and the curves were almost exactly the same.  The Na+ and
 * Cl- values were thus taken directly from the following:
 *
 * D.T. Cromer and J.B. Mann (1968) X-ray scattering factors computed from numerical Hartree-Fock wave functions. Acta Cryst., A24, 321.
 */
void AtomicInfo::buildAFFconstantMap() {

	//build the vector to insert, then put in the map, do for each element
	vector<double> hydrogen;
	hydrogen.push_back(0.493); //
	hydrogen.push_back(0.323);
	hydrogen.push_back(0.14);
	hydrogen.push_back(0.041);
	hydrogen.push_back(10.511);
	hydrogen.push_back(26.126);
	hydrogen.push_back(3.142);
	hydrogen.push_back(57.8);
	hydrogen.push_back(0.003);

	aff_constants.insert(AtomicFFPair(HYDROGEN, hydrogen));

	vector<double> deuterium;
	deuterium.push_back(0.493); //a1
	deuterium.push_back(0.323); //a2
	deuterium.push_back(0.14); //a3
	deuterium.push_back(0.041); //a4
	deuterium.push_back(10.511); //b1
	deuterium.push_back(26.126); //b2
	deuterium.push_back(3.142); //b3
	deuterium.push_back(57.8); //b4
	deuterium.push_back(0.003); //c

	aff_constants.insert(AtomicFFPair(DEUTERIUM, deuterium));

	vector<double> carbon;
	carbon.push_back(2.31); //a1
	carbon.push_back(1.02); //a2
	carbon.push_back(1.589); //a3
	carbon.push_back(0.865); //a4
	carbon.push_back(20.844); //b1
	carbon.push_back(10.208); //b2
	carbon.push_back(0.569); //b3
	carbon.push_back(51.651); //b4
	carbon.push_back(0.216); //c

	aff_constants.insert(AtomicFFPair(CARBON, carbon));

	vector<double> nitrogen;
	nitrogen.push_back(12.213); //a1
	nitrogen.push_back(3.132); //a2
	nitrogen.push_back(2.013); //a3
	nitrogen.push_back(1.166); //a4
	nitrogen.push_back(0.006); //b1
	nitrogen.push_back(9.893); //b2
	nitrogen.push_back(28.997); //b3
	nitrogen.push_back(0.583); //b4
	nitrogen.push_back(-11.524); //c

	aff_constants.insert(AtomicFFPair(NITROGEN, nitrogen));

	vector<double> oxygen;
	oxygen.push_back(3.049); //a1
	oxygen.push_back(2.287); //a2
	oxygen.push_back(1.546); //a3
	oxygen.push_back(0.867); //a4
	oxygen.push_back(13.277); //b1
	oxygen.push_back(5.701); //b2
	oxygen.push_back(0.324); //b3
	oxygen.push_back(32.909); //b4
	oxygen.push_back(0.251); //c

	aff_constants.insert(AtomicFFPair(OXYGEN, oxygen));

	vector<double> sodium_plus;
	sodium_plus.push_back(3.99479); //a1
	sodium_plus.push_back(3.37245); //a2
	sodium_plus.push_back(1.13877); //a3
	sodium_plus.push_back(0.65118); //a4
	sodium_plus.push_back(3.11047); //b1
	sodium_plus.push_back(7.14318); //b2
	sodium_plus.push_back(0.40692); //b3
	sodium_plus.push_back(15.7319); //b4
	sodium_plus.push_back(0.84267); //c

	aff_constants.insert(AtomicFFPair(SODIUM_PLUS, sodium_plus));

	vector<double> phosphorus;
	phosphorus.push_back(6.435); //a1
	phosphorus.push_back(4.179); //a2
	phosphorus.push_back(1.78); //a3
	phosphorus.push_back(1.491); //a4
	phosphorus.push_back(1.907); //b1
	phosphorus.push_back(27.157); //b2
	phosphorus.push_back(0.526); //b3
	phosphorus.push_back(68.164); //b4
	phosphorus.push_back(1.115); //c

	aff_constants.insert(AtomicFFPair(PHOSPHORUS, phosphorus));

	vector<double> chloride;
	chloride.push_back(18.0842); //a1
	chloride.push_back(7.47202); //a2
	chloride.push_back(6.46337); //a3
	chloride.push_back(2.43918); //a4
	chloride.push_back(0.00129); //b1
	chloride.push_back(1.12976); //b2
	chloride.push_back(19.3079); //b3
	chloride.push_back(59.0633); //b4
	chloride.push_back(-16.4654); //c

	aff_constants.insert(AtomicFFPair(CHLORIDE, chloride));
}
