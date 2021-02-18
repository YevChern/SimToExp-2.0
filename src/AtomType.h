//enum for different element types - most of the biologically relevent ones are here, could always just put in the entire periodic
//table just for fun.  Also includes united-atom naming and up to 9 user defined types.

#ifndef ATOMTYPE_H_
#define ATOMTYPE_H_ 

enum AtomType {

	//elements
	HYDROGEN, DEUTERIUM, CARBON, NITROGEN, OXYGEN, FLUORINE, PHOSPHORUS, SULPHUR, SODIUM, MAGNESIUM, CHLORINE, POTASSIUM, CALCIUM,

	//ions
	FLUORIDE, SODIUM_PLUS, MAGNESIUM_TWOPLUS, CHLORIDE, POTASSIUM_PLUS, CALCIUM_TWOPLUS,

	//united atom CG particles
	UA_METHYL, UA_METHYLENE, UA_METHINE, UA_WATER, UA_HEAVY_WATER,

	//martini CG particles
	//Lipids - headgroups and glycerides
	MART_NC3, MART_NH3, MART_PO4, MART_CNO,
	MART_GL0, MART_GL1, MART_GL2, MART_GL3, MART_GL4,

	//Lipids - tails, unfortunately are specific to each tail type, not each particle
	//Palmitoyl (e.g. DPPC) - 16:0 carbon
	MART_PC1, MART_PC2, MART_PC3, MART_PC4,

	//Hexanoyl (e.g. DHPC) - 6:0 carbon
	MART_HC1, MART_HC2,

	//Lauroyl (e.g. DLPC) - 12:0 carbon
	MART_LC1, MART_LC2, MART_LC3,

	//Myristoyl (e.g. DMPC) - 14:0 carbon
	MART_MC1, MART_MC2, MART_MC3,

	//Stearoyl (e.g. DSPC) - 18:0 carbon
	MART_SC1, MART_SC2, MART_SC3, MART_SC4, MART_SC5,

	//Oleoyl (e.g. POPC) - 18:1 carbon
	MART_OC1, MART_OC2, MART_OD3, MART_OC4, MART_OC5,

	//Arachidonoyl - 20:4 carbon
	MART_AD1, MART_AD2, MART_AD3, MART_AD4, MART_AC5,

	//Linoleyl - 18:2 carbon
	MART_UC1, MART_UD2, MART_UD3, MART_UC4,

	//BOLA - cyclic di-DP
	MART_BC1, MART_BC2, MART_BC3, MART_BC4,

	//martini water, heavy water and polarizable types
	MART_WATER, MART_POL, MART_D, MART_D2O,

	//user defined
	USER_1, USER_2, USER_3, USER_4, USER_5, USER_6, USER_7, USER_8, USER_9,

	//for calculating system totals
	SYSTEM_TOTAL,

	//default
	UNKNOWN
};

enum ResidueType {

	//At some point might include proteins?

	//Lipids; NOTE: DL = DM
	DHPC, DPPC, DLPC, DMPC, DSPC, DAPC, DUPC, POPC, DOPC,
	DPPE, DHPE, DLPE, DMPE, DSPE, POPE, DOPE,
	PPCS, //Sphingolipids
	DOPG, POPG, DOPS, POPS, //Negatively charged
	BOLA, BOLB, //cyclic and acyclic di-DPPC respectively

	WATER, HEAVY_WATER, ION,

	//default
	UNKNOWN_RESIDUE
};

#endif /* ATOMTYPE_H_ */
