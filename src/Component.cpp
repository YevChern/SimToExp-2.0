
#include "Component.h"

//Default constructor - not intended for use
Component::Component() {

	//initialize everything to zero
	cmp_name = "";
	number_groups = 0;
	total_electrons = 0;
	total_neutron_SL = 0;
	cmp_volume = 0;
	probs_calculated = false;
}

//Constructor that sets the name and number of groups.  The latter is
//necessary for updating the electron and nsld info
Component::Component(const QString& name, double number_groups) {
   
   cmp_name = name;
   this->number_groups = number_groups;
   
   //initialize values
   total_electrons = 0;
   total_neutron_SL = 0;
   cmp_volume = 0;
   probs_calculated = false;
}

//Destructor
Component::~Component() {}

//Adds a ParticleDensity object (i.e. an 'atom') to the component, and updates
//the electrons and neutrons
void Component::addParticle(const NormalDensity& particle) {
   
   particles.append(particle);
   
   total_electrons += particle.numberOfElectrons();
   total_neutron_SL += particle.neutronSL();
}

//Just sets the component volume
void Component::setVolume(double volume) {cmp_volume = volume;}

//Bunch of getters
QString Component::name() {return cmp_name;}
const QString& Component::name() const {return cmp_name;}
double Component::electrons() const {return total_electrons;}
double Component::neutronSL() const {return total_neutron_SL;}
double Component::numberOfGroups() const {return number_groups;}
double Component::volume() const {return cmp_volume;}
uint Component::numberOfAtoms() const {return particles.size();}

Points& Component::getNumberDensities() {return num_densities;}
const Points& Component::getNumberDensities() const {return num_densities;}
Points& Component::getElectronDensities() {return elec_densities;}
const Points& Component::getElectronDensities() const {return elec_densities;}
Points& Component::getNeutronSLDensities() {return neutron_SL_densities;}
const Points& Component::getNeutronSLDensities() const {return neutron_SL_densities;}

//Returns the probability density function for the component
Points& Component::getProbabilities() {

	//only calculate if not done yet
	if(not probs_calculated)
		calcProbabilities();

	return probabilities;
}

//Calculated values
double Component::atomsPerGroup() const {
   
   return particles.size() / number_groups;
}

double Component::electronsPerGroup() const {
   
   return total_electrons / number_groups;
}

double Component::neutronSLPerGroup() const {
   
   return total_neutron_SL / number_groups;
}

//Function to be called after all particles have been added; calculates all three of the density functions
void Component::calcDensities() {

	//clear all the container in case
	num_densities.clear();
	elec_densities.clear();
	neutron_SL_densities.clear();

	double atoms_per_group = atomsPerGroup(); //just so it isn't calculated numerous times

	//for the 'x' values
	double num_points = particles[0].densities().size(); //NOTE: should put an error check in for empty lists

	//initialize the points
	for(int i = 0; i < num_points; i++) {

		num_densities.append(QPointF ((particles[0].densities())[i].x(), 0));
		elec_densities.append(QPointF ((particles[0].densities())[i].x(), 0));
		neutron_SL_densities.append(QPointF ((particles[0].densities())[i].x(), 0));
	}

	//now go through all of the particles and sum the densities
	for(int i = 0; i < num_points; i++) {
		for(int j = 0; j < particles.size(); j++) {

			num_densities[i].ry() += (particles[j].densities())[i].y();
			elec_densities[i].ry() += (particles[j].electronDensities())[i].y();
			neutron_SL_densities[i].ry() += (particles[j].neutronDensities())[i].y();
		}

		//normalize the number density by the number of atoms per group
		num_densities[i].ry() /= atoms_per_group;
	}

	//probabilities now outdated
	probs_calculated = false;
}

//Private function that calculates the probabilities from the number density
void Component::calcProbabilities() {

	//clear first just in case already done
	probabilities.clear();

	for(int i = 0; i < num_densities.size(); i++) {

		double p = cmp_volume * num_densities[i].y();
		probabilities.append(QPointF (num_densities[i].x(), p));
	}

	probs_calculated = true;
}
