
#include <iostream>
#include <cmath>
#include <omp.h>

#include "FourierFormFactor.h"

using namespace std;

const int NUM_THREADS = 12;
const int MASTER_THREAD_ID = 0;

//Default constructor
FourierFormFactor::FourierFormFactor(ScatteringType type, AtomicInfo& info) {

	//set default values
	q_min = 0;
	q_max = 0;
	q_step = 1; //so won't divide by zero

	solvent_sd = 0;

	ff_type = type;
	this->info = &info;

	symmetrized = false;
}

//Destructor
FourierFormFactor::~FourierFormFactor() {}

//getters
double FourierFormFactor::qMin() const {return q_min;}
double FourierFormFactor::qMax() const {return q_max;}

//setters
void FourierFormFactor::setSymmetrized(bool value) {

	symmetrized = value;
}

//set the three constant Q values and calculates the q points
void FourierFormFactor::setQvalues(double q_min, double q_max, double q_step) {

	this->q_min = q_min;
	this->q_max = q_max;
	this->q_step = q_step;

	calcPoints();
}

//setter for the scattering density
void FourierFormFactor::setSolventScatteringDensity(double sd) {

	solvent_sd = sd;
}

//setter for the real space axis of points. Goes through the QPointF objects and keep
//the x values
void FourierFormFactor::setRealSpaceValues(const Points& points) {

	int size = points.size();
	for(int i = 0; i < size; i++) {

		x.append(points[i].x());
	}
}

Points& FourierFormFactor::getFormFactors() {return form_factors;}
const Points& FourierFormFactor::getFormFactors() const {return form_factors;}

//helper function for calculating 'q' points
void FourierFormFactor::calcPoints() {

	//figure out the number of q values to use
	//don't want to go over q_max (no data to integrate), so use floor function
	int N_q = floor((q_max - q_min)/q_step);

	//set all the points
	for (int i = 0; i <= N_q; ++i)
		q.push_back(q_min + i*q_step);
}

//creates the atomic form factors if the ff_type is set to XRAY, otherwise does nothing
//NOTE: since neutron AFFs are simply the neutron SL density, these values should be obtained
//directly from AtomicInfo
void FourierFormFactor::createAtomicFormFactors(const AtomTypes& types) {

	if(ff_type == XRAY) {

		int size = types.size();
		for(int i = 0; i < size; ++i) {

			XrayAtomicFormFactor xaff(types[i], *info);
			xray_atomic_FFs.insert(types[i],xaff);
		}
	}
}

//Calculates the form factors assuming that everything has been set already
void FourierFormFactor::calcFormFactors(NormalDensities densities, QProgressDialog& progress, uint init_progress) {

	//first if for xray, go through and calculate all of the xray atomic form factors and store them for use later
	AtomicFFMap xray_affs;
	if(ff_type == XRAY) {

		int pds_size = densities.size();
		for(int k = 0; k < pds_size; ++k) {

			AtomType cur_Type = densities[k].atomType();
			vector<double> cur_FFs;

			//only calculate and insert if the AtomType hasn't been done yet
			if(xray_affs.find(cur_Type) == xray_affs.end()) {

				AtomTypeXrayMap::iterator itr = xray_atomic_FFs.find(cur_Type);

				//if not found, output a debug error
				if(itr == xray_atomic_FFs.end()) {

					cerr << "FourierFormFactor::calcFormFactors: type not found in atomic FF map" << endl;
					return;
				}

				int q_size = q.size();
				for(int i = 0; i < q_size; ++i)
					cur_FFs.push_back(itr->calcAtomicFF(q[i]));

                itr->setAtomicFF_calculated(true);
				xray_affs.insert(AtomicFFPair(cur_Type, cur_FFs));
			}
		}
	}

	//For each 'q' point, do the integral over the bounds given.  The integral uses the trapezoid rule.

	//use OpenMP to break up the work among threads
#pragma omp parallel default(none) shared(xray_affs, densities, progress, init_progress) num_threads(NUM_THREADS)
	{
		int q_size = q.size();

#pragma omp for schedule(dynamic)
	for (int i = 0; i < q_size; ++i) {

		//if main thread, update the progress bar
		if(omp_get_thread_num() == MASTER_THREAD_ID) {

			progress.setValue(i + init_progress);
		}

		double Re = 0, Im = 0;
		double cur_Re = 0, prev_Re = 0, cur_Im = 0, prev_Im = 0;
		double dx = 0;

		//go through all points in real space
		int x_size = x.size();
		for(int j = 0; j < x_size; ++j) {

			double total_sd = 0; //total scattering density

			//sum according to whether neutron or xray, also taking into account symmetrized
			int pds_size = densities.size();
			for(int k = 0; k < pds_size; ++k) {

				AtomType cur_Type = densities[k].atomType();
                double cur_FF;

				if(ff_type == XRAY)
					cur_FF = xray_affs[cur_Type].at(i);

				else if(ff_type == NEUTRON)
					cur_FF = info->neutronSL(cur_Type);

				else
					cur_FF = 0; //so compiler doesn't complain

				//if not symmetrized, use original data, otherwise symmetric densities
				if(not symmetrized)
					total_sd += cur_FF * densities[k].densities().at(j).y(); //should be exact same x-axis as 'x'
				else
					total_sd += cur_FF * densities[k].symmetricDensities().at(j).y();

			}

			cur_Re = (total_sd - solvent_sd) * cos(q[i]*x[j]);

			if(not symmetrized)
				cur_Im = (total_sd - solvent_sd) * sin(q[i]*x[j]);

			//average the integrands for the trapezoid rule
			if(j > 0) {

				dx = x[j] - x[j-1];

				Re += (cur_Re + prev_Re) * dx / 2;
				Im += (cur_Im + prev_Im) * dx / 2;
			}

			//set up for next loop in real space
			prev_Re = cur_Re;
			prev_Im = cur_Im;
		}

		//integration for 'q' finished, store as a QPointF object
		QPointF point(q[i],sqrt(Re*Re + Im*Im));

		//Need to put a mutex lock around the container for appends
		ff_mutex.lock();
		form_factors.append(point);
		ff_mutex.unlock();
	}
}

	//end of parallel region - sort the form factor points so they display properly
	sortPoints();
}

//Sorts the points by their 'q' values
void FourierFormFactor::sortPoints() {

	//do a heap sort - slower, but easy to write and too many function calls in a recursive mergeSort

	//dump points into working container
	Points to_sort = form_factors;
	form_factors.clear();

	//special case for i = 0
	form_factors.append(to_sort[0]);

	int size = to_sort.size();
	for(int i = 1; i < size; i++) {

		bool inserted = false;

		//go through form_factors to find correct placement
		for(int j = 0; j < form_factors.size(); j++) {

			if(to_sort[i].x() < form_factors[j].x()) {

				form_factors.insert(j, to_sort[i]);
				inserted = true;
				break;
			}
		}

		//if bigger than everything, just put at the end
		if(not inserted)
			form_factors.append(to_sort[i]);
	}

	//and we're done
}

//simply returns the number of 'q' values in the container
uint FourierFormFactor::getNumberQValues() const {return q.size();}

//cleans up everything, i.e. cleans out all containers
void FourierFormFactor::clear() {

	x.clear();
	q.clear();
	xray_atomic_FFs.clear();
	form_factors.clear();
}
