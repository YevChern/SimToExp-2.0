
#include "ExpData.h"

#include <cmath>
#include <cfloat>

#include <QtCore/QDebug>

const int DIGITS_DISPLAYED = 5;
const double INITIAL_SCALING_DELTA = 0.1; //scaling factor is usually pretty close
const double CHI_SQUARED_TOLERANCE =  1E-8; //might be overkill, but it is very fast

//Default constructor
ExpData::ExpData(QWidget* parent) : QWidget(parent) {
   
   scaling_factor = 1.0;
   chi_squared = 0.0;
   autoscaled = false;
   aligned = false;
   use_uncertainty = true;
}

//Destructor
ExpData::~ExpData() {}

//Add a point to the collection
void ExpData::addPoint(double x, double y, double err) {
   
   data.push_back(QPointF(x,y));
   uncertainty.push_back(err);
}

//Getters
double ExpData::getScalingFactor() const {return scaling_factor;}
double ExpData::getChiSquared() const {return chi_squared;}

//Return the data
Points& ExpData::getData() {return data;}
const Points& ExpData::getData() const {return data;}

//Return the data aligned with any provided simulation data
Points& ExpData::getAlignedData() {return data_aligned;}
const Points& ExpData::getAlignedData() const {return data_aligned;}

//Uses equation 2 from SIMtoEXP paper to calculate the automatic scaling factor.
//Also scales the data using the scaling factor, both for convenience and to
//guarantee the values are always scaled to the current scaling factor
Points& ExpData::calcScalingFactor(const Points& sim_data) {
   
   //if scaling already performed, don't do again
   if(autoscaled)
      return data;

   //need to align the data in order to scale properly
   alignData(sim_data);

   //check to see if there is any uncertainty present at all
   double uncert_sum = 0;
   for(double uncert : uncertainty_aligned)
	   uncert_sum += uncert;

   //if the uncertainty is ~ 0, don't use it
   if(uncert_sum < 1E-14)
	   use_uncertainty = false;

   //calculate the denominator
   double denominator = 0, numerator = 0;
   
   int data_size = data_aligned.size();
   for(int i = 0; i < data_size; ++i) {

	   double exp_val = data_aligned[i].y();
	   double sim_val = sim_aligned[i].y();
	   double uncert_squared = pow(uncertainty_aligned[i],2);

	   //weighted by the uncertainty, but if all uncertainties are zero ignore this part
	   if(use_uncertainty) {

		   denominator += pow(exp_val,2) / uncert_squared;
		   numerator += fabs(exp_val * sim_val) / uncert_squared;

	   } else {

		   denominator += pow(exp_val,2);
		   numerator += fabs(exp_val * sim_val);
	   }
   }
   
   scaling_factor = numerator / denominator;
   
   //actually scale the data - need to scale the aligned and original data set
   for(int i = 0; i < data_size; ++i) {

      data_aligned[i].ry() *= scaling_factor;
      uncertainty_aligned[i] *= scaling_factor;
   }
   
   for(int i = 0; i < data.size(); ++i) {

      data[i].ry() *= scaling_factor;
      uncertainty[i] *= scaling_factor;
   }

   //now calculate chi_squared
   numerator = 0;
      
   for(int i = 0; i < data_size; ++i) {

	   if(use_uncertainty)
		   numerator += pow(fabs(sim_aligned[i].y()) - fabs(data_aligned[i].y()),2) / pow(uncertainty_aligned[i],2);
	   else
		   numerator += pow(fabs(sim_aligned[i].y()) - fabs(data_aligned[i].y()),2);
   }
   
   chi_squared = numerator / (data_size - 1);

   //now with a starting guess at the best fit, minimize chi squared numerically
   minimizeChiSquared(INITIAL_SCALING_DELTA);
   
   //Create parsed strings for displaying
   QString scale_string = QString::number(scaling_factor);
   scale_string.truncate(DIGITS_DISPLAYED);
   
   QString chi_string = QString::number(chi_squared);
   chi_string.truncate(DIGITS_DISPLAYED);
   
   emit scalingFactorChanged(scale_string);
   emit chiSquaredChanged(chi_string);
   autoscaled = true;
   
   return data;
}

//Scales the data by 'factor' of the ORIGINAL values. Returns the scaled values.  'sim_data' is
//optional, the default is a NULL pointer, so is the 'emitChange' boolean, the default is true.
Points& ExpData::scaleData(double factor, bool emitChange) {
   
   //current scaling is 'scaling_factor', so to get correct scaling of
   //original values, simply multiply by the ratio of the current over the previous
	double cur_scaling = factor / scaling_factor;

	int data_size = data.size();
	for(int i = 0; i < data_size; ++i) {

		data[i].ry() *= cur_scaling;
		uncertainty[i] *= cur_scaling;
	}

	//if aligned, also rescale the aligned data
	if(aligned) {

		data_size = data_aligned.size();
		for(int i = 0; i < data_size; ++i) {

			data_aligned[i].ry() *= cur_scaling;
			uncertainty_aligned[i] *= cur_scaling;
		}
	}

	//finally set the scaling factor
	scaling_factor = factor; //don't emit change, since this was provided by the user in first place

	//Now that the data is scaled, if aligned recalculate chi_squared
	if(aligned) {

		double numerator = 0;

		data_size = data_aligned.size();
		for(int i = 0; i < data_size; ++i) {

			if(use_uncertainty)
				numerator += pow(fabs(sim_aligned[i].y()) - fabs(data_aligned[i].y()),2) / pow(uncertainty_aligned[i],2);
			else
				numerator += pow(fabs(sim_aligned[i].y()) - fabs(data_aligned[i].y()),2);
		}

		chi_squared = numerator / (data_size - 1);

		if(emitChange) {

			QString chi_string = QString::number(chi_squared);
			chi_string.truncate(DIGITS_DISPLAYED);

			emit chiSquaredChanged(chi_string);
		}
	}
   
   return data;
}

//Private function that goes through the experimental data and creates a set of simulation data
//that corresponds to it - needed for proper calculation of the scaling factor
void ExpData::alignData(const Points& sim_data) {

	//make sure the storage containers are empty
	data_aligned.clear();
	uncertainty_aligned.clear();
	sim_aligned.clear();

	int sim_cursor = 0;

	//go through all experimental data
	int data_size = data.size();
	int sim_size = sim_data.size();

	//flag for whether the data sim data has been reached
	bool at_sim_data = false;

	for(int i = 0; i < data_size; ++i) {

		//skip ahead if not at the sim data yet (this can always be avoided in the program by
		//setting the proper q values in the Fourier transform).  Thanks to flag, will only skip ahead
		//until it reaches the first sim data, but not after
		if(data[i].x() < sim_data[sim_cursor].x() and not at_sim_data)
			continue;
		else if(not at_sim_data)
			at_sim_data = true;

		//loop through the sim data until similar q
		while(sim_cursor < sim_size and sim_data[sim_cursor].x() < data[i].x()) {

			++sim_cursor;
		}

		//if cursor at end of sim values, nothing left to do
		if(sim_cursor >= sim_size)
			break;

		//otherwise approximate sim value linearly and add both points to the containers
		if(sim_cursor > 0) {

			double m = (sim_data[sim_cursor].y() - sim_data[sim_cursor-1].y()) / (sim_data[sim_cursor].x() - sim_data[sim_cursor-1].x()); //calc slope
			double b = sim_data[sim_cursor].y() - m * sim_data[sim_cursor].x(); //calc y-intercept

			//approx sim value
			double sim_y = m * data[i].x() + b;

			//store values
			data_aligned.append(data[i]);
			uncertainty_aligned.append(uncertainty[i]);
			sim_aligned.append(QPointF(data[i].x(),sim_y));
		}
	}

	//set flag
	aligned = true;
}

//Private, iterative function that minimizes chi_squared using Newton's method. Used after
//the primary scaling factor has been calculated.
void ExpData::minimizeChiSquared(double scaling_delta) {
   
   double prev_chi = chi_squared;
   
   //adjust 'scaling_factor' by delta, but don't update member, it is done in scaleData(),
   //doing so here would mess it up
   double cur_scale = scaling_factor + scaling_delta;
   
   //scale data, but don't emit change
   scaleData(cur_scale, false);
   
   //if delta(chi) < tolerance, minimization is done
   if(fabs(chi_squared - prev_chi) < CHI_SQUARED_TOLERANCE)
      return;
   
   /*
    * If there is uncertainty, want to get a close to one as possible, this is the best possible fit without OVERfitting.
    * With no uncertainty, want to get as close to zero as possible.
    */
   if(use_uncertainty) {

	   //if current chi > prev_chi AND it is greater than one, or if it is less than one, want to move upward toward one!
	   //halve delta and change direction
	   if((chi_squared > 1 and chi_squared > prev_chi) or (chi_squared < 1 and chi_squared < prev_chi))
		   scaling_delta *= -0.5;

   } else {

	   //since want to go to zero, only change direction if chi_squared gets bigger
	   if(chi_squared > prev_chi)
		   scaling_delta *= -0.5;
   }

   //minimize again
   minimizeChiSquared(scaling_delta);
}

//Returns the data to the original values and resets scaling factor to 1 and chi squared to zero
void ExpData::resetData() {

   double scaling_inverse = 1 / scaling_factor; //to avoid multiple divisions
   
   for(int i = 0; i < data.size(); i++) {
      data[i].ry() *= scaling_inverse;
      uncertainty[i] *= scaling_inverse;
   }
   
   //get rid of all alignment data
   data_aligned.clear();
   uncertainty_aligned.clear();
   sim_aligned.clear();
   aligned = false;

   scaling_factor = 1;
   chi_squared = 0;
   autoscaled = false;

   //Create parsed strings for displaying
   QString scale_string = QString::number(scaling_factor);
   scale_string.truncate(DIGITS_DISPLAYED);

   QString chi_string = QString::number(chi_squared);
   chi_string.truncate(DIGITS_DISPLAYED);

   emit scalingFactorChanged(scale_string);
   emit chiSquaredChanged(chi_string);
}

//Returns the smallest value on the x-axis from this set of data
double ExpData::leftMargin() const {

	if(data.size() > 0) {

		double smallest = DBL_MAX;

		//go through and find the smallest x-value, make no assumptions about its location in the container
		for(int i = 0; i < data.size(); i++) {

			if(smallest > data[i].x())
				smallest = data[i].x();
		}

		return smallest;
	}

	//otherwise just return zero
	return 0;
}

//Returns the largest value on the x-axis from this set of data
double ExpData::rightMargin() const {

	if(data.size() > 0) {

		double largest = DBL_MIN;

		//go through and find the smallest x-value, make no assumptions about its location in the container
		for(int i = 0; i < data.size(); i++) {

			if(largest < data[i].x())
				largest = data[i].x();
		}

		return largest;
	}

	//otherwise just return zero
	return 0;
}
