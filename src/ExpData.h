/*
 * ExpData.h
 *
 *  Created on: Oct 19, 2012
 *      Author: bholland
 * 
 * Stores the data from a particular experimental file - same format for X-ray and neutron data.  Would
 * prefer not to use the Qt classes here, but it makes it much easier to pass them along later flor plotting
 * 
 * Currently does NOT store the error values from the experimental data
 */

#ifndef EXPDATA_H_
#define EXPDATA_H_

#include <QtCore/QVector>
#include <QtCore/QPointF>
#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QWidget>
#else
    #include <QtGui/QWidget>
#endif

typedef QVector<QPointF> Points;
typedef QVector<double> Doubles;

class ExpData : public QWidget {
   
   Q_OBJECT
   
   Points data;
   Doubles uncertainty;
   
   //For storing the points after having the simulation results aligned to the data
   //Needed for proper calculation of the scaling factor
   Points data_aligned;
   Doubles uncertainty_aligned;
   Points sim_aligned;

   double scaling_factor;
   double chi_squared;

   bool autoscaled;
   bool aligned;
   bool use_uncertainty;
   
   public:
   
      ExpData(QWidget* parent = 0);
      virtual ~ExpData();
      
      double getScalingFactor() const;
      double getChiSquared() const;
      
      void addPoint(double,double,double);  //only add for now
      
      Points& getData();
      const Points& getData() const;
      Points& getAlignedData();
      const Points& getAlignedData() const;
      Points& calcScalingFactor(const Points&);
      Points& scaleData(double,bool emitChange = true);
      
      //currently search through entire set every time - should consider sorting the list
      //upon creation, thereby guaranteeing the min and max margins will be at the extremities
      //fast enough not to care though
      double leftMargin() const;
      double rightMargin() const;

      void resetData();
      
   signals:
      
      void scalingFactorChanged(QString);
      void chiSquaredChanged(QString);
      
   private:
   
      void alignData(const Points&);
      void minimizeChiSquared(double);
};

#endif /* EXPDATA_H_ */
