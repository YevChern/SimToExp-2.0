/*
 * STEwindow.h
 *
 *  Created on: Oct 1, 2012
 *      Author: bholland
 * 
 * Main window for SIMtoEXP (STE)
 */

#ifndef STEWINDOW_H_
#define STEWINDOW_H_

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/qmainwindow.h>
#else
    #include <QtGui/qmainwindow.h>
#endif


#include "gui.h"
#include "SimManager.h"
#include "AtomicInfo.h"
#include "ExpManager.h"
#include "RGB.h"
#include "Tab.h"
#include "TrajectoryInfoWindow.h"
#include "OutputWriter.h"

#if QT_VERSION >= 0x050000
    #include <QtWidgets/QFileDialog>
#else
    #include <QtGui/QFileDialog>
#endif

#include <QtCore/QLocale>
#include <qwt/qwt_plot_curve.h>
#include <QtGui/QColor>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_legend.h>
#include <qwt/qwt_legend_label.h>

class STEwindow: public QMainWindow, private Ui::MainWindow {

   //qt macro
   Q_OBJECT

   //useful bools
   bool sim_open;
   bool sim_from_trajectory;
   bool form_factors_done;
   bool gro_file_open;

   //for opening new tabs
   bool scattering_tab_open;
   bool atomicff_tab_open;
   bool sim_specifics_tab_open;
   bool tab_area_maximized;

   //bounding rectangles for all the plots
   QMap<int,QRectF> bounding_rects;
   QMap<int,ScrollZoomer*>  scroll_zoomers;//the zoomer objects
   
   //legends for each plot
   QwtLegend num_densities_legend;
   QwtLegend e_densities_legend;
   QwtLegend neutron_SL_legend;
   QwtLegend xray_formfactor_legend;
   QwtLegend neutron_formfactor_legend;
   QwtLegend volume_legend;
   QMap<int,QwtPlot*> plots; //useful to have plots in mapped container

   //useful for converting numbers to strings - this can be used for changes to other
   //languages / number systems
   QLocale* converter;
   
   //For showing/hiding data deletion menu actions
   QList<QAction*> xray_deletes;
   QList<QAction*> neutron_deletes;

   //for keeping track of dynamically created objects (i.e. for deletion)
   QList<QwtPlotCurve*> sim_number_curves;
   QList<QwtPlotCurve*> sim_ed_curves;
   QList<QwtPlotCurve*> sim_nsld_curves;
   QwtPlotCurve* sim_FF_xray_curve;
   QwtPlotCurve* sim_FF_neutron_curve;
   QList<QwtPlotCurve*> xray_curves;
   QList<QwtPlotCurve*> neutron_curves;
   QList<QwtSymbol*> xray_symbols;
   QList<QwtSymbol*> neutron_symbols;
   QList<QwtPlotCurve*> exp_ed_curves;
   QList<QwtPlotCurve*> exp_nsld_curves;
   QList<QwtPlotCurve*> comp_number_curves;
   QList<QwtPlotCurve*> comp_ed_curves;
   QList<QwtPlotCurve*> comp_nsld_curves;
   QList<QwtPlotCurve*> volume_curves;
   QList<QwtLegendLabel*> sim_number_legend_items;

   //for extra tabs created for tool, etc.
   QList<Tab*> tabs;
      
   //useful for AtomType related info
   AtomicInfo info;
   
   //does much of the work
   SimManager sim;
   ExpManager exp;

   //for plotting
   QList<RGB> plot_colours;
   QList<RGB>::iterator xray_ff_colour_itr;
   QList<RGB>::iterator neutron_ff_colour_itr;
   QList<QwtSymbol::Style> plot_symbols;
   QList<QwtSymbol::Style>::iterator xray_ff_symbol_itr;
   QList<QwtSymbol::Style>::iterator neutron_ff_symbol_itr;
   
   //for connections
   QList<QLineEdit*> xray_scale;
   QList<QLineEdit*> xray_chi;
   QList<QLineEdit*> neutron_scale;
   QList<QLineEdit*> neutron_chi;
   
   //for writing output (duh)
   OutputWriter out_writer;

   public:

      STEwindow(QWidget* parent = 0);
      virtual ~STEwindow();
      
   public slots:

      void openSim(); //only reads in the data
      void displaySim(); //only displays available data
      void displayFourier();//calculates and displays the fourier transform of the input data
      void setFourierDefaults(); //simply resets the LineEdits to the default values
      void deleteSim();
      
      void openTrajectoryFile();
      void cancelTrajectoryOpen(TrajectoryInfoWindow*);
      void readTrajectoryFile(TrajectoryInfoWindow*);
      void readGromacsXTCFile(TrajectoryInfo);

      void openXray();
      void displayXray(unsigned int);
      void openNeutron();
      void displayNeutron(unsigned int);
      
      void openElecDensity();
      void displayElecDensity();
      void openNeutDensity();
      void displayNeutDensity();
      
      void openComponents();
      void displayComponents();
      void displayVolumes();
      
      void replotXrayFormFactor();
      void replotNeutronFormFactor();

      void scaleXrayUser();
      void scaleNeutronUser();
      
      void calcComponentVolumes();
      
      void writeElementalNumberDensities();
      void writeParticleNumberDensities();
      void writeComponentNumberDensities();
      void writeElementalElectronDensities();
      void writeComponentElectronDensities();
      void writeElementalNeutronDensities();
      void writeComponentNeutronDensities();
      void writeXrayFactors();
      void writeNeutronFactors();
      void writeVolumetricFit();
      void writeEverything();

      //can't figure out better way to do this, so whole bunch of slots just pass along the work to
      //private functions
      void deleteXraySample1();
      void deleteXraySample2();
      void deleteXraySample3();
      void deleteXraySample4();
      void deleteXraySample5();
      void deleteNeutronSample1();
      void deleteNeutronSample2();
      void deleteNeutronSample3();
      void deleteNeutronSample4();
      void deleteNeutronSample5();
      
      //for experimental scattering densities
      void deleteElecDensity();
      void deleteNeutDensity();
      
      void deleteComponents();
      
      //for setting visibility of a plot item - for QwtLegends
#if QT_VERSION >= 0x050000
      void showEDensities(const QVariant &itemInfo, bool on, int);
      void showNeutronSLDensities(const QVariant &itemInfo, bool on, int);
      void showXRayFF(const QVariant &itemInfo, bool on, int);
      void showNeutronFF(const QVariant &itemInfo, bool on, int);
      void showNumDensities(const QVariant &itemInfo, bool on, int);
      void showVolProbs(const QVariant &itemInfo, bool on, int);
#else
      void setVisibility(QwtPlotItem*,bool);
#endif

      //plot slots
      void setZoomToBounds();
      void togglePlotArea();

      //tab slots
      void closeTab(int);

      //for tools
      void displayScatteringDensities();

      //does what it says
      void setSimOpenTrue();

      //for debugging only
      void testButton();
      
   signals:
      
      void simOpened();
      void xrayOpened(unsigned int);
      void neutronOpened(unsigned int);
      void elecDensityOpened();
      void neutDensityOpened();
      void componentsOpened();

   private:

      //for reading trajectory data
      void readGromacsXTCFile();

      void makeConnections();

      //removes simulation data
      void removeSim();
      void removeComponents();
      void deleteXrayData(uint);
      void deleteNeutronData(uint);
      
      void updatePlot(int);
      QRectF setBoundingRect(int);
      QRectF setElectronDensityBoundingRect();
      QRectF setNeutronSLDensityBoundingRect();
      QRectF setXrayFormFactorBoundingRect();
      QRectF setNeutronFormFactorBoundingRect();
      QRectF setNumDensityBoundingRect();
      QRectF setVolumeBoundingRect();
      QRectF addMargins(const QRectF&);
      
      QColor getColour(int); //gets a random colour
      QColor getNextXrayColour();
      QColor getNextNeutronColour();
      QwtSymbol::Style getNextXraySymbol();
      QwtSymbol::Style getNextNeutronSymbol();
      
      void updateXrayCurves();
      void updateNeutronCurves();
      void cleanFourierCurves();
      
      void buildColourLibrary();
      void buildSymbolLibrary();
      void initializeBoundingRects();
      
      void buildLineEditContainers();
      void buildDeleteActionsContainers();
      void buildPlotContainers();
      void buildScrollZoomerContainer();

      void resetTabBoolean(TabName);
      void showTab(QWidget*);
      void hideTab(QWidget*);
      bool tabShown(QWidget*);
      bool tabHidden(QWidget*);
      int tabIndex(QWidget*);
      int currentTabIndex();
      QWidget* currentTab();
      void setCurrentTab(QWidget*);
      void goToFirstShownTab();
};

#endif /* STEWINDOW_H_ */
