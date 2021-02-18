/*
 * TrajectoryInfoWindow.h
 *
 *  Created on: Apr 12, 2013
 *      Author: bholland
 *
 *  A simple window class that acquires the information necessary to open, read and analyze an input
 *  trajectory.  Since all trajectories types (e.g. Gromacs, NAMD, AMBER, etc.) require the same info,
 *  this class can be used for all of them.
 *
 *  Is not allowed to exist without a parent, as it is intended to be called by the main window
 */

#ifndef TRAJECTORYINFOWINDOW_H_
#define TRAJECTORYINFOWINDOW_H_

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QWidget>
    #include <QtWidgets/QGroupBox>
    #include <QtWidgets/QLineEdit>
    #include <QtWidgets/QPushButton>
    #include <QtWidgets/QRadioButton>
    #include <QtWidgets/QLabel>
    #include <QtWidgets/QComboBox>
#else
    #include <QtGui/QWidget>
    #include <QtGui/QGroupBox>
    #include <QtGui/QLineEdit>
    #include <QtGui/QPushButton>
    #include <QtGui/QRadioButton>
    #include <QtGui/QLabel>
    #include <QtGui/QComboBox>
#endif

#include <QtGui/QCloseEvent>

#include "TrajectoryInfo.h"

class TrajectoryInfoWindow : public QWidget {

	//qt macro
	Q_OBJECT;

	TrajectoryInfo info;

	//Group boxes for the different regions of the window
	QGroupBox path_box;
	QGroupBox granularity_box;
	QGroupBox axis_box;
	QGroupBox heavy_fraction_box;
	QGroupBox normal_box;
	QGroupBox planar_box;

	//Parts of the boxes needed for information retrieval
	QComboBox traj_chooser;
	QLineEdit path_line_edit;
	QPushButton browse_button;
	QLabel opt_path_label;
	QLineEdit opt_path_line_edit;
	QPushButton opt_browse_button;

	QRadioButton atomistic_button;
	QRadioButton cg_button;

	QRadioButton x_button;
	QRadioButton y_button;
	QRadioButton z_button;

	QLabel heavy_fraction_label;
	QLineEdit heavy_fraction_lineEdit;

	QLabel normal_grid_label;
	QLineEdit normal_grid_lineEdit;
	QRadioButton normal_avg_button;
	QRadioButton normal_time_button;
	QLabel normal_num_label;
	QLineEdit normal_num_lineEdit;

	QLabel radial_grid_label;
	QLineEdit radial_grid_lineEdit;
	QLabel angular_grid_label;
	QLineEdit angular_grid_lineEdit;
	QRadioButton planar_avg_button;
	QRadioButton planar_time_button;
	QLabel planar_num_label;
	QLineEdit planar_num_lineEdit;

	//Command buttons
	QPushButton done_button;
	QPushButton cancel_button;

public:

	TrajectoryInfoWindow(QWidget*);
	virtual ~TrajectoryInfoWindow();

	TrajectoryInfo trajectoryInfo();
	const TrajectoryInfo& trajectoryInfo() const;

	void closeEvent(QCloseEvent*);

private slots:

	void updateFilePath(const QString&);
	void updateOptFilePath(const QString&);

	void openFileBrowser();
	void updateTrajInfo(int);

	void normalLineEditWritable(bool);
	void planarLineEditWritable(bool);

	void calcNormal(bool);
	void calcPlanar(bool);

	void dataEntryComplete();
	void cancel();

signals:

	void dataEntered(TrajectoryInfoWindow*);
	void canceled(TrajectoryInfoWindow*);

private:

	void setupWindow();
	void createPathBox();
	void createGranularityBox();
	void createAxisBox();
	void createHeavyFractionBox();
	void createNormalBox();
	void createPlanarBox();
};

#endif /* TRAJECTORYINFOWINDOW_H_ */
