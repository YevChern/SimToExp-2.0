/************************************************************************/
//	Simple program for the volume analysis of simulated profiles		//
//                                                                      //
//  v2.0 has switched to using Qt as a GUI library over Tcl/Tk/BLT      //
//                                                                      //
//  - original code and analysis routines by Norbert Kucerka at NRC     //
//  - 'Qt' code and trajectory reader added by Bryan W. Holland at      //
//     University of Calgary     									    //
//                                                                      //
//***********************************************************************/
// If you use this program, you MUST cite:                              //
// N. Kucerka, J. Katsaras and J.F. Nagle (2010) Comparing Membrane     //
// Simulations to Scattering Experiments: Introducing the SIMtoEXP      //
// Software. J. Membr. Biol., 235(1): 43-50                             //
//																		//
//  It is based on Horia's original routine written in fortran			//
//c Reference:															//
//c Determination of component volumes of lipid bilayers				//
//c from simulations. 1997.												//
//c H. I. Petrache, S. E. Feller, and J. F. Nagle						//
//c Biophys. J., 72, 2237--2242									        //
/************************************************************************/
//								****									//
//		****		Norbert.Kucerka@nrc-cnrc.gc.ca		****			//
//                     bryan.holland@ucalgary.ca                        //
//								****									//
/************************************************************************/

#include "STEwindow.h"

using namespace std;

//Main function... just call GUI
int main(int argc, char* argv[]){

	//create a QApplication (necessary) and QMainWindow object
	QApplication app(argc, argv);

	//populate the window with the code from Qt Designer
	STEwindow* window = new STEwindow;
	window->show();

	//run the application
	return app.exec();
}
