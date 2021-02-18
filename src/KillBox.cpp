
#include "KillBox.h"

//constructor that takes in a message. 'parent' should be the MainWindow only, so that
//the 'close()' function will kill the program
KillBox::KillBox(QWidget* parent, QString message) {
   
   connect(this, SIGNAL(buttonClicked(QAbstractButton*)), parent, SLOT(close()));
   
   //set up an angry looking box with the message
   setIcon(QMessageBox::Critical);
   setText(message);
   exec();
}

//Destructor
KillBox::~KillBox() {}
