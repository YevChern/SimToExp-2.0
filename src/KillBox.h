/*
 * QKillBox.h
 *
 *  Created on: Oct 10, 2012
 *      Author: bholland
 * 
 * Useful message box that allows a message to be displayed prior to exiting the program
 * due to a critical error
 * 
 */

#ifndef QKILLBOX_H_
#define QKILLBOX_H_

#include <QtGlobal>
#if QT_VERSION >= 0x050000
    #include <QtWidgets/QMessageBox>
#else
    #include <QtGui/QMessageBox>
#endif
#include <QtCore/QString>

class KillBox: public QMessageBox {

   //qt macro
   Q_OBJECT

   public:

      KillBox(QWidget* parent = 0, QString = "");
      virtual ~KillBox();
      
};

#endif /* QKILLBOX_H_ */
