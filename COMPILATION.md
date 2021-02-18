-----------------------
Compilation and running
-----------------------

To compile SimToExp you will need the following two libraries installed:
- Qt a popular C++ GUI library (likely already installed on most boxes) .
- Qwt - a plotting add-on for Qt .

You can install Qwt just by going to the terminal and type:
```
sudo apt-get install libqwt-dev
````
In some repositories qwt may go by name libqwt-qt5-dev, so if libqwt-dev package is not found, just substitute the libqwt-dev name in the command above with libqwt-qt5-dev.
After the installation of Qwt, go to the directory where you unpacked your SimToExp source code and just type:
```
qmake SIMtoEXP.pro
make
```
This should be it, although if you don't have qmake you can install it with:
```
sudo apt-get install qt5-qmake
```
Also, it might need to change Qwt lib name in SIMtoEXP.pro file to be properly detected. If Qwt is not automatically detected during compilation, try changing the following line in the SIMtoEXP.pro file:
```
LIBS +=  -lgomp -lqwt
```
to
```
LIBS +=  -lgomp -lqwt-qt5
```
NOTE: SimToExp can be compiled with either Qt4 or Qt5. Just make sure that you’re using the Qwt version that is compatible with the Qt version of your choice. You can obtain source code of Qwt version <6.1 and f
```
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:<path to Qwt>/lib"
```

At the moment, no attempt has been made to compile on either a Windows or Mac OS, so if anyone would like to try, please let me know how it goes. This may involve a bit of work, since Qt uses 'qmake' to generate

Windows: http://qt-project.org/wiki/Support_for_Windows

Mac: http://qt-project.org/wiki/Support_for_Mac_OS_X

