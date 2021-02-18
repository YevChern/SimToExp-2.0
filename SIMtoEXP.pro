DESTDIR = .
OBJECTS_DIR = build
MOC_DIR = src/moc

TEMPLATE = app
TARGET = SIMtoEXP
CONFIG += release
DEPENDPATH += .
INCLUDEPATH += ./src
LIBS +=  -lgomp -lqwt

QMAKE_CXXFLAGS += -std=c++0x -fopenmp

greaterThan(QT_MAJOR_VERSION, 4) {
    message( $$[QT_VERSION] )
    QT+=widgets
}

SRC_DIR = ./src

# Input
HEADERS += $${SRC_DIR}/AtomicInfo.h \
           $${SRC_DIR}/AtomType.h \
           $${SRC_DIR}/Axis.h \
           $${SRC_DIR}/Colours.h \
           $${SRC_DIR}/Component.h \
           $${SRC_DIR}/DensityFrame.h \
           $${SRC_DIR}/ExpData.h \
           $${SRC_DIR}/ExpManager.h \
           $${SRC_DIR}/ExpReader.h \
           $${SRC_DIR}/FileCode.h \
           $${SRC_DIR}/FourierFormFactor.h \
           $${SRC_DIR}/gui.h \
           $${SRC_DIR}/HSVColour.h \
           $${SRC_DIR}/KillBox.h \
           $${SRC_DIR}/LinearSolver.h \
           $${SRC_DIR}/MathUtils.h \
           $${SRC_DIR}/Matrix3D.h \
           $${SRC_DIR}/MDBox.h \
           $${SRC_DIR}/NormalDensity.h \
           $${SRC_DIR}/NormalGrid.h \
	   $${SRC_DIR}/OutputWriter.h \
           $${SRC_DIR}/Particle.h \
           $${SRC_DIR}/ParticleIntensities.h \
           $${SRC_DIR}/ParticleLabel.h \
           $${SRC_DIR}/PolarFourierFormFactor.h \
           $${SRC_DIR}/PolarGrid.h \
           $${SRC_DIR}/Residue.h \
           $${SRC_DIR}/RGB.h \
           $${SRC_DIR}/ScatteringTabFrame.h \
           $${SRC_DIR}/ScatteringType.h \
           $${SRC_DIR}/scrollbar.h \
           $${SRC_DIR}/scrollzoomer.h \
           $${SRC_DIR}/SimManager.h \
           $${SRC_DIR}/SimReader.h \
           $${SRC_DIR}/STEwindow.h \
           $${SRC_DIR}/System.h \
           $${SRC_DIR}/Tab.h \
           $${SRC_DIR}/TabWidgetCloseable.h \
           $${SRC_DIR}/TrajectoryData.h \
           $${SRC_DIR}/TrajectoryInfo.h \
           $${SRC_DIR}/TrajectoryInfoWindow.h \
           $${SRC_DIR}/Util.h \
           $${SRC_DIR}/Vector3D.h \
           $${SRC_DIR}/XrayAtomicFormFactor.h \
           $${SRC_DIR}/XtcFile.h
SOURCES += $${SRC_DIR}/AtomicInfo.cpp \
           $${SRC_DIR}/Colours.cpp \
           $${SRC_DIR}/Component.cpp \
           $${SRC_DIR}/DensityFrame.cpp \
           $${SRC_DIR}/ExpData.cpp \
           $${SRC_DIR}/ExpManager.cpp \
           $${SRC_DIR}/ExpReader.cpp \
           $${SRC_DIR}/FourierFormFactor.cpp \
           $${SRC_DIR}/HSVColour.cpp \
           $${SRC_DIR}/KillBox.cpp \
           $${SRC_DIR}/LinearSolver.cpp \
           $${SRC_DIR}/MathUtils.cpp \
           $${SRC_DIR}/Matrix3D.cpp \
           $${SRC_DIR}/MDBox.cpp \
           $${SRC_DIR}/NormalDensity.cpp \
           $${SRC_DIR}/NormalGrid.cpp \
	   $${SRC_DIR}/OutputWriter.cpp \
           $${SRC_DIR}/Particle.cpp \
           $${SRC_DIR}/ParticleIntensities.cpp \
           $${SRC_DIR}/ParticleLabel.cpp \
           $${SRC_DIR}/PolarFourierFormFactor.cpp \
           $${SRC_DIR}/PolarGrid.cpp \
           $${SRC_DIR}/Residue.cpp \
           $${SRC_DIR}/RGB.cpp \
           $${SRC_DIR}/ScatteringTabFrame.cpp \
           $${SRC_DIR}/scrollbar.cpp \
           $${SRC_DIR}/scrollzoomer.cpp \
           $${SRC_DIR}/SimManager.cpp \
           $${SRC_DIR}/SimReader.cpp \
           $${SRC_DIR}/SIMtoEXP_qt.cpp \
           $${SRC_DIR}/STEwindow.cpp \
           $${SRC_DIR}/System.cpp \
           $${SRC_DIR}/Tab.cpp \
           $${SRC_DIR}/TabWidgetCloseable.cpp \
           $${SRC_DIR}/TrajectoryData.cpp \
           $${SRC_DIR}/TrajectoryInfo.cpp \
           $${SRC_DIR}/TrajectoryInfoWindow.cpp \
           $${SRC_DIR}/Util.cpp \
           $${SRC_DIR}/Vector3D.cpp \
           $${SRC_DIR}/XrayAtomicFormFactor.cpp \
           $${SRC_DIR}/XtcFile.cpp
