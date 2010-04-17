TEMPLATE = app

TARGET = DetectionEfficiencyPlot

CONFIG += qt

QT -= gui

CONFIG += boinc_include boinc_lib gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib zlib_include

win32 {
	CONFIG += osg_lib
}

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

unix:!macx {
	INCLUDEPATH += /home/tricaric/dislin/
	LIBS += -L/home/tricaric/dislin/ -ldislin
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lorsaUtil -lOpenThreads
}

HEADERS = skycoverage.h   eta.h   fit.h
SOURCES = skycoverage.cpp DetectionEfficiencyPlotMain.cpp
