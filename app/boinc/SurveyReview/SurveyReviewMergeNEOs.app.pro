TEMPLATE = app

TARGET = SurveyReviewMergeNEOs

CONFIG += qt

QT -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib zlib_include

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
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaUtil -L/home/tricaric/sqlite -lsqlite3
}

HEADERS += grain.h   eta.h   fit.h   skycoverage.h
SOURCES += grain.cpp SurveyReviewMergeNEOsMain.cpp 

HEADERS += sqlite3.h sqlite3ext.h
