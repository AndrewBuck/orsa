TEMPLATE = app

TARGET   = SurveySimulator

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
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads /home/buck/Code/orsa/boinc/api/libboinc_api.a /home/buck/Code/orsa/boinc/lib/libboinc.a
}

macx {
        QMAKE_LFLAGS += -static-libgcc
        LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads /Users/buck/Code/orsa/boinc/api/libboinc_api.a /Users/buck/Code/orsa/boinc/lib/libboinc.a
}

win32 {
        LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE
}

HEADERS = parameter.h   telescope.h   SurveySimulator.h   boinc_util.h
SOURCES = parameter.cpp telescope.cpp SurveySimulator.cpp SurveySimulatorMain.cpp
