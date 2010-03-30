TEMPLATE = app

TARGET = SurveyReview

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
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaUtil /home/tricaric/boinc/api/libboinc_api.a /home/tricaric/boinc/lib/libboinc.a -L/home/tricaric/sqlite -lsqlite3
}

macx {
        QMAKE_LFLAGS += -static-libgcc
        LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads /Users/tricaric/boinc/api/libboinc_api.a /Users/tricaric/boinc/lib/libboinc.a
}

win32 {
        LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE
}

HEADERS = skycoverage.h   SurveyReview.h   grain.h
SOURCES = skycoverage.cpp SurveyReview.cpp grain.cpp SurveyReviewMain.cpp

HEADERS += sqlite3.h sqlite3ext.h
