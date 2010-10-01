TEMPLATE = app

TARGET   = SurveyReviewMultipleJobsSubmission

CONFIG += qt

QT -= gui

CONFIG += boinc_include boinc_lib gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib zlib_include

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

SOURCES += SurveyReviewMultipleJobsSubmission.cpp SurveyReview.cpp grain.cpp

unix:!macx {
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaUtil /home/buck/Code/orsa/boinc/api/libboinc_api.a /home/buck/Code/orsa/boinc/lib/libboinc.a -L/home/buck/Code/orsa/sqlite -lsqlite3
}
