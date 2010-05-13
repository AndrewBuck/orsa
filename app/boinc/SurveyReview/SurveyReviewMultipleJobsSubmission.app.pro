TEMPLATE = app

TARGET   = SurveyReviewMultipleJobsSubmission

CONFIG -= qt
#QT -= gui

CONFIG += gmp_include gmp_lib

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

SOURCES += SurveyReviewMultipleJobsSubmission.cpp grain.cpp
