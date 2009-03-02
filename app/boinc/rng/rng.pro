TEMPLATE = app

TARGET   = rng

CONFIG -= qt

CONFIG += boinc_include boinc_lib gsl_include gsl_lib

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += 
}

macx {
        LIBS += 
}

win32 {
        LIBS += 
}

HEADERS = 
SOURCES = rng.cpp
