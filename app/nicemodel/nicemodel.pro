TEMPLATE = app

TARGET   = nicemodel

CONFIG += qt

QT -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib tbb_include tbb_lib

include(../../orsa.pri)

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaTBB -ltbb -ltbbmalloc
}

macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads
}

win32 {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaSPICE
}

HEADERS = nicemodel.h
SOURCES = nicemodel.cpp
