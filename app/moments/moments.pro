TEMPLATE = app

CONFIG += qt
QT     -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib tbb_include tbb_lib

include(../../orsa.pri)

TARGET   = moments

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:!macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaEssentialOSG -lorsaSPICE -lOpenThreads
}

macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

win32 {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

HEADERS = 
SOURCES = moments.cpp
