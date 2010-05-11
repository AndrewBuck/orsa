TEMPLATE = app

CONFIG += qt
QT     -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib

include(../../orsa.pri)

TARGET   = mp.$${PLATFORM_NAME}

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

unix:!macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSPICE -lorsaSolarSystem -lorsaEssentialOSG -lOpenThreads 
}

macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

win32 {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

HEADERS = vesta.h
SOURCES = mp.cpp
